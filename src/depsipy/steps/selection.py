"""PSC / PSP selection (port of ``ps_selection.m``).

Two passes over the data:

1. **Amplitude pass** -- stream the SLCs in row blocks and accumulate
   ``sum(amp)`` and ``sum(amp^2)`` over the selected scenes.  This yields the
   multi-reflectivity map (mean amplitude), the amplitude dispersion
   ``std/mean`` and the PSC/PSP candidates.  Streaming keeps the memory
   footprint at one block per scene instead of the whole stack.
2. **Phase pass** -- read the phase and the height-to-phase factor *only at the
   candidate pixels*.  The MATLAB version loads complete interferogram buffers
   for every candidate grid cell, which is the dominant cost of
   ``ps_selection.m`` for large crops; here the number of values read scales
   with the number of candidates, not with the raster size.

Side-lobe suppression follows the modified (``_hhf``) branch of the toolbox: a
3x3 local maximum above the 95th amplitude percentile is a candidate, and
candidates that are coherent with a stronger neighbour in the same row/column
are suppressed (they are the side lobes of a point target).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import ndimage

from ..io.raw import RawFile, read_pixels, write_raw

__all__ = [
    "PSCCandidate",
    "SelectionResult",
    "grid_arrays",
    "sidelobe_mask",
    "select_psc_psp",
    "write_psc_files",
]


# ---------------------------------------------------------------------------
# grid definition (ps_selection.m, "Determine borders gridcells")
# ---------------------------------------------------------------------------


def grid_arrays(Nlines: int, Npixels: int, az_spacing: float, r_spacing: float, gridsize: float):
    """Return ``(grid_az, grid_r)``: (Ng, 2) arrays of 1-based cell borders.

    The first cell absorbs the remainder so that the last cell always ends on the
    last line/pixel, exactly as in the MATLAB code.
    """
    Nlines_g = int(round(gridsize / az_spacing))
    Npixels_g = int(round(gridsize / r_spacing))
    if Nlines_g < 1 or Npixels_g < 1:
        raise ValueError("psc_selection_gridsize is smaller than one pixel")

    def _edges(total, size):
        """Cell borders (1-based, inclusive) exactly as ps_selection.m builds them."""
        Ng = total // size
        offset = (total % size) // 2
        if Ng == 0:  # raster smaller than one cell: a single cell
            return np.array([[1, total]], dtype=int)
        stops = offset + np.arange(1, Ng + 1) * size
        stops[-1] = total
        starts = np.empty(Ng, dtype=int)
        starts[0] = 1
        starts[1:] = stops[:-1] + 1
        return np.column_stack([starts, stops])

    return _edges(Nlines, Nlines_g), _edges(Npixels, Npixels_g)


# ---------------------------------------------------------------------------
# side-lobe mask
# ---------------------------------------------------------------------------


def _dependent(theta, mask, threshold=0.95):
    """Coherence of each candidate with the reference, over the scenes."""
    # |mean(exp(i(theta_cand - theta_ref)))|, Nifgs columns
    d = np.exp(1j * (theta - theta[..., :1] if theta.ndim > 1 else theta - theta[0]))
    if d.ndim == 1:
        return 1.0
    return np.abs(d.mean(axis=-1))


def sidelobe_mask(
    mean_amp: np.ndarray,
    theta_at=None,
    *,
    peak_tolerance: float = 0.70,
    strong_percentile: float = 95.0,
    coherence_threshold: float = 0.95,
    enabled: bool = True,
) -> np.ndarray:
    """Boolean map of *suppressed* pixels (``True`` = side lobe).

    Parameters
    ----------
    mean_amp
        Multi-reflectivity map.
    theta_at
        Callable ``theta_at(az, r) -> (Ncand, Nifgs)`` returning the SLC phases
        of the given (0-based) pixels, used for the coherence de-clustering.
        When ``None`` only the amplitude criterion is used.
    """
    mask = np.zeros(mean_amp.shape, dtype=bool)
    if not enabled:
        return mask
    amp = np.nan_to_num(mean_amp, nan=0.0)
    finite = amp > 0.0
    if not finite.any():
        return mask

    maxima = ndimage.maximum_filter(amp, size=3)
    local_max = (maxima * peak_tolerance <= amp) & finite
    strong_thr = np.percentile(amp[finite], strong_percentile)
    strong = local_max & (amp >= strong_thr)
    if not strong.any():
        return mask

    r, c = np.nonzero(strong)  # r = row index (0-based)
    kept = np.arange(r.size)
    if theta_at is not None and r.size > 1:
        theta = theta_at(r, c)
        kept = _decluster(r, c, theta, amp[r, c], coherence_threshold)
    # everything that is a strong local maximum but did not survive is masked
    suppressed = np.ones(r.size, dtype=bool)
    suppressed[kept] = False
    mask[r[suppressed], c[suppressed]] = True
    return mask


def _decluster(rows, cols, theta, cand_amp, threshold):
    """Suppress strong candidates that are coherent with a stronger neighbour.

    Port of the ``GexYuu`` block of ``ps_selection.m``: a range sweep (per row),
    an azimuth sweep (per column) and a nearest-neighbour pass.  Each sweep takes
    the **strongest remaining** candidate as the reference and removes the
    contiguous run of coherent candidates around it; the nearest-neighbour pass
    drops the weaker member of every coherent pair (simultaneously per pass, as
    the reference does).

    ``theta`` holds the interferometric phases ``(Ncand, Nifgs)``; ``threshold`` is
    the coherence threshold (0.95 in the reference).
    """
    cand_amp = np.asarray(cand_amp, dtype=float)
    theta_c = np.exp(1j * np.asarray(theta, dtype=float))
    n_ifg = theta_c.shape[1]
    if n_ifg == 0:
        return np.arange(rows.size)

    def sweep(kept, group_coord, sort_coord):
        """One directional sweep over the candidates still alive."""
        groups = group_coord[kept]
        idx_alive = kept
        out = []
        for gv in np.unique(groups):
            members = np.nonzero(groups == gv)[0]
            if members.size == 1:
                out.append(idx_alive[members[0]])
                continue
            order = members[np.argsort(sort_coord[idx_alive[members]], kind="stable")]
            idx = idx_alive[order]
            amp_eff = cand_amp[idx].astype(float).copy()
            ph = theta_c[idx]
            coh = np.abs(ph @ np.conj(ph).T) / n_ifg
            removed = np.zeros(idx.size, dtype=bool)
            while True:
                m = int(np.argmax(amp_eff))
                out.append(idx[m])  # the reference max is always kept
                row = coh[m]
                not_coh = (row <= threshold) | removed
                pos = np.nonzero(not_coh[m + 1 :])[0]
                max_i = m + int(pos[0]) if pos.size else idx.size - 1
                neg = np.nonzero(not_coh[:m])[0]
                min_i = int(neg[-1]) + 1 if neg.size else 0
                removed[min_i : max_i + 1] = True
                amp_eff[min_i : max_i + 1] = -np.inf
                if removed.all():
                    break
        return np.asarray(sorted(out), dtype=np.int64)

    kept = np.arange(rows.size, dtype=np.int64)
    kept = sweep(kept, rows, cols)  # range direction: fixed row, order by column
    kept = sweep(kept, cols, rows)  # azimuth direction: fixed column, order by row
    if kept.size < 2:
        return kept

    # nearest-neighbour pass: drop the weaker member of every coherent pair
    alive = np.zeros(rows.size, dtype=bool)
    alive[kept] = True
    alive_idx = kept.copy()
    while alive_idx.size > 1:
        near = _nearest_positive(alive_idx, alive_idx, rows, cols)
        has = near >= 0
        if not has.any():
            break
        ii, jj = alive_idx[has], near[has]
        coh = np.abs(np.einsum("ij,ij->i", theta_c[ii], np.conj(theta_c[jj]))) / n_ifg
        hit = coh > threshold
        if not hit.any():
            break
        weaker = np.where(cand_amp[ii] < cand_amp[jj], ii, jj)
        drop = np.unique(weaker[hit])
        alive[drop] = False
        alive_idx = np.nonzero(alive)[0]
    return alive_idx


def _nearest_positive(query, pool, rows, cols, *, block: int = 512):
    """Nearest member of ``pool`` with strictly positive distance, per query.

    Returns absolute candidate indices, ``-1`` where no such neighbour exists.
    Ties resolve to the lowest position in ``pool`` (matching ``np.argmin`` over
    the original ordered scan); entries at identical coordinates are excluded,
    exactly like the reference loop's ``d2 > 0`` mask.

    Large pools go through a k-d tree; the blocked quadratic path is kept as the
    exact fallback (and is used directly for small pools).
    """
    query = np.asarray(query, dtype=np.int64)
    pool = np.asarray(pool, dtype=np.int64)
    if pool.size > 512 and query.size > 32:
        tree_based = _nearest_positive_kdtree(query, pool, rows, cols)
        if tree_based is not None:
            return tree_based
    return _nearest_positive_blocked(query, pool, rows, cols, block=block)


def _nearest_positive_kdtree(query, pool, rows, cols):
    """k-d tree path; returns ``None`` when ties cannot be resolved safely.

    Distances are sorted by the tree, so the first strictly positive neighbour is
    the nearest one; only rows whose nearest distance is attained by several pool
    entries need the exact tie-break (lowest pool position, as ``argmin`` over the
    reference order gives).  Duplicate coordinates at distance zero are skipped.
    """
    from scipy.spatial import cKDTree

    pts = np.column_stack([cols[pool].astype(np.float64), rows[pool].astype(np.float64)])
    tree = cKDTree(pts)
    k = min(16, pool.size)
    dist, idx = tree.query(pts[query] if query is not pool else pts, k=k)
    dist = np.atleast_2d(dist)
    idx = np.atleast_2d(idx)
    n_rows = dist.shape[0]
    positive = dist > 0
    has = positive.any(axis=1)
    if not has.all():
        return None  # all k neighbours coincide with the query point
    first = np.argmax(positive, axis=1)
    d_min = dist[np.arange(n_rows), first]
    tol = d_min * (1 + 16 * np.finfo(float).eps)
    tied = positive & (dist <= tol[:, None])
    out = pool[idx[np.arange(n_rows), first]].astype(np.int64)
    ambiguous = np.nonzero(tied.sum(axis=1) > 1)[0]
    for row in ambiguous:
        cand = [int(pos) for pos in tree.query_ball_point(pts[query[row]], tol[row])
                if dist2(pts[query[row]], pts[pos]) > 0]
        if not cand:  # pragma: no cover - defensive
            return None
        out[row] = pool[min(cand)]
    return out


def dist2(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return dx * dx + dy * dy


def _nearest_positive_blocked(query, pool, rows, cols, *, block: int = 512):
    """Exact quadratic path with the reference tie-breaking rule."""
    query = np.asarray(query, dtype=np.int64)
    pool = np.asarray(pool, dtype=np.int64)
    out = np.full(query.size, -1, dtype=np.int64)
    if pool.size == 0:
        return out
    rq = rows[query].astype(np.float64)
    cq = cols[query].astype(np.float64)
    rp = rows[pool].astype(np.float64)
    cp = cols[pool].astype(np.float64)
    for start in range(0, query.size, block):
        stop = min(start + block, query.size)
        d2 = (cq[start:stop, None] - cp[None, :]) ** 2 + (rq[start:stop, None] - rp[None, :]) ** 2
        finite = np.isfinite(d2) & (d2 > 0)
        d2 = np.where(finite, d2, np.inf)
        best = np.argmin(d2, axis=1)
        ok = finite[np.arange(stop - start), best]
        out[start:stop] = np.where(ok, pool[best], -1)
    return out


# ---------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------


@dataclass
class PSCCandidate:
    """One selected persistent scatterer candidate."""

    grid_az: int
    grid_r: int
    az: int  # 1-based
    r: int  # 1-based
    phase: np.ndarray  # (Nifgs,) interferometric phase [rad]
    h2ph: np.ndarray  # (Nifgs,) height-to-phase factors [1/m]
    amp: np.ndarray  # (Nslc,) amplitudes
    amp_disp: float


@dataclass
class SelectionResult:
    grid_az: np.ndarray
    grid_r: np.ndarray
    Npsc: np.ndarray  # per selection
    Npsp: np.ndarray
    psc: list[list[PSCCandidate]] = field(default_factory=list)
    psp_az: list[np.ndarray] = field(default_factory=list)
    psp_r: list[np.ndarray] = field(default_factory=list)
    psp_amp_disp: list[np.ndarray] = field(default_factory=list)
    psp_phase: list[np.ndarray] = field(default_factory=list)
    psp_h2ph: list[np.ndarray] = field(default_factory=list)
    psp_amp: list[np.ndarray] = field(default_factory=list)
    mean_amp: np.ndarray | None = None
    amp_disp: np.ndarray | None = None
    side_lobe: np.ndarray | None = None
    sigma_amp: np.ndarray | None = None

    def ps_coords(self, z: int = 0):
        cand = self.psc[z]
        az = np.array([c.az for c in cand])
        r = np.array([c.r for c in cand])
        return az, r


def select_psc_psp(
    *,
    slc_files: list[str],
    ifg_files: list[str],
    h2ph_files: list[str],
    Nlines: int,
    Npixels: int,
    az_spacing: float,
    r_spacing: float,
    crop,
    crop_in,
    slc_selection,
    calfactors,
    gridsize: float,
    psc_threshold: float,
    psp_threshold1: float,
    psp_threshold2: float | None,
    psc_distribution: str = "nonuniform",
    do_apriori_sidelobe_mask: str = "yes",
    livetime_threshold: float = 0.2,
    peak_tolerance: float = 0.99,
    Npsc_selections: int = 1,
    ref_cn=None,
    workspace: str | os.PathLike | None = None,
    block_rows: int | None = None,
    progress=None,
) -> SelectionResult:
    """Select PSC (and PSP) candidates from a stack of SLC/interferogram rasters.

    Parameters are the ones of ``ps_selection.m``; ``crop``/``crop_in`` are
    ``(az0, azN, r0, rN)`` in 1-based global coordinates.
    """
    ifg_files = [os.fspath(f) for f in ifg_files]
    slc_files = [os.fspath(f) for f in slc_files]
    Nifgs = len(ifg_files)
    Nslc = len(slc_files)
    sel = np.asarray(slc_selection, dtype=int)
    cal = np.asarray(calfactors, dtype=float)

    crop = tuple(int(v) for v in crop)
    crop_in = np.asarray(crop_in, dtype=int).reshape(-1, 4)

    # block size: enough rows to amortise the per-block overhead, small enough
    # that one block of every scene fits comfortably in memory
    if block_rows is None:
        bytes_per_row = 8 * (2 * Nslc + 2 * Nifgs + 1)  # complex64 + accumulators
        block_rows = max(1, int(min(512, 2e8 / max(bytes_per_row * Npixels, 1))))
    block_rows = int(min(block_rows, Nlines))

    grid_az, grid_r = grid_arrays(Nlines, Npixels, az_spacing, r_spacing, gridsize)

    sum_amp = np.zeros((Nlines, Npixels), dtype=np.float64)
    sum_amp2 = np.zeros((Nlines, Npixels), dtype=np.float64)

    # ---- pass 1: amplitudes ------------------------------------------------
    for b0, b1 in _blocks(Nlines, block_rows):
        acc1 = np.zeros((b1 - b0 + 1, Npixels), dtype=np.float64)
        acc2 = np.zeros_like(acc1)
        for k in range(sel.size):
            cin = crop_in[int(sel[k]) - 1]  # MATLAB indices are 1-based
            nlines_file = int(cin[1] - cin[0] + 1)
            loff = crop[0] - int(cin[0])
            poff = crop[2] - int(cin[2])
            block = _read_window(
                slc_files[int(sel[k]) - 1], nlines_file, "cpxfloat32",
                b0 + loff, b1 + loff, 1 + poff, Npixels + poff,
            )
            a = np.abs(block) / cal[k]
            acc1 += a
            acc2 += a * a
        n = sel.size
        mean_amp = acc1 / n
        var = np.maximum(acc2 / n - mean_amp**2, 0.0) * (n / max(n - 1, 1))
        std_amp = np.sqrt(var)
        with np.errstate(divide="ignore", invalid="ignore"):
            amp_disp = np.where(mean_amp > 0, std_amp / mean_amp, np.inf)
        sum_amp[b0 - 1 : b1] = mean_amp
        sum_amp2[b0 - 1 : b1] = amp_disp
        if progress:
            progress("selection/amplitude", b1 / Nlines)

    mrm = np.log10(np.maximum(sum_amp, 1e-6)).astype(np.float32)
    amp_disp_map = sum_amp2
    del sum_amp2

    # ---- side-lobe mask ----------------------------------------------------
    def theta_at(rows, cols):
        """Phases of the candidate pixels, shape (Ncand, Nifgs) [rad]."""
        out = np.empty((rows.size, Nifgs), dtype=float)
        for i in range(Nifgs):
            cin = crop_in[i]
            nlines_file = int(cin[1] - cin[0] + 1)
            loff = crop[0] - int(cin[0])
            poff = crop[2] - int(cin[2])
            out[:, i] = np.angle(
                read_pixels(
                    ifg_files[i], nlines_file, "cpxfloat32", rows + 1 + loff, cols + 1 + poff
                )
            )
        return out

    sl = sidelobe_mask(
        sum_amp,
        theta_at,
        peak_tolerance=peak_tolerance,
        enabled=str(do_apriori_sidelobe_mask).lower() in ("yes", "piers"),
    )
    amp_disp_sel = amp_disp_map.copy()
    amp_disp_sel[sl] = np.nan
    amp_disp_sel[~np.isfinite(amp_disp_sel)] = np.nan

    # ---- pass 2: grid selection -------------------------------------------
    psc_sel = _grid_select(
        amp_disp_sel,
        grid_az,
        grid_r,
        Npsc_selections,
        psc_threshold,
        ref_cn=ref_cn,
    )
    if psc_distribution == "uniform":
        psc_sel = _uniform_distribution(
            psc_sel, amp_disp_sel, grid_az, grid_r, Npsc_selections, psc_threshold, az_spacing, r_spacing, ref_cn
        )
    elif psc_distribution != "nonuniform":
        raise ValueError("psc_distribution must be 'uniform' or 'nonuniform'")

    # PSP candidates: amplitude dispersion below the first threshold
    psp_mask = np.isfinite(amp_disp_sel) & (amp_disp_sel < float(psp_threshold1))
    psp_r_idx, psp_c_idx = np.nonzero(psp_mask)
    psp_disp = amp_disp_sel[psp_r_idx, psp_c_idx]

    # ---- read the observations at the PSP pixels --------------------------
    Npsp = psp_r_idx.size
    psp_phase = np.empty((Npsp, Nifgs))
    psp_h2ph = np.empty((Npsp, Nifgs))
    psp_amp = np.empty((Npsp, Nslc))
    for i in range(Nifgs):
        cin = crop_in[i]
        nlines_file = int(cin[1] - cin[0] + 1)
        loff = crop[0] - int(cin[0])
        poff = crop[2] - int(cin[2])
        psp_phase[:, i] = np.angle(
            read_pixels(ifg_files[i], nlines_file, "cpxfloat32",
                        psp_r_idx + 1 + loff, psp_c_idx + 1 + poff)
        )
        psp_h2ph[:, i] = read_pixels(
            h2ph_files[i], nlines_file, "float32",
            psp_r_idx + 1 + loff, psp_c_idx + 1 + poff,
        )
    for k in range(sel.size):
        cin = crop_in[int(sel[k]) - 1]
        nlines_file = int(cin[1] - cin[0] + 1)
        loff = crop[0] - int(cin[0])
        poff = crop[2] - int(cin[2])
        psp_amp[:, k] = np.abs(
            read_pixels(slc_files[int(sel[k]) - 1], nlines_file, "cpxfloat32",
                        psp_r_idx + 1 + loff, psp_c_idx + 1 + poff)
        )

    # ---- read phases at the selected pixels -------------------------------
    psc_lists: list[list[PSCCandidate]] = []
    for z in range(Npsc_selections):
        pts = psc_sel[z]  # list of (row, col) 0-based
        if not pts:
            psc_lists.append([])
            continue
        rows = np.array([p[0] for p in pts])
        cols = np.array([p[1] for p in pts])
        phase = np.empty((rows.size, Nifgs))
        h2ph = np.empty((rows.size, Nifgs))
        for i in range(Nifgs):
            cin = crop_in[i]
            nlines_file = int(cin[1] - cin[0] + 1)
            loff = crop[0] - int(cin[0])
            poff = crop[2] - int(cin[2])
            phase[:, i] = np.angle(
                read_pixels(ifg_files[i], nlines_file, "cpxfloat32", rows + 1 + loff, cols + 1 + poff)
            )
            h2ph[:, i] = read_pixels(h2ph_files[i], nlines_file, "float32", rows + 1 + loff, cols + 1 + poff)
        amps = np.empty((rows.size, Nslc))
        for k in range(sel.size):
            cin = crop_in[int(sel[k]) - 1]
            nlines_file = int(cin[1] - cin[0] + 1)
            loff = crop[0] - int(cin[0])
            poff = crop[2] - int(cin[2])
            amps[:, k] = np.abs(
                read_pixels(
                    slc_files[int(sel[k]) - 1], nlines_file, "cpxfloat32",
                    rows + 1 + loff, cols + 1 + poff,
                )
            )
        disp = amp_disp_sel[rows, cols]
        order = np.lexsort((cols, rows))
        cand = [
            PSCCandidate(
                grid_az=int(grid_az[0, 0]),
                grid_r=int(grid_r[0, 0]),
                az=int(rows[k]) + 1,
                r=int(cols[k]) + 1,
                phase=phase[k],
                h2ph=h2ph[k],
                amp=amps[k],
                amp_disp=float(disp[k]),
            )
            for k in order
        ]
        psc_lists.append(cand)

    return SelectionResult(
        grid_az=grid_az,
        grid_r=grid_r,
        Npsc=np.array([len(c) for c in psc_lists], dtype=int),
        Npsp=np.array([psp_r_idx.size] * Npsc_selections, dtype=int),
        psc=psc_lists,
        psp_az=[psp_r_idx + 1] * Npsc_selections,
        psp_r=[psp_c_idx + 1] * Npsc_selections,
        psp_amp_disp=[psp_disp] * Npsc_selections,
        psp_phase=[psp_phase] * Npsc_selections,
        psp_h2ph=[psp_h2ph] * Npsc_selections,
        psp_amp=[psp_amp] * Npsc_selections,
        mean_amp=sum_amp,
        amp_disp=amp_disp_map,
        side_lobe=sl,
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _blocks(total, block):
    start = 0
    while start < total:
        stop = min(start + block, total)
        yield start + 1, stop
        start = stop


def _read_window(path, lines_file, fmt, r0, rN, c0, cN):
    with RawFile(path, lines_file, fmt) as f:
        return f.read(r0, rN, c0, cN)


def _grid_select(amp_disp, grid_az, grid_r, Npsc_selections, threshold, ref_cn=None):
    """Per grid cell, pick the ``Npsc_selections`` minima below ``threshold``."""
    out = [[] for _ in range(Npsc_selections)]
    for a0, a1 in grid_az:
        for r0, r1 in grid_r:
            cell = amp_disp[a0 - 1 : a1, r0 - 1 : r1]
            flat = cell.ravel()
            ok = np.isfinite(flat) & (flat < threshold)
            if not ok.any():
                continue
            idx = np.nonzero(ok)[0]
            vals = flat[idx]
            for z, zi in enumerate(np.argsort(vals)[:Npsc_selections]):
                p = int(idx[zi])
                out[z].append((a0 - 1 + p // cell.shape[1], r0 - 1 + p % cell.shape[1]))

    # the reference point is always kept, even if its cell is above threshold
    if ref_cn is not None and np.all(np.isfinite(ref_cn)):
        ra, rr = int(ref_cn[0]) - 1, int(ref_cn[1]) - 1
        if 0 <= ra < amp_disp.shape[0] and 0 <= rr < amp_disp.shape[1]:
            for z in range(Npsc_selections):
                if (ra, rr) not in out[z]:
                    out[z].append((ra, rr))
    return out


def _uniform_distribution(pts, amp_disp, grid_az, grid_r, Npsc_selections, threshold,
                          az_spacing, r_spacing, ref_cn):
    """Keep one PSC per *shifted* grid cell (``psc_distribution='uniform'``).

    The grid is shifted by half a cell in range and in azimuth; applying the
    keep-best rule on both shifted grids leaves at most one PSC per
    unshifted cell while removing the dependence on where the cell borders fell.
    """
    def shifted(edges, total):
        size = int(edges[0, 1] - edges[0, 0] + 1)
        half = size // 2
        out = [(int(a) - half, int(b) - half) for a, b in edges]
        out[0] = (int(edges[0, 0]), int(edges[0, 1]) - half)
        out.append((total - half + 1, total))
        return out

    current = [list(p) for p in pts]
    grids = (
        (shifted(grid_az, amp_disp.shape[0]), [(int(a), int(b)) for a, b in grid_r]),
        ([(int(a), int(b)) for a, b in grid_az], shifted(grid_r, amp_disp.shape[1])),
    )
    for edges_a, edges_r in grids:
        for z in range(Npsc_selections):
            kept, used = [], set()
            for a0, a1 in edges_a:
                for r0, r1 in edges_r:
                    best = None
                    for i, (ra, rr) in enumerate(current[z]):
                        if i in used:
                            continue
                        if a0 <= ra + 1 <= a1 and r0 <= rr + 1 <= r1:
                            val = amp_disp[ra, rr]
                            if best is None or val < best[0]:
                                best = (val, i)
                    if best is not None:
                        used.add(best[1])
                        kept.append(current[z][best[1]])
            current[z] = kept

    if ref_cn is not None and np.all(np.isfinite(ref_cn)):
        ra, rr = int(ref_cn[0]) - 1, int(ref_cn[1]) - 1
        for z in range(Npsc_selections):
            if (ra, rr) not in current[z]:
                current[z].append((ra, rr))
    return current


def write_psc_files(result: SelectionResult, workspace, project_id: str, Nifgs: int = 0):
    """Write ``_psc_2orig_sel*.raw`` / ``_psp_2orig_sel*.raw`` (DePSI layout)."""
    workspace = Path(workspace)
    for z, cand in enumerate(result.psc):
        n = len(cand)
        Nifgs = cand[0].phase.size if n else Nifgs
        data = np.empty((n, 2 * Nifgs + 7))
        for k, c in enumerate(cand):
            data[k, 0] = k + 1
            data[k, 1] = 1
            data[k, 2] = c.grid_az
            data[k, 3] = c.grid_r
            data[k, 4] = c.az
            data[k, 5] = c.r
            data[k, 6 : 6 + Nifgs] = c.phase
            data[k, 6 + Nifgs : 6 + 2 * Nifgs] = c.h2ph
            data[k, 6 + 2 * Nifgs] = c.amp_disp
        write_raw(workspace / f"{project_id}_psc_2orig_sel{z + 1}.raw", data, "float64")
        amp = np.array([c.amp for c in cand]) if n else np.zeros((0, 0))
        write_raw(workspace / f"{project_id}_psc_validation_sel{z + 1}.raw", amp, "float64")

    for z in range(len(result.psp_az)):
        az = np.asarray(result.psp_az[z])
        r = np.asarray(result.psp_r[z])
        n = az.size
        Nifgs = result.psc[z][0].phase.size if result.psc[z] else Nifgs
        phase = (
            np.asarray(result.psp_phase[z])
            if z < len(result.psp_phase) and result.psp_phase[z] is not None
            else np.zeros((n, Nifgs))
        )
        h2ph = (
            np.asarray(result.psp_h2ph[z])
            if z < len(result.psp_h2ph) and result.psp_h2ph[z] is not None
            else np.zeros((n, Nifgs))
        )
        if phase.shape != (n, Nifgs) or h2ph.shape != (n, Nifgs):
            raise ValueError("PSP phase/h2ph shape does not match az/r")
        cell_a = np.searchsorted(result.grid_az[:, 0], az, side="right") - 1
        cell_r = np.searchsorted(result.grid_r[:, 0], r, side="right") - 1
        cell_a = np.clip(cell_a, 0, result.grid_az.shape[0] - 1)
        cell_r = np.clip(cell_r, 0, result.grid_r.shape[0] - 1)
        data = np.empty((n, 2 * Nifgs + 4))
        data[:, 0] = cell_a
        data[:, 1] = cell_r
        data[:, 2] = az
        data[:, 3] = r
        data[:, 4 : 4 + Nifgs] = phase
        data[:, 4 + Nifgs :] = h2ph
        write_raw(workspace / f"{project_id}_psp_2orig_sel{z + 1}.raw", data, "float64")

        amp_disp = np.asarray(result.psp_amp_disp[z]).ravel()
        amp = (
            np.asarray(result.psp_amp[z])
            if z < len(result.psp_amp) and result.psp_amp[z] is not None
            else np.zeros((n, 0))
        )
        validation = np.column_stack([amp_disp, amp]) if amp.size else amp_disp[:, None]
        write_raw(workspace / f"{project_id}_psp_validation_sel{z + 1}.raw", validation, "float64")
