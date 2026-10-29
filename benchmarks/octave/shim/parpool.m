function pool = parpool(varargin)
% Octave shim: no parallel pool, ``parfor`` runs serially anyway.
pool = 1;
