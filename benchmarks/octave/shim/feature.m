function out = feature(name, varargin)
% Octave shim for MATLAB's feature('numCores') used by the CPU scheduler block.
switch lower(name)
  case 'numcores'
    out = nproc();
  otherwise
    error('feature(%s) is not available in this shim', name);
end
