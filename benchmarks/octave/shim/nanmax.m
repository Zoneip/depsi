function [m, idx] = nanmax(x, varargin)
x(~isfinite(x)) = -Inf;
[m, idx] = max(x, varargin{:});
