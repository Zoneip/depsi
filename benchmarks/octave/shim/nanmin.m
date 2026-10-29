function [m, idx] = nanmin(x, varargin)
x(~isfinite(x)) = Inf;
[m, idx] = min(x, varargin{:});
