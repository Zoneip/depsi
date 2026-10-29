function m = nanmean(x, dim)
if nargin < 2 || isempty(dim), dim = find(size(x) > 1, 1); end
if nargin < 2 && isempty(dim), dim = 1; end
m = sum(x .* ~isnan(x), dim) ./ sum(~isnan(x), dim);
