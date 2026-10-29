function s = nanstd(x, flag, dim)
if nargin < 2 || isempty(flag), flag = 0; end
if nargin < 3 || isempty(dim), dim = find(size(x) > 1, 1); end
if nargin < 3 && isempty(dim), dim = 1; end
ok = ~isnan(x);
n = sum(ok, dim);
mu = sum(x .* ok, dim) ./ n;
d = (x - mu) .^ 2 .* ok;
if flag
  s = sqrt(sum(d, dim) ./ n);
else
  s = sqrt(sum(d, dim) ./ max(n - 1, 1));
end
