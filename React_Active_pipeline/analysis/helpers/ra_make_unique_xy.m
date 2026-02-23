function [t_u, y_u] = ra_make_unique_xy(t_s, y)
% ra_make_unique_xy
% Make time vector unique and strictly increasing for interp1.
% Robust to:
%   - row/column inputs
%   - y being NxM
%   - length mismatch between t_s and y (cropped to min length)
% Duplicate timestamps are collapsed by averaging y across duplicates.

t_s = double(t_s(:));

y = double(y);
if isvector(y)
    y = y(:);
end

% If lengths mismatch, crop to common length
n = min(numel(t_s), size(y,1));
t_s = t_s(1:n);
y   = y(1:n,:);

% Remove NaNs/infs using a mask on t and all y cols
ok = isfinite(t_s);
if ~isempty(y)
    ok = ok & all(isfinite(y), 2);
end

t_s = t_s(ok);
y   = y(ok,:);

if numel(t_s) < 2
    t_u = t_s;
    y_u = y;
    return
end

% Sort
[t_s, ord] = sort(t_s);
y = y(ord,:);

% Unique time points and group indices
[t_u, ~, ic] = unique(t_s, 'sorted');

% Average duplicates per column
y_u = nan(numel(t_u), size(y,2));
for c = 1:size(y,2)
    y_u(:,c) = accumarray(ic, y(:,c), [], @mean);
end

% Final finite guard
ok = isfinite(t_u) & all(isfinite(y_u),2);
t_u = t_u(ok);
y_u = y_u(ok,:);

end
