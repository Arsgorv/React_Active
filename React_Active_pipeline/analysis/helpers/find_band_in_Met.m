function bb = find_band_in_Met(Met, gi, br)
% find_band_in_Met
% Find band index bb in Met.group(gi).band whose .range matches br.
% Uses tolerance to avoid float/shape issues.

bb = [];
if nargin < 3 || isempty(Met) || isempty(gi) || isempty(br), return, end
if ~isfield(Met,'group') || gi > numel(Met.group), return, end
if ~isfield(Met.group(gi),'band'), return, end

br = double(br(:))';
tol = 1e-9;

for k = 1:numel(Met.group(gi).band)
    if ~isfield(Met.group(gi).band(k),'range'), continue, end
    r = Met.group(gi).band(k).range;
    if isempty(r), continue, end
    r = double(r(:))';
    if numel(r) ~= 2, continue, end
    if all(abs(r - br) < tol)
        bb = k;
        return
    end
end
end
