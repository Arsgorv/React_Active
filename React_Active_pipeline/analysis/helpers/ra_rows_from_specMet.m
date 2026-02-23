function rows = ra_rows_from_specMet(sessname, tag, Met, bandNames)
% ra_rows_from_specMet
% Convert Met (from eventlocked_spec_grid_plot) into long rows:
% session, tag, modality='SPEC', group, band, window, metric, value

rows = {};
if nargin < 4, bandNames = {'delta','theta','gamma','highgamma'}; end
if ~isfield(Met,'group') || isempty(Met.group), return, end
if ~isfield(Met,'metricBands') || isempty(Met.metricBands), return, end

for g = 1:numel(Met.group)
    if ~isfield(Met.group(g),'name'), continue, end
    gname = string(Met.group(g).name);

    for bb = 1:numel(Met.group(g).band)
        if ~isfield(Met.group(g).band(bb),'range'), continue, end
        br = Met.group(g).band(bb).range;
        bidx = find_band_idx_tol(br, Met.metricBands);
        if isempty(bidx)
            bname = '';
        else
            bname = bandNames{min(bidx, numel(bandNames))};
        end

        if ~isfield(Met.group(g).band(bb),'win'), continue, end
        for ww = 1:numel(Met.group(g).band(bb).win)
            W = Met.group(g).band(bb).win(ww);
            if ~isfield(W,'name'), continue, end
            wname = string(W.name);

            if isfield(W,'mean')
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, wname, 'mean', W.mean}; %#ok<AGROW>
            end
            if isfield(W,'sem')
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, wname, 'sem', W.sem}; %#ok<AGROW>
            end
            if isfield(W,'slope_mean')
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, wname, 'slope_mean', W.slope_mean}; %#ok<AGROW>
            end
            if isfield(W,'slope_sem')
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, wname, 'slope_sem', W.slope_sem}; %#ok<AGROW>
            end
            if isfield(W,'n')
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, wname, 'n', W.n}; %#ok<AGROW>
            end
        end

        % optional post-pre summary
        if isfield(Met.group(g).band(bb),'postpre_mean')
            rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, 'postpre', 'mean', Met.group(g).band(bb).postpre_mean}; %#ok<AGROW>
        end
        if isfield(Met.group(g).band(bb),'postpre_sem')
            rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, 'postpre', 'sem', Met.group(g).band(bb).postpre_sem}; %#ok<AGROW>
        end
        if isfield(Met.group(g).band(bb),'postpre_n')
            rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, 'postpre', 'n', Met.group(g).band(bb).postpre_n}; %#ok<AGROW>
        end
    end
end

end

function idx = find_band_idx_tol(br, metricBands)
idx = [];
if isempty(br) || numel(br) ~= 2, return, end
br = double(br(:))';
tol = 1e-9;
for k = 1:size(metricBands,1)
    mb = double(metricBands(k,:));
    if all(abs(mb - br) < tol)
        idx = k;
        return
    end
end
end
