function BRES = RAA_collect_behaviour_marker_AUROC_across_sessions(sessions, opts)
% Collect per-session AUROC for one behaviour marker from behaviour_analysis outputs.
%
% Expected files (per session):
%   <datapath>/analysis/behaviour/<bodypart>/*_metrics.csv
% We pick rows matching opts.subset / opts.window / opts.marker.
%
% Outputs:
%   BRES.T : table with session, AUROC, nTar, nRef, marker, subset, window

if nargin < 2, opts = struct(); end
if ~isfield(opts,'bodypart'), opts.bodypart = 'Pupil-EyeCam'; end
if ~isfield(opts,'marker'),   opts.marker   = 'pupil_area_007'; end
if ~isfield(opts,'subset'),   opts.subset   = 'regular'; end
if ~isfield(opts,'window'),   opts.window   = 'stimoff_to_arrival'; end

rows = {};
for s = 1:numel(sessions)
    dp = sessions{s};
    [~,sessname] = fileparts(dp);

    behDir = fullfile(dp,'analysis','behaviour', opts.bodypart);
    if ~exist(behDir,'dir'), continue, end

    % behaviour_analysis writes multiple condition-specific metrics.csv files.
    F = dir(fullfile(behDir, '*_metrics.csv'));
    if isempty(F), continue, end

    % scan all csvs in this bodypart folder, take the first matching row
    found = false;
    for fi = 1:numel(F)
        T = readtable(fullfile(F(fi).folder, F(fi).name));
        if isempty(T), continue, end

        % normalize strings
        if ismember('marker',T.Properties.VariableNames), T.marker = string(T.marker); end
        if ismember('subset',T.Properties.VariableNames), T.subset = string(T.subset); end
        if ismember('window',T.Properties.VariableNames), T.window = string(T.window); end

        m = (strcmpi(T.marker, opts.marker) & strcmpi(T.subset, opts.subset) & strcmpi(T.window, opts.window));
        if any(m)
            r = T(find(m,1,'first'),:);
            rows(end+1,:) = {string(sessname), string(opts.bodypart), string(opts.marker), ...
                string(opts.subset), string(opts.window), ...
                r.AUROC, r.nTar, r.nRef}; %#ok<AGROW>
            found = true;
            break
        end
    end

    if ~found
        % keep a NaN row so joins still work if you want
        rows(end+1,:) = {string(sessname), string(opts.bodypart), string(opts.marker), ...
            string(opts.subset), string(opts.window), ...
            nan, nan, nan}; %#ok<AGROW>
    end
end

if isempty(rows)
    BRES.T = table();
    return
end

BRES.T = cell2table(rows, 'VariableNames', ...
    {'session','bodypart','marker','subset','window','AUROC','nTar','nRef'});
end
