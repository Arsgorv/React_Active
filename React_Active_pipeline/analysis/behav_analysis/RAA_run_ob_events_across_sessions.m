function RES = RAA_run_ob_events_across_sessions(sessions, opts)
% RAA_run_ob_events_across_sessions
% Collect OB event metrics across sessions from ob_event_analysis outputs.
%
% OUTPUT:
%   RES.T : long table (session-level rows)
%   RES.S : summary table across sessions (mean/sem across sessions)
%   RES.bySession : per-session bookkeeping (paths/files/loaded snapshots)
%
% SAVES:
%   <saveRoot>/analysis/OB_events_across_sessions.mat
%   <saveRoot>/analysis/OB_events_across_sessions.csv
%   <saveRoot>/analysis/OB_events_across_sessions_summary.csv

if nargin < 2, opts = struct(); end
if ~iscell(sessions), error('sessions must be a cell array'); end

% canonical bands for summary views (PLV bandmean + optional checks)
bands = [0.5 4; 4 6; 40 60; 60 80];
bandNames = {'delta','theta','gamma','highgamma'};

rows = {};
byS = struct([]);

for s = 1:numel(sessions)
    dp = sessions{s};
    if ~exist(dp,'dir')
        warning('Missing session folder: %s', dp);
        continue
    end

    [~,sessname] = fileparts(dp);
    outDir = fullfile(dp,'analysis','ob_events');
    if ~exist(outDir,'dir')
        warning('Missing outDir: %s', outDir);
        continue
    end

    % ---- find files robustly (don’t assume exact prefix = sessname)
    F = struct();
    F.lowSpec = pick_one_file(outDir, '*_B_LowEvent_Spectrum_*obSpecMetrics.mat');
    F.midSpec = pick_one_file(outDir, '*_B_Middle_Spectrum_*obSpecMetrics.mat');

    F.itpcLow = pick_one_file(outDir, '*_ITPC_*LOW*_metrics.mat');
    F.itpcMid = pick_one_file(outDir, '*_ITPC_*MID*_metrics.mat');

    F.plvStim = pick_one_file(outDir, '*_PLV_*metrics.mat');

    S = struct();
    S.sessname = sessname;
    S.datapath = dp;

    if ~isempty(F.lowSpec),  S.lowSpec  = load(F.lowSpec);  else, S.lowSpec  = []; end
    if ~isempty(F.midSpec),  S.midSpec  = load(F.midSpec);  else, S.midSpec  = []; end
    if ~isempty(F.itpcLow),  S.itpcLow  = load(F.itpcLow);  else, S.itpcLow  = []; end
    if ~isempty(F.itpcMid),  S.itpcMid  = load(F.itpcMid);  else, S.itpcMid  = []; end
    if ~isempty(F.plvStim),  S.plvStim  = load(F.plvStim);  else, S.plvStim  = []; end

    byS(end+1).sessname = sessname; %#ok<AGROW>
    byS(end).datapath = dp;
    byS(end).files = F;
    byS(end).loaded = S;

    % ---- SPEC metrics: per-window means/sem/slopes + postpre
    if ~isempty(S.lowSpec) && isfield(S.lowSpec,'Met')
        rows = [rows; ra_rows_from_specMet(sessname, 'LOW', S.lowSpec.Met, bandNames)]; %#ok<AGROW>
    end
    if ~isempty(S.midSpec) && isfield(S.midSpec,'Met')
        rows = [rows; ra_rows_from_specMet(sessname, 'MID', S.midSpec.Met, bandNames)]; %#ok<AGROW>
    end

    % ---- ITPC metrics
    if ~isempty(S.itpcLow) && isfield(S.itpcLow,'MetITPC')
        rows = [rows; ra_rows_from_itpcMet(sessname, 'ITPC_LOW', S.itpcLow.MetITPC)]; %#ok<AGROW>
    end
    if ~isempty(S.itpcMid) && isfield(S.itpcMid,'MetITPC')
        rows = [rows; ra_rows_from_itpcMet(sessname, 'ITPC_MID', S.itpcMid.MetITPC)]; %#ok<AGROW>
    end

    % ---- PLV band means (canonical bands)
    if ~isempty(S.plvStim) && isfield(S.plvStim,'MetPLV')
        rows = [rows; ra_rows_from_plvMet(sessname, 'PLV_stimOn', S.plvStim.MetPLV, bands, bandNames)]; %#ok<AGROW>
    end
end

T = ra_rows_to_table(rows);

% ---- across-session summaries (mean/sem over sessions)
Ssummary = ra_summarize_across_sessions(T);

RES = struct();
RES.T = T;
RES.S = Ssummary;
RES.bySession = byS;

% ---- choose save root
saveRoot = '';
if isfield(opts,'saveRoot') && ~isempty(opts.saveRoot)
    saveRoot = opts.saveRoot;
else
    % default: parent folder of first session
    saveRoot = fileparts(sessions{1});
end

saveDir = fullfile(saveRoot,'analysis');
if ~exist(saveDir,'dir'), mkdir(saveDir); end

save(fullfile(saveDir,'OB_events_across_sessions.mat'), 'RES', 'opts');

try, writetable(T,        fullfile(saveDir,'OB_events_across_sessions.csv'));         catch, end
try, writetable(Ssummary, fullfile(saveDir,'OB_events_across_sessions_summary.csv')); catch, end

end

% ========================= helpers =========================

function fn = pick_one_file(outDir, pattern)
d = dir(fullfile(outDir, pattern));
if isempty(d)
    fn = '';
    return
end
% choose the newest (or the largest; newest is usually safer)
[~,ix] = max([d.datenum]);
fn = fullfile(d(ix).folder, d(ix).name);
end

function rows = ra_rows_from_specMet(sessname, tag, Met, bandNames)

rows = {};
if ~isfield(Met,'group') || isempty(Met.group), return, end

hasMetricBands = isfield(Met,'metricBands') && ~isempty(Met.metricBands) && isnumeric(Met.metricBands);

for g = 1:numel(Met.group)
    if ~isfield(Met.group(g),'name'), continue, end
    gname = char(Met.group(g).name);

    if ~isfield(Met.group(g),'band') || isempty(Met.group(g).band), continue, end
    for bb = 1:numel(Met.group(g).band)

        br = double(Met.group(g).band(bb).range(:)'); % [f1 f2]
        if numel(br) >= 2
            br = br(1:2);
        else
            continue
        end

        if hasMetricBands
            bidx = find_band_idx_tol(br, Met.metricBands, 1e-6);
        else
            bidx = [];
        end

        if isempty(bidx) || bidx > numel(bandNames)
            bname = sprintf('%.3g_%.3gHz', br(1), br(2));
        else
            bname = bandNames{bidx};
        end

        % per-window metrics
        if isfield(Met.group(g).band(bb),'win') && ~isempty(Met.group(g).band(bb).win)
            for ww = 1:numel(Met.group(g).band(bb).win)
                W = Met.group(g).band(bb).win(ww);

                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, char(W.name), 'mean',       double(W.mean)}; %#ok<AGROW>
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, char(W.name), 'sem',        double(W.sem)}; %#ok<AGROW>
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, char(W.name), 'slope_mean', double(W.slope_mean)}; %#ok<AGROW>
                rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, char(W.name), 'n',          double(W.n)}; %#ok<AGROW>
            end
        end

        % post-pre (W3-W2)
        if isfield(Met.group(g).band(bb),'postpre_mean')
            rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, 'postpre', 'mean', double(Met.group(g).band(bb).postpre_mean)}; %#ok<AGROW>
            rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, 'postpre', 'sem',  double(Met.group(g).band(bb).postpre_sem)};  %#ok<AGROW>
            rows(end+1,:) = {sessname, tag, 'SPEC', gname, bname, 'postpre', 'n',    double(Met.group(g).band(bb).postpre_n)};    %#ok<AGROW>
        end
    end
end

end


function rows = ra_rows_from_itpcMet(sessname, tag, MetITPC)

rows = {};
if ~isfield(MetITPC,'group') || isempty(MetITPC.group), return, end

for g = 1:numel(MetITPC.group)
    if ~isfield(MetITPC.group(g),'name'), continue, end
    gname = char(MetITPC.group(g).name);

    if isfield(MetITPC.group(g),'mean_post_0_1s')
        rows(end+1,:) = {sessname, tag, 'ITPC', gname, '', 'post_0_1s', 'mean', double(MetITPC.group(g).mean_post_0_1s)}; %#ok<AGROW>
    end

    if isfield(MetITPC.group(g),'band1_post_0_1s')
        rows(end+1,:) = {sessname, tag, 'ITPC', gname, 'band1', 'post_0_1s', 'mean', double(MetITPC.group(g).band1_post_0_1s)}; %#ok<AGROW>
        rows(end+1,:) = {sessname, tag, 'ITPC', gname, 'band1', 'peak',      'mean', double(MetITPC.group(g).band1_peak)};      %#ok<AGROW>
    end
    if isfield(MetITPC.group(g),'band2_post_0_1s')
        rows(end+1,:) = {sessname, tag, 'ITPC', gname, 'band2', 'post_0_1s', 'mean', double(MetITPC.group(g).band2_post_0_1s)}; %#ok<AGROW>
        rows(end+1,:) = {sessname, tag, 'ITPC', gname, 'band2', 'peak',      'mean', double(MetITPC.group(g).band2_peak)};      %#ok<AGROW>
    end

    if isfield(MetITPC.group(g),'nEvents')
        rows(end+1,:) = {sessname, tag, 'ITPC', gname, '', 'nEvents', 'n', double(MetITPC.group(g).nEvents)}; %#ok<AGROW>
    end
end

end

function rows = ra_rows_from_plvMet(sessname, tag, MetPLV, bands, bandNames)

rows = {};
if ~isfield(MetPLV,'group') || isempty(MetPLV.group), return, end
if ~isfield(MetPLV,'fCenters') || isempty(MetPLV.fCenters), return, end

f = double(MetPLV.fCenters(:)'); % IMPORTANT: often uint8

for g = 1:numel(MetPLV.group)
    if ~isfield(MetPLV.group(g),'name'), continue, end
    gname = char(MetPLV.group(g).name);

    if ~isfield(MetPLV.group(g),'plv') || isempty(MetPLV.group(g).plv)
        continue
    end
    plv = double(MetPLV.group(g).plv(:)');

    for bb = 1:size(bands,1)
        br = bands(bb,:);
        m = nanmean(plv(f>=br(1) & f<=br(2)));
        rows(end+1,:) = {sessname, tag, 'PLV', gname, bandNames{bb}, 'bandmean', 'mean', double(m)}; %#ok<AGROW>
    end

    if isfield(MetPLV.group(g),'nEvents')
        rows(end+1,:) = {sessname, tag, 'PLV', gname, '', 'nEvents', 'n', double(MetPLV.group(g).nEvents)}; %#ok<AGROW>
    end
end

end

function T = ra_rows_to_table(rows)
if isempty(rows)
    T = table();
    return
end
T = cell2table(rows, 'VariableNames', ...
    {'session','tag','modality','group','band','window','metric','value'});
end

function idx = find_band_idx_tol(br, metricBands, tol)
% Robust match of a [f1 f2] band range against Met.metricBands.
% Handles metricBands as:
%   - n x 2 numeric
%   - vector [f1 f2 f1 f2 ...]
%   - n x m numeric with m>2 (uses first 2 cols)
% Returns [] if cannot interpret.

idx = [];

% br -> 1x2
br = double(br(:)');
if numel(br) < 2, return, end
br = br(1:2);

% metricBands may be missing / non-numeric
if nargin < 2 || isempty(metricBands) || ~isnumeric(metricBands)
    return
end
mb = double(metricBands);

% Normalize mb to n x 2
if isvector(mb)
    mb = mb(:);
    if mod(numel(mb),2) ~= 0
        return
    end
    mb = reshape(mb, 2, []).';
else
    % if it has >=2 cols, keep first 2; if 1 col, try reshape
    if size(mb,2) >= 2
        mb = mb(:,1:2);
    elseif size(mb,2) == 1
        if mod(numel(mb),2) ~= 0
            return
        end
        mb = reshape(mb(:), 2, []).';
    else
        return
    end
end

% Now mb is n x 2
for k = 1:size(mb,1)
    if all(abs(br - mb(k,:)) <= tol)
        idx = k;
        return
    end
end
end


function S = ra_summarize_across_sessions(T)
% Summary over sessions for each unique (tag, modality, group, band, window, metric)
% Only summarizes metric=='mean' rows; ignores 'sem' rows (those are per-session SEMs).

if isempty(T) || height(T)==0
    S = table();
    return
end

% keep only mean-like rows to summarize across sessions
keep = strcmp(T.metric,'mean') | strcmp(T.metric,'slope_mean');
T0 = T(keep,:);

% group keys
K = {'tag','modality','group','band','window','metric'};
[G,keys] = findgroups(T0(:,K));

mSess = splitapply(@(x) mean(x,'omitnan'), T0.value, G);
sdSess = splitapply(@(x) std(x,'omitnan'),  T0.value, G);
nSess  = splitapply(@(x) sum(~isnan(x)),    T0.value, G);
semSess = sdSess ./ sqrt(max(nSess,1));

S = keys;
S.nSessions = nSess;
S.meanAcrossSessions = mSess;
S.semAcrossSessions  = semSess;

end
