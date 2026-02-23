function RES = RAA_collect_resp_coupling_across_sessions(sessions, opts)
% RAA_collect_resp_coupling_across_sessions
% Collect respiration–OB coupling metrics across sessions from ob_event_analysis outputs.
%
% Expects per session:
%   <datapath>/analysis/ob_events/*RESP_overlay_stimOn_summary.csv
%   <datapath>/analysis/ob_events/*RESP_coupling_trials.csv
%
% Output:
%   RES.Ttrial : long per-trial table
%   RES.Ssess  : per-session summary (band x window)
%   RES.SessInfo : session bookkeeping table

if nargin < 2, opts = struct(); end
if ~isfield(opts,'alignName'), opts.alignName = 'stimOn'; end
if ~isfield(opts,'tryLoadBaphyForTrialType'), opts.tryLoadBaphyForTrialType = false; end

Ttrial = table();
Ssess  = table();
SessInfo = table();

for s = 1:numel(sessions)
    dp = sessions{s};
    outDir = fullfile(dp,'analysis','ob_events');
    if ~exist(outDir,'dir'), continue, end
    [~,sessname] = fileparts(dp);

    fSum = dir(fullfile(outDir, sprintf('*RESP_overlay_%s_summary.csv', opts.alignName)));
    fTr  = dir(fullfile(outDir, sprintf('*RESP_coupling_trials.csv')));

    if ~isempty(fSum)
        Ts = readtable(fullfile(fSum(1).folder, fSum(1).name));
        Ts.session = string(Ts.session);
        Ts.band    = string(Ts.band);
        Ts.window  = string(Ts.window);
        Ssess = [Ssess; Ts]; %#ok<AGROW>
    end

    if ~isempty(fTr)
        Tt = readtable(fullfile(fTr(1).folder, fTr(1).name));
        % normalize types
        if ismember('session', Tt.Properties.VariableNames)
            Tt.session = string(Tt.session);
        else
            Tt.session = repmat(string(sessname), height(Tt), 1);
        end
        if ismember('band', Tt.Properties.VariableNames),   Tt.band   = string(Tt.band);   end
        if ismember('window', Tt.Properties.VariableNames), Tt.window = string(Tt.window); end

        % add trialType if missing and user wants it (loads Baphy once per session)
        if ~ismember('trialType', Tt.Properties.VariableNames) && opts.tryLoadBaphyForTrialType
            B = ra_load_baphy(dp);
            if ~isempty(B) && isfield(B,'trial') && isfield(B.trial,'type')
                tt = string(B.trial.type(:));
                tri = Tt.trialIdx;
                tri = double(tri(:));
                good = tri>=1 & tri<=numel(tt);
                Tt.trialType = strings(height(Tt),1);
                Tt.trialType(good) = tt(tri(good));
            end
        end

        Ttrial = [Ttrial; Tt]; %#ok<AGROW>
    end

    SessInfo = [SessInfo; table(string(sessname), string(dp), ...
        'VariableNames', {'session','datapath'})]; %#ok<AGROW>
end

RES = struct();
RES.Ttrial = Ttrial;
RES.Ssess  = Ssess;
RES.SessInfo = SessInfo;
RES.opts_snapshot = opts;

% If summary CSV was missing, compute session summaries from trial table.
if isempty(Ssess) && ~isempty(Ttrial)
    RES.Ssess = ra_compute_resp_session_summary_from_trials(Ttrial);
end

end

% ---------------- helpers ----------------
function Baphy = ra_load_baphy(datapath)
Baphy = [];
f1 = fullfile(datapath,'Master_sync.mat');
f2 = fullfile(datapath,'Baphy_RA.mat');
try
    if exist(f1,'file')
        S = load(f1,'Baphy'); Baphy = S.Baphy;
    elseif exist(f2,'file')
        S = load(f2,'Baphy'); Baphy = S.Baphy;
    end
catch
    Baphy = [];
end
end

function S = ra_compute_resp_session_summary_from_trials(Ttrial)
% Summarize by session x band x window
Ttrial.session = string(Ttrial.session);
if ismember('piezoMode', Ttrial.Properties.VariableNames), Ttrial.piezoMode = string(Ttrial.piezoMode); end
if ismember('group',    Ttrial.Properties.VariableNames), Ttrial.group    = string(Ttrial.group);    end
if ismember('spec',     Ttrial.Properties.VariableNames), Ttrial.spec     = string(Ttrial.spec);     end
T.band   = string(T.band);
T.window = string(T.window);

[G, keys] = findgroups(Ttrial.session, Ttrial.band, Ttrial.window);

S = table();
S.session = keys(:,1);
S.band    = keys(:,2);
S.window  = keys(:,3);

nUsed = splitapply(@(x) sum(isfinite(x)), Ttrial.peakR, G);
medLag = splitapply(@(x) median(x,'omitnan'), Ttrial.peakLag_s, G);
iqrLag = splitapply(@(x) iqr(x(isfinite(x))), Ttrial.peakLag_s, G);
meanR  = splitapply(@(x) mean(x,'omitnan'), Ttrial.peakR, G);
medR   = splitapply(@(x) median(x,'omitnan'), Ttrial.peakR, G);

S.nUsed = nUsed;
S.medianLag_s = medLag;
S.iqrLag_s = iqrLag;
S.meanR = meanR;
S.medianR = medR;

end
