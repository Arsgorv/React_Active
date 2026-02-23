function RES = RAA_collect_ob_trialmetrics_across_sessions(sessions, opts)
% Collect OB per-trial metrics tables across sessions (LOW+MID),
% and compute session-level summaries + heatmap-ready matrices.

if nargin < 2, opts = struct(); end
if ~isfield(opts,'alignName'), opts.alignName = 'stimOn'; end

Tall = table();
SessInfo = table();

for s = 1:numel(sessions)
    dp = sessions{s};
    outDir = fullfile(dp,'analysis','ob_events');
    if ~exist(outDir,'dir'), continue, end
    [~,sessname] = fileparts(dp);

    Fl = dir(fullfile(outDir, sprintf('*_OB_trialMetrics_%s_LOW.csv', opts.alignName)));
    Fm = dir(fullfile(outDir, sprintf('*_OB_trialMetrics_%s_MID.csv', opts.alignName)));

    T = table();
    if ~isempty(Fl), T = [T; readtable(fullfile(Fl(1).folder, Fl(1).name))]; end %#ok<AGROW>
    if ~isempty(Fm), T = [T; readtable(fullfile(Fm(1).folder, Fm(1).name))]; end %#ok<AGROW>
    if isempty(T), continue, end

    % enforce types
    T.session = string(T.session);
    T.align   = string(T.align);
    T.trialType = string(T.trialType);
    T.band = string(T.band);
    T.window = string(T.window);

    Tall = [Tall; T]; %#ok<AGROW>

    SessInfo = [SessInfo; table(string(sessname), string(dp), 'VariableNames',{'session','datapath'})]; %#ok<AGROW>
end

Tall.window = lower(string(Tall.window));
Tall.window = replace(Tall.window, "stimon_to_stimoff", "stimon_to_stimoff");
Tall.window = replace(Tall.window, "stimon_to_stimoff", "stimon_to_stimoff"); % keep
Tall.window = replace(Tall.window, "stimon_to_stimoff", "stimon_to_stimoff"); % no-op

Tall.window = replace(Tall.window, "stimon_to_stimoff", "stimon_to_stimoff");
Tall.window = replace(Tall.window, "stimoff_to_arrival", "stimoff_to_arrival");
Tall.window = replace(Tall.window, "stimoff_to_arr", "stimoff_to_arrival");
Tall.window = replace(Tall.window, "arr_to_stop", "arr_to_stop");
Tall.window = replace(Tall.window, "stimon_to_arrival", "stimon_to_arrival");
Tall.window = replace(Tall.window, "arrival_m500_to_arrival", "arrival_m500_to_arrival");


RES = struct();
RES.Tall = Tall;
RES.SessInfo = SessInfo;
RES.opts_snapshot = opts;

if isempty(Tall)
    RES.Sess = table();
    RES.Heat = struct();
    return
end

% ---- define groups from trialType (robust)
tt = lower(Tall.trialType);
Tall.isTar = strcmp(tt,'target');
Tall.isRef = strcmp(tt,'reference') | strcmp(tt,'ref');
Tall.isNoSound = strcmp(tt,'nosound') | strcmp(tt,'no_sound');
Tall.isNoMotor = strcmp(tt,'nomotor') | strcmp(tt,'no_motor');

% ---- session-level tar/ref effect metrics for each band x window
[Gid, keySess, keyBand, keyWin] = findgroups(Tall.session, Tall.band, Tall.window);

S = table();
S.session = keySess;
S.band    = keyBand;
S.window  = keyWin;

nGroups = max(Gid);
nTar = zeros(nGroups,1);
nRef = zeros(nGroups,1);
AUROC = nan(nGroups,1);
CohenD = nan(nGroups,1);
meanDiff = nan(nGroups,1);

for k = 1:nGroups
    Tk = Tall(Gid==k,:);
    xT = Tk.meanPower(Tk.isTar);
    xR = Tk.meanPower(Tk.isRef);

    xT = xT(isfinite(xT));
    xR = xR(isfinite(xR));

    nTar(k) = numel(xT);
    nRef(k) = numel(xR);

    if nTar(k) < 5 || nRef(k) < 5
        continue
    end

    M = effect_metrics(xT, xR);   % must match behaviour_analysis convention
    AUROC(k)   = M.AUROC;
    CohenD(k)  = M.CohensD;
    meanDiff(k)= M.meanDiff;
end

S.nTar = nTar;
S.nRef = nRef;
S.AUROC = AUROC;
S.CohenD = CohenD;
S.meanDiff = meanDiff;

RES.Sess = S;
% ---- heatmap-ready matrices (meanPower per session for each group)
RES.Heat = ra_build_heat_struct(Tall);

end

function Heat = ra_build_heat_struct(Tall)
% Heatmaps: session x timeWindow for each band and group.
% Windows are categorical names; you can map to fixed order.

winOrder = ["stimOn_to_stimOff","stimOff_to_arr","arr_to_stop"];
bands = unique(Tall.band);

sessions = unique(Tall.session,'stable');

Heat = struct();
Heat.sessions = sessions;
Heat.winOrder = winOrder;
Heat.bands = bands;

% group masks
gNames = ["all","tar","ref","nosound","nomotor"];
for bi = 1:numel(bands)
    b = bands(bi);
    for gi = 1:numel(gNames)
        g = gNames(gi);
        M = nan(numel(sessions), numel(winOrder));
        for si = 1:numel(sessions)
            sess = sessions(si);
            Ts = Tall(Tall.session==sess & Tall.band==b,:);
            for wi = 1:numel(winOrder)
                w = winOrder(wi);
                Tw = Ts(Ts.window==w,:);
                if isempty(Tw), continue, end

                switch g
                    case "all"
                        x = Tw.meanPower;
                    case "tar"
                        x = Tw.meanPower(Tw.isTar);
                    case "ref"
                        x = Tw.meanPower(Tw.isRef);
                    case "nosound"
                        x = Tw.meanPower(Tw.isNoSound);
                    case "nomotor"
                        x = Tw.meanPower(Tw.isNoMotor);
                end
                if sum(isfinite(x)) >= 3
                    M(si,wi) = mean(x,'omitnan');
                end
            end
        end
        Heat.(char(b)).(char(g)) = M;
    end
end
end
