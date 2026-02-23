function OUT = ob_bands_behaviour_analysis(datapath, opts)
% ob_bands_behaviour_analysis
% Behaviour_analysis-style per-session OB band figure + metrics across:
%   windows x subsets, Target vs Reference
%
% Depends on:
%   - load_baphy(datapath)
%   - B_Low_Spectrum.mat (required)
%   - B_Middle_Spectrum.mat (optional, needed for gamma/highgamma)
%   - MakeSpreadAndBoxPlot3_SB
%
% Saves into: <datapath>/analysis/ob_events/OB-bands/

if nargin < 2, opts = struct(); end
if ~isfield(opts,'smoothing_win_s'),      opts.smoothing_win_s = 0.10; end
if ~isfield(opts,'spec_twin_whole_s'),    opts.spec_twin_whole_s = [-0.5 7.0]; end
if ~isfield(opts,'spec_baseline_s'),      opts.spec_baseline_s = [-0.5 -0.1]; end
if ~isfield(opts,'spec_freq_xlim_low'),   opts.spec_freq_xlim_low = [0.5 12]; end
if ~isfield(opts,'spec_freq_xlim_mid'),   opts.spec_freq_xlim_mid = [20 100]; end
if ~isfield(opts,'min_trials_per_group'), opts.min_trials_per_group = 3; end

if ~isfield(opts,'ob_band_delta'),        opts.ob_band_delta = [0.5 4]; end
if ~isfield(opts,'ob_band_theta'),        opts.ob_band_theta = [4 12]; end
if ~isfield(opts,'ob_band_gamma'),        opts.ob_band_gamma = [40 60]; end
if ~isfield(opts,'ob_band_highgamma'),    opts.ob_band_highgamma = [60 80]; end

% output folders
outBase = fullfile(datapath,'analysis','ob_events','OB-bands');
if ~exist(outBase,'dir'), mkdir(outBase); end
[~,sessname] = fileparts(datapath);

% Load Baphy
Baphy = load_baphy(datapath);
nTr = Baphy.n_trials;

% Absolute stim-on per trial
stim_on_abs_s = Baphy.trial.abs_stim_start_s(:);

% Trial masks (same logic as behaviour_analysis)
mskGood = true(nTr,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'goodTrials')
    mskGood = false(nTr,1); mskGood(Baphy.idx.goodTrials(:)) = true;
end

mskTar = false(nTr,1); mskRef = false(nTr,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'Target')
    mskTar(Baphy.idx.Target(:)) = true;
    mskRef(Baphy.idx.Reference(:)) = true;
elseif isfield(Baphy.trial,'type')
    mskTar(strcmpi(Baphy.trial.type,'Target')) = true;
    mskRef(strcmpi(Baphy.trial.type,'Reference')) = true;
end

mskNS = false(nTr,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'NoSound')
    mskNS(Baphy.idx.NoSound(:)) = true;
elseif isfield(Baphy.trial,'is_nosound')
    mskNS = logical(Baphy.trial.is_nosound(:));
end

mskNM = false(nTr,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'NomotorExclusive')
    mskNM(Baphy.idx.NomotorExclusive(:)) = true;
elseif isfield(Baphy.trial,'is_nomotor')
    mskNM = logical(Baphy.trial.is_nomotor(:));
end

mskReg = ~(mskNS | mskNM);

subsets = struct();
subsets(1).name = 'regular';
subsets(1).idxA = find(mskTar & mskGood & mskReg);
subsets(1).idxB = find(mskRef & mskGood & mskReg);

subsets(2).name = 'nosound';
subsets(2).idxA = find(mskTar & mskGood & mskNS);
subsets(2).idxB = find(mskRef & mskGood & mskNS);

subsets(3).name = 'nomotor';
subsets(3).idxA = find(mskTar & mskGood & mskNM);
subsets(3).idxB = find(mskRef & mskGood & mskNM);

% Windows (match behaviour_analysis)
winNames = {'stimon_to_stimoff','stimon_to_arrival','stimoff_to_arrival','arrival_m500_to_arrival'};

% Event times relative to stimOn (for lines/patch)
EvRel = ra_get_trial_rel_times(Baphy);
[winList, w0_all, w1_all] = ra_windows_from_EvRel(EvRel);
assert(numel(winList)==numel(winNames));

% Load spectrograms
lowFile = fullfile(datapath,'ephys','B_Low_Spectrum.mat');
if ~exist(lowFile,'file'), error('Missing: %s', lowFile); end
Sl = load(lowFile,'Spectro'); SpectroLow = Sl.Spectro;

[tL, fL, PLraw] = ra_extract_spectro_triplet(SpectroLow);
[t_abs_L, P_tf_L] = ra_standardize_spectro_time_and_orient(tL, fL, PLraw);
P_tf_L = log10(max(double(P_tf_L), eps));
dtL = median(diff(t_abs_L));
tvec = (opts.spec_twin_whole_s(1):dtL:opts.spec_twin_whole_s(2))';
bmask = (tvec >= opts.spec_baseline_s(1)) & (tvec <= opts.spec_baseline_s(2));
swSamp = max(1, round(opts.smoothing_win_s / dtL));

fmaskL = (fL >= opts.spec_freq_xlim_low(1)) & (fL <= opts.spec_freq_xlim_low(2));
f_plot_L = fL(fmaskL);
V_L = P_tf_L(:, fmaskL);

% MID
hasMID = false;
V_M = []; t_abs_M = []; f_plot_M = [];
midFile = fullfile(datapath,'ephys','B_Middle_Spectrum.mat');
if exist(midFile,'file')
    Sm = load(midFile,'Spectro');
    if isfield(Sm,'Spectro') && ~isempty(Sm.Spectro)
        [tM, fM, PMraw] = ra_extract_spectro_triplet(Sm.Spectro);
        [t_abs_M, P_tf_M] = ra_standardize_spectro_time_and_orient(tM, fM, PMraw);
        P_tf_M = log10(max(double(P_tf_M), eps));
        fmaskM = (fM >= opts.spec_freq_xlim_mid(1)) & (fM <= opts.spec_freq_xlim_mid(2));
        f_plot_M = fM(fmaskM);
        V_M = P_tf_M(:, fmaskM);
        hasMID = true;
    end
end

% Marker list
markers = {'OB_delta','OB_theta','OB_gamma','OB_highgamma'};
bodypart = 'OB-bands';

% Table vars (match behaviour_analysis pattern + include raw scalars in MAT)
varNames = {'bodypart','marker','subset','window', ...
    'nTar','nRef', ...
    'AshmanD','AUROC','CohensD','meanDiff', ...
    'meanRef','meanTar','stdRef','stdTar', ...
    'rocFPR','rocTPR', ...
    'trialIdxTar','trialIdxRef','xTar','xRef'};

% Colors consistent everywhere
COL = struct();
COL.Ref = [0.10 0.55 0.90];
COL.Tar = [0.95 0.45 0.10];
COL.StimPatch = [0.70 0.70 0.70];
COL.AntPatch  = [0.00 1.00 0.00];
COL.SpoutLine = [0.00 0.00 0.00];
COL.RewLine   = [0.85 0.10 0.10];
COL.DiagROC   = [0.40 0.40 0.40];
ColsRT = {COL.Ref, COL.Tar};
LegRT  = {'Ref','Tar'};

OUT = struct();
OUT.session = sessname;
OUT.saved = {};

for ss = 1:numel(subsets)
    subsetName = subsets(ss).name;
    idxTar = subsets(ss).idxA(:);
    idxRef = subsets(ss).idxB(:);

    % only trials with stim_on
    idxTar = idxTar(isfinite(stim_on_abs_s(idxTar)));
    idxRef = idxRef(isfinite(stim_on_abs_s(idxRef)));

    for w = 1:numel(winNames)
        wName = winNames{w};
        wi = find(strcmpi(winList, wName), 1);
        if isempty(wi), error('Window not found: %s', wName); end

        f = figure('Color','w','Visible','on');
        set(f,'Position',[50 50 1900 max(420, 220*numel(markers))]);
        sgtitle(sprintf('%s | %s | %s | sw=%.3fs | Tar=%d Ref=%d', ...
            bodypart, subsetName, wName, opts.smoothing_win_s, numel(idxTar), numel(idxRef)), ...
            'Interpreter','none');

        rows = cell(0, numel(varNames));

        pooled = [idxTar; idxRef];
        t_stimOff = median(EvRel.stimOff_rel(pooled),'omitnan');
        t_arr     = median(EvRel.arr_rel(pooled),'omitnan');
        t_reward  = median(EvRel.reward_rel(pooled),'omitnan');
        w0_med    = median(w0_all(pooled,wi),'omitnan');
        w1_med    = median(w1_all(pooled,wi),'omitnan');

        for m = 1:numel(markers)
            mk = markers{m};

            [useMid, br] = ra_band_range_from_marker(mk, opts);
            if useMid && ~hasMID
                rows(end+1,:) = {bodypart, mk, subsetName, wName, ...
                    numel(idxTar), numel(idxRef), nan,nan,nan,nan,nan,nan,nan,nan,{[]},{[]}, ...
                    {idxTar},{idxRef},{[]},{[]}};
                continue
            end

            if useMid
                V = V_M; t_abs = t_abs_M; f_plot = f_plot_M;
            else
                V = V_L; t_abs = t_abs_L; f_plot = f_plot_L;
            end

            % band traces per trial (baseline-sub per trial)
            BTar = ra_bandtrace_from_spec_to_grid(V, t_abs, f_plot, stim_on_abs_s(idxTar), tvec, bmask, br);
            BRef = ra_bandtrace_from_spec_to_grid(V, t_abs, f_plot, stim_on_abs_s(idxRef), tvec, bmask, br);

            % scalar per trial using trial-specific window edges
            xTar = ra_trial_scalar_in_window(BTar, tvec, w0_all(idxTar,wi), w1_all(idxTar,wi));
            xRef = ra_trial_scalar_in_window(BRef, tvec, w0_all(idxRef,wi), w1_all(idxRef,wi));

            doStats = (numel(xTar) >= opts.min_trials_per_group) && (numel(xRef) >= opts.min_trials_per_group);
            if doStats
                met = ra_effect_metrics_simple(xTar, xRef);
                met.meanRef = mean(xRef,'omitnan');
                met.meanTar = mean(xTar,'omitnan');
                met.stdRef  = std(xRef,[],'omitnan');
                met.stdTar  = std(xTar,[],'omitnan');
            else
                met = struct('AshmanD',nan,'AUROC',nan,'CohensD',nan,'meanDiff',nan, ...
                    'rocFPR',[],'rocTPR',[], 'meanRef',nan,'meanTar',nan,'stdRef',nan,'stdTar',nan);
            end

            % PSTH mean±SEM (smoothed)
            meanTar = movmean(mean(BTar,1,'omitnan'), swSamp, 'omitnan');
            meanRef = movmean(mean(BRef,1,'omitnan'), swSamp, 'omitnan');
            semTar  = movmean(std(BTar,[],1,'omitnan'), swSamp, 'omitnan') ./ max(1,sqrt(size(BTar,1)));
            semRef  = movmean(std(BRef,[],1,'omitnan'), swSamp, 'omitnan') ./ max(1,sqrt(size(BRef,1)));

            lo = min([meanTar-semTar, meanRef-semRef]); hi = max([meanTar+semTar, meanRef+semRef]);
            if ~isfinite(lo), lo = -1; end
            if ~isfinite(hi), hi =  1; end
            if lo == hi, lo = lo-1; hi = hi+1; end
            pad = 0.05*(hi-lo);
            yl = [lo-pad hi+pad];

            row0 = (m-1)*6;

            % 1) PSTH
            subplot(numel(markers),6,row0+1); hold on
            hStim=[]; hAnt=[]; hSpout=[]; hRew=[];

            if isfinite(t_stimOff) && t_stimOff > 0
                hStim = patch([0 t_stimOff t_stimOff 0], [yl(1) yl(1) yl(2) yl(2)], ...
                    COL.StimPatch,'EdgeColor','none','FaceAlpha',0.18);
                set(hStim,'HandleVisibility','on');
            end
            if isfinite(w0_med) && isfinite(w1_med) && w1_med > w0_med
                hAnt = patch([w0_med w1_med w1_med w0_med], yl([1 1 2 2]), ...
                    COL.AntPatch,'EdgeColor','none','FaceAlpha',0.08);
                set(hAnt,'HandleVisibility','on');
            end

            hRef = ra_plot_mean_sem_col(tvec, meanRef, semRef, COL.Ref);
            hTar = ra_plot_mean_sem_col(tvec, meanTar, semTar, COL.Tar);

            xline(0,'k-','LineWidth',1.2);
            if isfinite(t_arr),    hSpout = xline(t_arr,'--','LineWidth',1.2,'Color',COL.SpoutLine); set(hSpout,'HandleVisibility','on'); end
            if isfinite(t_reward), hRew   = xline(t_reward,'--','LineWidth',1.2,'Color',COL.RewLine);   set(hRew,'HandleVisibility','on'); end
            yline(0,'k:');

            ylabel(mk,'Interpreter','none');
            if m==numel(markers), xlabel('Time from stimOn (s)'); end
            xlim([tvec(1) tvec(end)]); ylim(yl);
            title(sprintf('%s (Tar=%d Ref=%d)', mk, numel(xTar), numel(xRef)), 'Interpreter','none');
            box off

            if m==1
                hh = gobjects(0); ll = {};
                if ~isempty(hStim) && isgraphics(hStim),  hh(end+1)=hStim;  ll{end+1}='Stimulus'; end
                if ~isempty(hAnt)  && isgraphics(hAnt),   hh(end+1)=hAnt;   ll{end+1}='Anticip window'; end
                if isgraphics(hRef), hh(end+1)=hRef; ll{end+1}='Reference'; end
                if isgraphics(hTar), hh(end+1)=hTar; ll{end+1}='Target'; end
                if ~isempty(hSpout) && isgraphics(hSpout), hh(end+1)=hSpout; ll{end+1}='Spout'; end
                if ~isempty(hRew)   && isgraphics(hRew),   hh(end+1)=hRew;   ll{end+1}='Reward'; end
                lgd = legend(hh,ll,'Location','westoutside','Box','off');
                set(lgd,'Interpreter','none');
            end
            hold off

            % 2) Histogram+KDE
            subplot(numel(markers),6,row0+2); hold on
            if ~isempty(xTar), histogram(xTar,30,'Normalization','pdf','FaceColor',COL.Tar,'FaceAlpha',0.25,'EdgeColor','none'); end
            if ~isempty(xRef), histogram(xRef,30,'Normalization','pdf','FaceColor',COL.Ref,'FaceAlpha',0.25,'EdgeColor','none'); end
            if numel(xTar) >= 7 && numel(unique(xTar)) > 1, [f1,xi1]=ksdensity(xTar); plot(xi1,f1,'Color',COL.Tar,'LineWidth',1.4); end
            if numel(xRef) >= 7 && numel(unique(xRef)) > 1, [f2,xi2]=ksdensity(xRef); plot(xi2,f2,'Color',COL.Ref,'LineWidth',1.4); end
            title(sprintf('AshmanD=%.2f', met.AshmanD), 'Interpreter','none');
            box off; set(gca,'YAxisLocation','right');
            hold off

            % 3) meanDiff box
            subplot(numel(markers),6,row0+3); cla; hold on
            if doStats
                MakeSpreadAndBoxPlot3_SB({xRef, xTar}, ColsRT, [1 2], LegRT, ...
                    'newfig',0,'paired',0,'optiontest','ranksum','showpoints',1,'ShowSigstar','sig');
                title(sprintf('meanDiff (Tar-Ref)=%.3g', met.meanDiff), 'Interpreter','none');
            else
                axis off; title('meanDiff skipped', 'Interpreter','none');
            end
            box off

            % 4) CohenD box
            subplot(numel(markers),6,row0+4); cla; hold on
            if doStats
                MakeSpreadAndBoxPlot3_SB({xRef, xTar}, ColsRT, [1 2], LegRT, ...
                    'newfig',0,'paired',0,'optiontest','ranksum','showpoints',1,'ShowSigstar','sig');
                title(sprintf('Cohen''s d=%.2f', met.CohensD), 'Interpreter','none');
            else
                axis off; title('CohenD skipped', 'Interpreter','none');
            end
            box off

            % 5) AUROC box
            subplot(numel(markers),6,row0+5); cla; hold on
            if doStats
                MakeSpreadAndBoxPlot3_SB({xRef, xTar}, ColsRT, [1 2], LegRT, ...
                    'newfig',0,'paired',0,'optiontest','ranksum','showpoints',1,'ShowSigstar','sig');
                title(sprintf('AUROC=%.2f', met.AUROC), 'Interpreter','none');
            else
                axis off; title('AUROC skipped', 'Interpreter','none');
            end
            box off

            % 6) ROC
            subplot(numel(markers),6,row0+6); cla; hold on
            plot([0 1],[0 1],':','LineWidth',1.2,'Color',COL.DiagROC);
            if doStats && ~isempty(met.rocFPR)
                plot(met.rocFPR, met.rocTPR, 'k', 'LineWidth',1.6);
                title(sprintf('AUC=%.2f', met.AUROC), 'Interpreter','none');
            else
                title('ROC skipped', 'Interpreter','none');
            end
            xlim([0 1]); ylim([0 1]); axis square
            xlabel('FPR'); ylabel('TPR'); box off
            hold off

            % store row
            rows(end+1,:) = {bodypart, mk, subsetName, wName, ...
                numel(xTar), numel(xRef), ...
                met.AshmanD, met.AUROC, met.CohensD, met.meanDiff, ...
                met.meanRef, met.meanTar, met.stdRef, met.stdTar, ...
                {met.rocFPR}, {met.rocTPR}, ...
                {idxTar},{idxRef},{xTar},{xRef}};
        end

        if isempty(rows)
            close(f);
            continue
        end

        T_full = cell2table(rows,'VariableNames',varNames);

        % Save condition-specific files (like behaviour_analysis)
        baseName = sprintf('%s_%s_%s_%s_sw%.3f', sessname, bodypart, subsetName, wName, opts.smoothing_win_s);
        saveas(f, fullfile(outBase, [baseName '.png']));
        saveas(f, fullfile(outBase, [baseName '.svg']));
        close(f);

        S = struct();
        S.table = T_full;
        S.session = sessname;
        S.bodypart = bodypart;
        S.subset = subsetName;
        S.window = wName;
        save(fullfile(outBase,[baseName '_metrics.mat']), '-struct','S');

        % CSV-safe (drop ROC + raw scalars)
        T_csv = T_full;
        drop = intersect({'rocFPR','rocTPR','trialIdxTar','trialIdxRef','xTar','xRef'}, T_csv.Properties.VariableNames);
        if ~isempty(drop), T_csv = removevars(T_csv, drop); end
        writetable(T_csv, fullfile(outBase,[baseName '_metrics.csv']));

        % Session-level overwrite pointers (one per session)
        save(fullfile(outBase, sprintf('%s_metrics_full.mat', sessname)), 'T_full');
        writetable(T_csv, fullfile(outBase, sprintf('%s_metrics.csv', sessname)));

        OUT.saved{end+1} = baseName; %#ok<AGROW>
    end
end

OUT.outBase = outBase;
end

% ================= helpers (same as before) =================
function [useMid, br] = ra_band_range_from_marker(mk, opts)
useMid = false;
mk = lower(string(mk));
if contains(mk,'delta')
    br = opts.ob_band_delta;
elseif contains(mk,'theta')
    br = opts.ob_band_theta;
elseif contains(mk,'highgamma') || contains(mk,'high_gamma')
    br = opts.ob_band_highgamma; useMid = true;
elseif contains(mk,'gamma')
    br = opts.ob_band_gamma; useMid = true;
else
    error('Unknown marker: %s', mk);
end
end

function hLine = ra_plot_mean_sem_col(t, m, s, col)
t = t(:); m = m(:); s = s(:);
ok = isfinite(t) & isfinite(m) & isfinite(s);
t = t(ok); m = m(ok); s = s(ok);
if numel(t) < 3
    hLine = gobjects(1);
    return
end
x = [t; flipud(t)];
y = [m-s; flipud(m+s)];
hFill = fill(x,y,col);
set(hFill,'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
hLine = plot(t,m,'LineWidth',2,'Color',col);
end

function [winNames, w0, w1] = ra_windows_from_EvRel(EvRel)
winNames = {'stimon_to_stimoff','stimon_to_arrival','stimoff_to_arrival','arrival_m500_to_arrival'};
n = numel(EvRel.stimOff_rel);
w0 = nan(n,numel(winNames));
w1 = nan(n,numel(winNames));
w0(:,1) = 0;                w1(:,1) = EvRel.stimOff_rel;
w0(:,2) = 0;                w1(:,2) = EvRel.arr_rel;
w0(:,3) = EvRel.stimOff_rel;w1(:,3) = EvRel.arr_rel;
w0(:,4) = EvRel.arr_rel-0.5;w1(:,4) = EvRel.arr_rel;
end

function EvRel = ra_get_trial_rel_times(Baphy)
Ttr = Baphy.trial;
n = Baphy.n_trials;
EvRel = struct();
EvRel.stimOff_rel = nan(n,1);
EvRel.arr_rel     = nan(n,1);
EvRel.reward_rel  = nan(n,1);
EvRel.stop_rel    = nan(n,1);

if isfield(Ttr,'rel_stim_start') && isfield(Ttr,'rel_stim_stop')
    EvRel.stimOff_rel = Ttr.rel_stim_stop(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_stim_stop')
    EvRel.stimOff_rel = Ttr.rel_stim_stop(:);
end

if isfield(Ttr,'spout_arrival_rel_s') && isfield(Ttr,'rel_stim_start')
    EvRel.arr_rel = Ttr.spout_arrival_rel_s(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'spout_arrival_rel_s')
    EvRel.arr_rel = Ttr.spout_arrival_rel_s(:);
end

if isfield(Ttr,'rel_reward') && isfield(Ttr,'rel_stim_start')
    EvRel.reward_rel = Ttr.rel_reward(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_reward')
    EvRel.reward_rel = Ttr.rel_reward(:);
end

if isfield(Ttr,'rel_trialend') && isfield(Ttr,'rel_stim_start')
    EvRel.stop_rel = Ttr.rel_trialend(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_trialend')
    EvRel.stop_rel = Ttr.rel_trialend(:);
end
end

function Xmat = ra_bandtrace_from_spec_to_grid(V, t_abs, f_plot, ev_all, t_grid, bmask_t_grid, bandRange)
nEv = numel(ev_all);
Xmat = nan(nEv, numel(t_grid));
bmask = (f_plot >= bandRange(1)) & (f_plot <= bandRange(2));
for k = 1:nEv
    tq = ev_all(k) + t_grid(:);
    if tq(1) < t_abs(1) || tq(end) > t_abs(end), continue, end
    Pk = interp1(t_abs, V, tq, 'linear', NaN);
    basek = nanmean(Pk(bmask_t_grid,:),1);
    Pkb = Pk - basek;
    Xmat(k,:) = nanmean(Pkb(:, bmask), 2);
end
end

function x = ra_trial_scalar_in_window(B, t_rel, t1, t2)
x = nan(size(B,1),1);
for k = 1:size(B,1)
    if ~isfinite(t1(k)) || ~isfinite(t2(k)) || t2(k) <= t1(k), continue, end
    m = (t_rel >= t1(k)) & (t_rel <= t2(k));
    if sum(m) < 5, continue, end
    x(k) = mean(B(k,m), 'omitnan');
end
x = x(isfinite(x));
end

function [t, f, Praw] = ra_extract_spectro_triplet(Spectro)
t = []; f = []; Praw = [];
if isstruct(Spectro)
    if isfield(Spectro,'t'), t = Spectro.t; end
    if isempty(t) && isfield(Spectro,'time'), t = Spectro.time; end
    if isfield(Spectro,'f'), f = Spectro.f; end
    if isempty(f) && isfield(Spectro,'freq'), f = Spectro.freq; end
    if isfield(Spectro,'P'), Praw = Spectro.P; end
    if isempty(Praw) && isfield(Spectro,'S'), Praw = Spectro.S; end
    if isempty(t) || isempty(f) || isempty(Praw)
        error('Spectro struct: missing t/f/P (or S).');
    end
    t = double(t(:)); f = double(f(:)); Praw = double(Praw);
    return
end
if iscell(Spectro)
    A = double(Spectro{1}); B = double(Spectro{2}); C = double(Spectro{3});
    isMat = [~isvector(A), ~isvector(B), ~isvector(C)];
    if sum(isMat) ~= 1, error('Spectro cell: need 1 matrix + 2 vectors.'); end
    if isMat(1), P = A; v1 = B(:); v2 = C(:);
    elseif isMat(2), P = B; v1 = A(:); v2 = C(:);
    else, P = C; v1 = A(:); v2 = B(:);
    end
    [nR,nC] = size(P);
    if numel(v1)==nR && numel(v2)==nC
        t = v1; f = v2; Praw = P;
    elseif numel(v1)==nC && numel(v2)==nR
        t = v2; f = v1; Praw = P;
    else
        error('Spectro cell dims mismatch.');
    end
    return
end
error('Unknown Spectro type');
end

function [t_s, P_tf] = ra_standardize_spectro_time_and_orient(t, f, Praw)
nT = numel(t); nF = numel(f);
if isequal(size(Praw), [nT nF])
    P = Praw;
elseif isequal(size(Praw), [nF nT])
    P = Praw.';
else
    if size(Praw,1)==nT, P=Praw;
    elseif size(Praw,2)==nT, P=Praw.';
    else, error('Spectro dims mismatch'); end
end
good = isfinite(t);
t = t(good); P = P(good,:);
[t, ord] = sort(t); P = P(ord,:);
[t, ia] = unique(t,'stable'); P = P(ia,:);
rowok = all(isfinite(P),2);
t = t(rowok); P = P(rowok,:);
dt0 = median(diff(t));
if dt0 > 5
    t_s = t/1e4;
else
    t_s = t;
end
P_tf = P;
end

function M = ra_effect_metrics_simple(xTar, xRef)
xTar = double(xTar(:)); xRef = double(xRef(:));
xTar = xTar(isfinite(xTar)); xRef = xRef(isfinite(xRef));
M = struct('AUROC',nan,'CohensD',nan,'meanDiff',nan,'AshmanD',nan,'rocFPR',[],'rocTPR',[]);
if numel(xTar) < 3 || numel(xRef) < 3, return, end
M.meanDiff = mean(xTar) - mean(xRef);
m1 = mean(xTar); m0 = mean(xRef);
s1 = var(xTar,1); s0 = var(xRef,1);
sp = sqrt(0.5*(s1+s0));
if sp > 0, M.CohensD = (m1 - m0) / sp; else, M.CohensD = 0; end
if s1 > 0 && s0 > 0
    M.AshmanD = sqrt(2) * abs(m1 - m0) / sqrt(s1 + s0);
else
    M.AshmanD = 0;
end
y = [zeros(numel(xRef),1); ones(numel(xTar),1)];
s = [xRef; xTar];
[~,~,r] = unique(s); r = double(r);
n0 = numel(xRef); n1 = numel(xTar);
R1 = sum(r(n0+1:end));
U1 = R1 - n1*(n1+1)/2;
M.AUROC = U1 / (n0*n1);

thr = unique(s); thr = sort(thr,'ascend');
FPR = nan(numel(thr)+2,1); TPR = nan(numel(thr)+2,1);
FPR(1)=0; TPR(1)=0;
for i = 1:numel(thr)
    t = thr(i);
    yhat = s >= t;
    tp = sum(yhat==1 & y==1);
    fp = sum(yhat==1 & y==0);
    fn = sum(yhat==0 & y==1);
    tn = sum(yhat==0 & y==0);
    TPR(i+1) = tp / max(1,(tp+fn));
    FPR(i+1) = fp / max(1,(fp+tn));
end
FPR(end)=1; TPR(end)=1;
[FF,ord] = sort(FPR);
M.rocFPR = FF;
M.rocTPR = TPR(ord);
end
