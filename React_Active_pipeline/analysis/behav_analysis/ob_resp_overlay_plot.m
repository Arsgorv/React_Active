function Out = ob_resp_overlay_plot(datapath, SpectroLow, Baphy, stim_on_abs_s, G, opts, ttl)
% ob_resp_overlay_plot
% StimOn-aligned spectrograms (LOW + MID) with respiration/sniff piezo + OB band overlays + coupling metrics.
%
% OUTPUTS:
% - Per-group (all/tar/ref/nosound/nomotor) LOW+MID mean spectrograms (baseline-corrected) + overlays.
% - Coupling: trial-wise xcorr peak lag/r in 3 canonical windows for delta/theta/gamma/highgamma, PER GROUP.
% - Exports:
%   * summaryTable: session x piezoMode x group x spec x band x window x {nUsed, medianLag_s, iqrLag_s, meanR, medianR}
%   * trialTable:   session x piezoMode x group x spec x band x window x trialIdx x trialType x {peakLag_s, peakR}
% - Optional behaviour-style figure for OB power classification (Tar vs Ref) in a chosen window (opts.make_behaviour_style).
%
% PIEZO MODES:
% opts.piezo_mode = 'resp' | 'sniff' | 'both'  (DEFAULT = 'both')
% If 'both', the whole pipeline is duplicated: separate figures + separate outputs per mode.

if nargin < 6 || isempty(opts), opts = struct(); end
if nargin < 7, ttl = ''; end

[~, sessname] = fileparts(datapath);

% canonical windows (trial-relative)
winDefs = {'stimOn_to_stimOff','stimOff_to_arr','arr_to_stop'};

% ---------------- load piezo channel ----------------
respFile = fullfile(datapath,'ephys','LFPData', sprintf('LFP%d.mat', opts.resp_lfp_chan));
if ~exist(respFile,'file')
    error('Piezo LFP file not found: %s', respFile);
end
S = load(respFile,'LFP');
Resp = S.LFP;

% ---------------- events / rel times ----------------
EvRel = ra_get_trial_rel_times(Baphy);

nTr = Baphy.n_trials;
stim_on_abs_s = stim_on_abs_s(:);
if numel(stim_on_abs_s) ~= nTr
    error('stim_on_abs_s must be n_trials x 1');
end

% ---------------- sampling rate estimate ----------------
rt = double(Range(Resp,'s'));
fs = round(1/median(diff(rt)));
if ~isfinite(fs) || fs <= 0, fs = 1024; end

% ---------------- compute piezo bands ----------------
[respSig, sniffSig] = ra_piezo_two_bands(Resp, fs, opts);

% ---------------- load MID spectrogram automatically ----------------
SpectroMid = [];
midFile = fullfile(datapath,'ephys','B_Middle_Spectrum.mat');
if exist(midFile,'file')
    Sm = load(midFile,'Spectro');
    if isfield(Sm,'Spectro'), SpectroMid = Sm.Spectro; end
end

% ---------------- parse LOW spectrogram triplet ----------------
[tL, fL, PLraw] = ra_extract_spectro_triplet(SpectroLow);
[t_abs_L, P_tf_L] = ra_standardize_spectro_time_and_orient(tL, fL, PLraw); % time x freq
P_tf_L = log10(max(double(P_tf_L), eps));

% ---------------- parse MID spectrogram triplet ----------------
t_abs_M = []; fM = []; P_tf_M = [];
hasMID = ~isempty(SpectroMid);
if hasMID
    [tM, fM, PMraw] = ra_extract_spectro_triplet(SpectroMid);
    [t_abs_M, P_tf_M] = ra_standardize_spectro_time_and_orient(tM, fM, PMraw);
    P_tf_M = log10(max(double(P_tf_M), eps));
end

% ---------------- build eventlocked grids ----------------
dtL = median(diff(t_abs_L));
t_rel_L = (opts.spec_twin_whole_s(1):dtL:opts.spec_twin_whole_s(2))';
bmask_t_L = (t_rel_L >= opts.spec_baseline_s(1)) & (t_rel_L <= opts.spec_baseline_s(2));

if hasMID
    dtM = median(diff(t_abs_M));
    t_rel_M = (opts.spec_twin_whole_s(1):dtM:opts.spec_twin_whole_s(2))';
    bmask_t_M = (t_rel_M >= opts.spec_baseline_s(1)) & (t_rel_M <= opts.spec_baseline_s(2));
else
    dtM = []; t_rel_M = []; bmask_t_M = [];
end

% ---------------- restrict frequency axes ----------------
fmaskL = (fL >= opts.spec_freq_xlim_low(1)) & (fL <= opts.spec_freq_xlim_low(2));
f_plot_L = fL(fmaskL);
V_L = P_tf_L(:, fmaskL);

if hasMID
    fmaskM = (fM >= opts.spec_freq_xlim_mid(1)) & (fM <= opts.spec_freq_xlim_mid(2));
    f_plot_M = fM(fmaskM);
    V_M = P_tf_M(:, fmaskM);
else
    f_plot_M = []; V_M = [];
end

% ---------------- choose modes to run ----------------
mode = lower(string(opts.piezo_mode));
if mode == "both"
    modeList = {'resp','sniff'};
elseif mode == "sniff"
    modeList = {'sniff'};
else
    modeList = {'resp'};
end

Out = struct();
Out.opts_snapshot = opts;
Out.groups = G.names;
Out.winDefs = winDefs;
Out.hasMID = hasMID;

tmpl = struct( ...
    'modeName', '', ...
    'piezo_band_hz', [], ...
    'piezo_smooth_s', [], ...
    'low', struct('t_rel',[],'f_plot',[],'group',[]), ...
    'mid', struct('t_rel',[],'f_plot',[],'group',[]), ...
    'coupling', struct('group',[]), ...
    'figure1_handle', [], ...
    'figure2_handle', [], ...
    'figure3_handle', [], ...
    'figure4_handle', [], ...
    'figure5_handle', [], ...
    'figureB_handle', [], ...
    'phase', struct('phCenters',[],'meanDeltaPow',[],'meanThetaPow',[],'nPerBin',[]), ...
    'export', struct('summaryTable',[],'trialTable',[]) ...
);

Out.mode = repmat(tmpl, 1, numel(modeList));

% Precompute window edges per trial (relative to stimOn=0)
[w1_all, w2_all] = ra_get_window_edges_all(winDefs, EvRel); % [nTr x 3]

for mi = 1:numel(modeList)
    modeName = modeList{mi};

    if strcmpi(modeName,'sniff')
        piezoSig = sniffSig;
        piezoBand = opts.sniff_band_hz;
        piezoSmooth = opts.sniff_smooth_s;
    else
        piezoSig = respSig;
        piezoBand = opts.resp_band_hz;
        piezoSmooth = opts.resp_smooth_s;
    end

    t_piezo_abs = piezoSig.t_abs(:);
    env_piezo   = piezoSig.env(:);
    ph_piezo    = piezoSig.ph(:);

    MOut = struct();
    MOut.modeName = modeName;
    MOut.piezo_band_hz = piezoBand;
    MOut.piezo_smooth_s = piezoSmooth;

    % containers
    nG = numel(G.names);

    MOut.low = struct();
    MOut.mid = struct();
    MOut.low.t_rel = t_rel_L; MOut.low.f_plot = f_plot_L;
    if hasMID
        MOut.mid.t_rel = t_rel_M; MOut.mid.f_plot = f_plot_M;
    else
        MOut.mid.t_rel = []; MOut.mid.f_plot = [];
    end

    MOut.low.group = repmat(struct(), nG, 1);
    MOut.mid.group = repmat(struct(), nG, 1);

    % coupling per group x band x window
    bandList = {'delta','theta','gamma','highgamma'};
    specOfBand = {'LOW','LOW','MID','MID'};

    MOut.coupling = struct();
    MOut.coupling.group = repmat(struct(), nG, 1);
    for gi = 1:nG
        MOut.coupling.group(gi).name = string(G.names{gi});
        MOut.coupling.group(gi).band = repmat(struct('name','','spec','','win',[]), numel(bandList), 1);
        for bi = 1:numel(bandList)
            MOut.coupling.group(gi).band(bi).name = string(bandList{bi});
            MOut.coupling.group(gi).band(bi).spec = string(specOfBand{bi});
            MOut.coupling.group(gi).band(bi).win = repmat(struct('name','','trial',[],'summary',[]), numel(winDefs), 1);
            for wi = 1:numel(winDefs)
                MOut.coupling.group(gi).band(bi).win(wi).name = string(winDefs{wi});
                MOut.coupling.group(gi).band(bi).win(wi).trial = struct('peakLag_s', nan(nTr,1), 'peakR', nan(nTr,1));
                MOut.coupling.group(gi).band(bi).win(wi).summary = struct('nUsed',0,'medianLag_s',nan,'iqrLag_s',nan,'meanR',nan,'medianR',nan);
            end
        end
    end

    % =========================
    % FIGURE 1: spectrograms + overlays (LOW row, MID row)
    % =========================
    ttl1 = sprintf('%s | %s piezo', ttl, upper(modeName));
    fig1 = figure('Color','w','Units','pixels','Position',[30 60 1900 900], 'Name', ttl1);

    climL = [inf -inf];
    climM = [inf -inf];

    for gi = 1:nG
        mG = G.mask(:,gi) & isfinite(stim_on_abs_s);
        idxG = find(mG);

        if isempty(idxG)
            MOut.low.group(gi).idxUsed = [];
            MOut.low.group(gi).nUsed = 0;
            MOut.low.group(gi).Lb_mean = nan(numel(t_rel_L), numel(f_plot_L));
            MOut.low.group(gi).piezo_z = nan(numel(t_rel_L),1);
            MOut.low.group(gi).delta_z = nan(numel(t_rel_L),1);
            MOut.low.group(gi).theta_z = nan(numel(t_rel_L),1);
            if hasMID
                MOut.mid.group(gi).idxUsed = [];
                MOut.mid.group(gi).nUsed = 0;
                MOut.mid.group(gi).Lb_mean = nan(numel(t_rel_M), numel(f_plot_M));
                MOut.mid.group(gi).piezo_z = nan(numel(t_rel_M),1);
                MOut.mid.group(gi).gamma_z = nan(numel(t_rel_M),1);
                MOut.mid.group(gi).highgamma_z = nan(numel(t_rel_M),1);
            end
            continue
        end

        ev = stim_on_abs_s(idxG);

        % LOW mean spec
        [Lb_mean_L, usedMask_L] = ra_eventlocked_mean_spec(V_L, t_abs_L, ev, t_rel_L, bmask_t_L);
        idxUsed_L = idxG(usedMask_L);
        MOut.low.group(gi).idxUsed = idxUsed_L;
        MOut.low.group(gi).nUsed = numel(idxUsed_L);
        MOut.low.group(gi).Lb_mean = Lb_mean_L;

        if any(isfinite(Lb_mean_L(:)))
            climL(1) = min(climL(1), nanmin(Lb_mean_L(:)));
            climL(2) = max(climL(2), nanmax(Lb_mean_L(:)));
        end

        [piezo_z_L, delta_z_L, theta_z_L] = ra_group_overlay_traces_low_onepiezo( ...
            V_L, t_abs_L, f_plot_L, ev(usedMask_L), t_rel_L, bmask_t_L, ...
            t_piezo_abs, env_piezo, opts);
        MOut.low.group(gi).piezo_z = piezo_z_L;
        MOut.low.group(gi).delta_z = delta_z_L;
        MOut.low.group(gi).theta_z = theta_z_L;

        % MID mean spec
        if hasMID
            [Lb_mean_M, usedMask_M] = ra_eventlocked_mean_spec(V_M, t_abs_M, ev, t_rel_M, bmask_t_M);
            idxUsed_M = idxG(usedMask_M);
            MOut.mid.group(gi).idxUsed = idxUsed_M;
            MOut.mid.group(gi).nUsed = numel(idxUsed_M);
            MOut.mid.group(gi).Lb_mean = Lb_mean_M;

            if any(isfinite(Lb_mean_M(:)))
                climM(1) = min(climM(1), nanmin(Lb_mean_M(:)));
                climM(2) = max(climM(2), nanmax(Lb_mean_M(:)));
            end

            [piezo_z_M, gamma_z_M, highgamma_z_M] = ra_group_overlay_traces_mid_onepiezo( ...
                V_M, t_abs_M, f_plot_M, ev(usedMask_M), t_rel_M, bmask_t_M, ...
                t_piezo_abs, env_piezo, opts);
            MOut.mid.group(gi).piezo_z     = piezo_z_M;
            MOut.mid.group(gi).gamma_z     = gamma_z_M;
            MOut.mid.group(gi).highgamma_z = highgamma_z_M;
        end
    end

    if ~all(isfinite(climL)) || climL(2) <= climL(1), climL = []; end
    if ~all(isfinite(climM)) || climM(2) <= climM(1), climM = []; end

    % plot LOW row
    for gi = 1:nG
        ax = subplot(2, nG, gi);
        imagesc(t_rel_L, f_plot_L, MOut.low.group(gi).Lb_mean'); axis xy
        box off
        try colormap(ax,'viridis'); catch, colormap(ax,parula); end
        if ~isempty(climL), caxis(ax, climL); end
        xlabel('t from stimOn (s)'); ylabel('Hz');

        idxUsed = MOut.low.group(gi).idxUsed;
        [t_stimOff, t_arr, t_reward] = ra_group_event_lines(EvRel, idxUsed);
        xline(ax,0,'k-','LineWidth',1);
        if isfinite(t_stimOff), xline(ax,t_stimOff,'k--','LineWidth',1); end
        if isfinite(t_arr),     xline(ax,t_arr,'k--','LineWidth',1); end
        if isfinite(t_reward),  xline(ax,t_reward,'k--','LineWidth',1); end

        title(sprintf('%s LOW (n=%d)', G.names{gi}, MOut.low.group(gi).nUsed), 'Interpreter','none');

        yyaxis(ax,'right'); hold(ax,'on');
        hP = plot(ax, t_rel_L, runmean(MOut.low.group(gi).piezo_z, 5), ...
            'LineStyle','-','Marker','none','LineWidth',1.5,'Color','w');
        hD = plot(ax, t_rel_L, runmean(MOut.low.group(gi).delta_z, 5), ...
            'LineStyle','-','Marker','none','LineWidth',1.5,'Color',[0.75 0.65 0.20]);
        hT = plot(ax, t_rel_L, runmean(MOut.low.group(gi).theta_z, 5), ...
            'LineStyle','-','Marker','none','LineWidth',1.5,'Color',[0.55 0.15 0.15]);

        ylabel(ax,'z'); ylim(ax,[-3 3]);
        if gi == 1
            legend(ax, [hP hD hT], {sprintf('%s env', modeName),'delta','theta'}, ...
                'Location','northwest','Box','off');
        end
        hold(ax,'off');
        yyaxis(ax,'left');
    end

    % plot MID row
    if hasMID
        for gi = 1:nG
            ax = subplot(2, nG, nG + gi);
            imagesc(t_rel_M, f_plot_M, MOut.mid.group(gi).Lb_mean'); axis xy
            box off
            try colormap(ax,'viridis'); catch, colormap(ax,parula); end
            if ~isempty(climM), caxis(ax, climM); end
            xlabel('t from stimOn (s)'); ylabel('Hz');

            idxUsed = MOut.mid.group(gi).idxUsed;
            [t_stimOff, t_arr, t_reward] = ra_group_event_lines(EvRel, idxUsed);
            xline(ax,0,'k-','LineWidth',1);
            if isfinite(t_stimOff), xline(ax,t_stimOff,'k--','LineWidth',1); end
            if isfinite(t_arr),     xline(ax,t_arr,'k--','LineWidth',1); end
            if isfinite(t_reward),  xline(ax,t_reward,'k--','LineWidth',1); end

            title(sprintf('%s MID (n=%d)', G.names{gi}, MOut.mid.group(gi).nUsed), 'Interpreter','none');

            yyaxis(ax,'right'); hold(ax,'on');
            hP = plot(ax, t_rel_M, runmean(MOut.mid.group(gi).piezo_z, 75), 'LineStyle','-', 'Marker','none','LineWidth',1.5, 'Color','w');
            hG = plot(ax, t_rel_M, runmean(MOut.mid.group(gi).gamma_z, 75), 'LineStyle','-', 'Marker','none','LineWidth',1.5, 'Color',[0.75 0.65 0.20]);
            hH = plot(ax, t_rel_M, runmean(MOut.mid.group(gi).highgamma_z, 75), 'LineStyle','-', 'Marker','none','LineWidth',1.5, 'Color',[0.55 0.15 0.15]);

            ylabel(ax,'z'); ylim(ax,[-3 3]);
            if gi == 1
                legend(ax, [hP hG hH], {sprintf('%s env', modeName),'gamma','high gamma'}, ...
                    'Location','northwest','Box','off');
            end
            hold(ax,'off');
            yyaxis(ax,'left');
        end
    else
        subplot(2,1,2); axis off
        text(0.1,0.5,'No B_Middle_Spectrum.mat found -> MID row skipped','Interpreter','none');
    end

    sgtitle(ttl1,'Interpreter','none');
    MOut.figure1_handle = fig1;

    % =========================
    % Coupling metrics PER GROUP (trial-wise xcorr peaks) + exports
    % =========================
    t_grid = t_rel_L(:);
    maxLagSamp = round(opts.xcorr_maxLag_s / dtL);
    lags_s = (-maxLagSamp:maxLagSamp) * dtL;

    for gi = 1:nG
        idxG = find(G.mask(:,gi) & isfinite(stim_on_abs_s));
        if isempty(idxG), continue, end

        % piezo aligned matrix for this group (rows correspond to idxG)
        PiezoMat = nan(numel(idxG), numel(t_grid));
        for k = 1:numel(idxG)
            ii = idxG(k);
            tq = stim_on_abs_s(ii) + t_grid;
            PiezoMat(k,:) = interp1(t_piezo_abs, env_piezo, tq, 'linear', NaN);
        end

        % band traces onto LOW grid (use LOW baseline mask)
        DeltaMat = ra_bandtrace_from_spec_to_grid(V_L, t_abs_L, f_plot_L, stim_on_abs_s(idxG), t_grid, bmask_t_L, opts.ob_band_delta);
        ThetaMat = ra_bandtrace_from_spec_to_grid(V_L, t_abs_L, f_plot_L, stim_on_abs_s(idxG), t_grid, bmask_t_L, opts.ob_band_theta);

        if hasMID
            % for MID, use MID baseline mask but still map to LOW time grid (OK; baseline uses bmask_t_M with same window)
            GammaMat = ra_bandtrace_from_spec_to_grid(V_M, t_abs_M, f_plot_M, stim_on_abs_s(idxG), t_grid, (t_grid>=opts.spec_baseline_s(1) & t_grid<=opts.spec_baseline_s(2)), opts.ob_band_gamma);
            HGMat    = ra_bandtrace_from_spec_to_grid(V_M, t_abs_M, f_plot_M, stim_on_abs_s(idxG), t_grid, (t_grid>=opts.spec_baseline_s(1) & t_grid<=opts.spec_baseline_s(2)), opts.ob_band_highgamma);
        else
            GammaMat = nan(size(PiezoMat));
            HGMat    = nan(size(PiezoMat));
        end

        mats = {DeltaMat, ThetaMat, GammaMat, HGMat};

        % compute per-window peaks, store into nTr-sized vectors (by original trialIdx)
        for bi = 1:numel(bandList)
            Xband = mats{bi};
            for wi = 1:numel(winDefs)
                pkLag = nan(nTr,1);
                pkR   = nan(nTr,1);

                for k = 1:numel(idxG)
                    ii = idxG(k);
                    t1 = w1_all(ii,wi); t2 = w2_all(ii,wi);
                    if ~isfinite(t1) || ~isfinite(t2) || t2 <= t1, continue, end

                    m = (t_grid >= t1) & (t_grid <= t2);
                    if sum(m) < 10, continue, end

                    p = ra_zscore_safe(PiezoMat(k,m));
                    b = ra_zscore_safe(Xband(k,m));
                    good = isfinite(p) & isfinite(b);
                    if sum(good) < 10, continue, end

                    [cc, lags] = xcorr(b(good), p(good), maxLagSamp, 'coeff');
                    [pk, im] = max(cc);

                    pkR(ii) = pk;
                    pkLag(ii) = lags(im) * dtL;
                end

                MOut.coupling.group(gi).band(bi).win(wi).trial.peakLag_s = pkLag;
                MOut.coupling.group(gi).band(bi).win(wi).trial.peakR     = pkR;

                ok = isfinite(pkLag) & isfinite(pkR);
                lagv = pkLag(ok); rv = pkR(ok);

                Sx = struct();
                Sx.nUsed = sum(ok);
                Sx.medianLag_s = median(lagv,'omitnan');
                Sx.iqrLag_s = iqr(lagv);
                Sx.meanR = mean(rv,'omitnan');
                Sx.medianR = median(rv,'omitnan');

                MOut.coupling.group(gi).band(bi).win(wi).summary = Sx;
            end
        end
    end

    % =========================
    % FIGURE 2: coupling summary for ALL group + overlay (ALL)
    % =========================
    fig2 = figure('Color','w','Units','pixels','Position',[80 80 1700 900], 'Name', [ttl1 ' | coupling']);
    gi_all = find(strcmpi(string(G.names),'all'),1,'first');
    if isempty(gi_all), gi_all = 1; end

    plotBands = {'delta','theta'};
    for bi2 = 1:numel(plotBands)
        bname = plotBands{bi2};
        bix = find(strcmpi(cellstr(string(bandList)), bname), 1, 'first');
        if isempty(bix), continue, end

        allLag = []; allR = [];
        for wi = 1:numel(winDefs)
            L = MOut.coupling.group(gi_all).band(bix).win(wi).trial.peakLag_s;
            R = MOut.coupling.group(gi_all).band(bix).win(wi).trial.peakR;
            ok = isfinite(L) & isfinite(R);
            allLag = [allLag; L(ok)]; %#ok<AGROW>
            allR   = [allR; R(ok)]; %#ok<AGROW>
        end

        r0 = (bi2-1)*4;

        subplot(3,4,r0+1);
        histogram(allLag, 30);
        xlabel('peak lag (s)'); ylabel('count');
        title(sprintf('%s: peak lag (all windows)', bname), 'Interpreter','none'); box off

        subplot(3,4,r0+2);
        histogram(allR, 30);
        xlabel('peak r'); ylabel('count');
        title(sprintf('%s: peak r (all windows)', bname), 'Interpreter','none'); box off

        subplot(3,4,r0+3);
        plot(allLag, allR, '.', 'MarkerSize', 8);
        xlabel('peak lag (s)'); ylabel('peak r');
        title(sprintf('%s: lag vs r', bname), 'Interpreter','none'); box off
        xline(0,'k:'); yline(0,'k:');

        subplot(3,4,r0+4);
        hold on
        for wi = 1:numel(winDefs)
            Sx = MOut.coupling.group(gi_all).band(bix).win(wi).summary;
            plot(wi, Sx.meanR, 'o', 'MarkerSize', 6);
            txt2 = sprintf('n=%d%smedLag=%.2f', Sx.nUsed, char(10), Sx.medianLag_s);
            text(wi+0.05, Sx.meanR, txt2, 'FontSize', 8, 'Interpreter','none');
        end
        xlim([0.5 numel(winDefs)+0.5]);
        set(gca,'XTick',1:numel(winDefs),'XTickLabel',winDefs,'XTickLabelRotation',20);
        ylabel('mean r');
        title(sprintf('%s: window summaries (ALL)', bname), 'Interpreter','none');
        box off
        hold off
    end

    subplot(3,4,9:12);
    piezo_z = MOut.low.group(gi_all).piezo_z(:);
    delta_z = MOut.low.group(gi_all).delta_z(:);
    theta_z = MOut.low.group(gi_all).theta_z(:);

    hold on
    plot(t_rel_L, piezo_z,  'LineWidth', 1);
    plot(t_rel_L, delta_z,  'LineWidth', 1);
    plot(t_rel_L, theta_z,  'LineWidth', 1);
    xline(0,'k-'); yline(0,'k:');

    idxUsedAll = MOut.low.group(gi_all).idxUsed;
    [t_stimOff, t_arr, t_reward] = ra_group_event_lines(EvRel, idxUsedAll);
    if isfinite(t_stimOff), xline(t_stimOff,'k--'); end
    if isfinite(t_arr),     xline(t_arr,'k--'); end
    if isfinite(t_reward),  xline(t_reward,'k--'); end

    xlabel('t from stimOn (s)'); ylabel('z');
    title(sprintf('Overlay (ALL): %s env, delta, theta', modeName), 'Interpreter','none');
    legend({sprintf('%s env', modeName),'delta','theta'},'Location','northwest','Box','off');
    box off
    hold off

    sgtitle([ttl1 ' | piezo–OB coupling'], 'Interpreter','none');
    MOut.figure2_handle = fig2;

    % =========================
    % FIGURE 3/4: trial-wise and windowed xcorr heatmaps (ALL only, diagnostic)
    % =========================
    use_win = opts.xcorr_use_win;
    wmask = (t_grid >= use_win(1)) & (t_grid <= use_win(2));

    idxAll = find(G.mask(:,gi_all) & isfinite(stim_on_abs_s));
    PiezoAll = nan(numel(idxAll), numel(t_grid));
    for k = 1:numel(idxAll)
        ii = idxAll(k);
        tq = stim_on_abs_s(ii) + t_grid;
        PiezoAll(k,:) = interp1(t_piezo_abs, env_piezo, tq, 'linear', NaN);
    end
    DeltaAll = ra_bandtrace_from_spec_to_grid(V_L, t_abs_L, f_plot_L, stim_on_abs_s(idxAll), t_grid, bmask_t_L, opts.ob_band_delta);
    ThetaAll = ra_bandtrace_from_spec_to_grid(V_L, t_abs_L, f_plot_L, stim_on_abs_s(idxAll), t_grid, bmask_t_L, opts.ob_band_theta);

    [XC_D, pkLagD_samp, ~] = ra_trialwise_xcorr_mat(DeltaAll, PiezoAll, wmask, maxLagSamp);
    [XC_T, pkLagT_samp, ~] = ra_trialwise_xcorr_mat(ThetaAll, PiezoAll, wmask, maxLagSamp);
    pkLagD = pkLagD_samp * dtL;
    pkLagT = pkLagT_samp * dtL;

    fig3 = figure('Color','w','Units','pixels','Position',[120 80 1400 700], 'Name', [ttl1 ' | trialwise xcorr']);
    subplot(1,2,1);
    [~,ord] = sort(pkLagD,'ascend');
    imagesc(lags_s, 1:size(XC_D,1), XC_D(ord,:)); axis xy
    try colormap('viridis'); catch, colormap(parula); end
    caxis([-0.5 0.5]); colorbar
    xlabel('Lag (s): band leads (+) piezo'); ylabel('Trials (sorted by peak lag)');
    title(sprintf('delta vs %s: trial-wise xcorr (ALL)', modeName), 'Interpreter','none'); box off

    subplot(1,2,2);
    [~,ord] = sort(pkLagT,'ascend');
    imagesc(lags_s, 1:size(XC_T,1), XC_T(ord,:)); axis xy
    try colormap('viridis'); catch, colormap(parula); end
    caxis([-0.5 0.5]); colorbar
    xlabel('Lag (s): band leads (+) piezo'); ylabel('Trials (sorted by peak lag)');
    title(sprintf('theta vs %s: trial-wise xcorr (ALL)', modeName), 'Interpreter','none'); box off
    MOut.figure3_handle = fig3;

    % windowed mean xcorr
    win_s  = 1.0;
    step_s = 0.10;
    tCenters = (use_win(1):step_s:use_win(2))';
    half = win_s/2;

    XCwD = nan(numel(tCenters), numel(lags_s));
    XCwT = nan(numel(tCenters), numel(lags_s));

    for it = 1:numel(tCenters)
        c0 = tCenters(it);
        m = (t_grid >= (c0-half)) & (t_grid <= (c0+half));
        if sum(m) < 10, continue, end

        XD = []; XT = [];
        for k = 1:size(PiezoAll,1)
            p = ra_zscore_safe(PiezoAll(k,m));
            d = ra_zscore_safe(DeltaAll(k,m));
            th = ra_zscore_safe(ThetaAll(k,m));

            okD = isfinite(p) & isfinite(d);
            okT = isfinite(p) & isfinite(th);

            if sum(okD) >= 10
                cc = xcorr(d(okD)-mean(d(okD)), p(okD)-mean(p(okD)), maxLagSamp, 'coeff');
                XD = [XD; cc(:)']; %#ok<AGROW>
            end
            if sum(okT) >= 10
                cc = xcorr(th(okT)-mean(th(okT)), p(okT)-mean(p(okT)), maxLagSamp, 'coeff');
                XT = [XT; cc(:)']; %#ok<AGROW>
            end
        end

        if ~isempty(XD), XCwD(it,:) = mean(XD,1,'omitnan'); end
        if ~isempty(XT), XCwT(it,:) = mean(XT,1,'omitnan'); end
    end

    fig4 = figure('Color','w','Units','pixels','Position',[160 80 1400 700], 'Name', [ttl1 ' | windowed xcorr']);
    subplot(1,2,1);
    imagesc(lags_s, tCenters, XCwD); axis xy
    try colormap('viridis'); catch, colormap(parula); end
    caxis([-1 1]); colorbar
    xlabel('Lag (s): band leads (+) piezo'); ylabel('Time from stimOn (s)');
    title(sprintf('delta vs %s: windowed mean xcorr (ALL)', modeName), 'Interpreter','none'); box off

    subplot(1,2,2);
    imagesc(lags_s, tCenters, XCwT); axis xy
    try colormap('viridis'); catch, colormap(parula); end
    caxis([-0.75 0.75]); colorbar
    xlabel('Lag (s): band leads (+) piezo'); ylabel('Time from stimOn (s)');
    title(sprintf('theta vs %s: windowed mean xcorr (ALL)', modeName), 'Interpreter','none'); box off
    MOut.figure4_handle = fig4;

    % =========================
    % Phase -> power (delta/theta vs piezo phase), ALL group
    % =========================
    PiezoPhaseMat = nan(numel(idxAll), numel(t_grid));
    for k = 1:numel(idxAll)
        ii = idxAll(k);
        tq = stim_on_abs_s(ii) + t_grid;
        PiezoPhaseMat(k,:) = ra_interp_phase_circ(t_piezo_abs, ph_piezo, tq);
    end

    bmask_grid = (t_grid >= opts.spec_baseline_s(1)) & (t_grid <= opts.spec_baseline_s(2));
    DeltaBL = ra_subtract_baseline_per_trial(DeltaAll, bmask_grid);
    ThetaBL = ra_subtract_baseline_per_trial(ThetaAll, bmask_grid);

    [phC, dPow, nBin] = ra_phase_to_power(DeltaBL, PiezoPhaseMat, t_grid, opts.phase_win_s, opts.phase_nBins);
    [~,  tPow, ~]     = ra_phase_to_power(ThetaBL, PiezoPhaseMat, t_grid, opts.phase_win_s, opts.phase_nBins);

    fig5 = figure('Color','w','Units','pixels','Position',[250 120 1200 450], 'Name', [ttl1 ' | phase->power']);
    subplot(1,2,1);
    yyaxis left; plot(phC, dPow, 'k-','LineWidth',2); ylabel('Mean \delta power (BL-corr)');
    yyaxis right; plot(phC, nBin, '-','LineWidth',1); ylabel('N samples/bin');
    xlabel(sprintf('%s phase (rad)', modeName)); xlim([-pi pi]);
    title(sprintf('%s phase -> \\delta power (ALL)', modeName),'Interpreter','none'); box off

    subplot(1,2,2);
    yyaxis left; plot(phC, tPow, 'k-','LineWidth',2); ylabel('Mean \theta power (BL-corr)');
    yyaxis right; plot(phC, nBin, '-','LineWidth',1); ylabel('N samples/bin');
    xlabel(sprintf('%s phase (rad)', modeName)); xlim([-pi pi]);
    title(sprintf('%s phase -> \\theta power (ALL)', modeName),'Interpreter','none'); box off

    MOut.figure5_handle = fig5;
    MOut.phase.phCenters = phC;
    MOut.phase.meanDeltaPow = dPow;
    MOut.phase.meanThetaPow = tPow;
    MOut.phase.nPerBin = nBin;

    % =========================
    % EXPORT TABLES (group-aware)
    % =========================
    % summary: session x piezoMode x group x spec x band x window x stats
    rowsS = {};
    for gi = 1:nG
        gname = string(G.names{gi});
        for bi = 1:numel(bandList)
            bname = string(bandList{bi});
            spec  = string(specOfBand{bi});
            for wi = 1:numel(winDefs)
                Sx = MOut.coupling.group(gi).band(bi).win(wi).summary;
                rowsS(end+1,:) = {string(sessname), string(modeName), gname, spec, bname, string(winDefs{wi}), ...
                    Sx.nUsed, Sx.medianLag_s, Sx.iqrLag_s, Sx.meanR, Sx.medianR}; %#ok<AGROW>
            end
        end
    end
    if isempty(rowsS), rowsS = cell(0,11); end
    MOut.export.summaryTable = cell2table(rowsS, 'VariableNames', ...
        {'session','piezoMode','group','spec','band','window','nUsed','medianLag_s','iqrLag_s','meanR','medianR'});

    % trial: session x piezoMode x group x spec x band x window x trialIdx x trialType x peaks
    tt = strings(nTr,1);
    if isfield(Baphy,'trial') && isstruct(Baphy.trial)
        if isfield(Baphy.trial,'trialType')
            try tt = string({Baphy.trial.trialType})'; catch, end %#ok<CTCH>
        elseif isfield(Baphy.trial,'type')
            try tt = string({Baphy.trial.type})'; catch, end %#ok<CTCH>
        end
    end
    if all(tt=="")
        % fallback from G masks: tar/ref
        tt(:) = "other";
        ciT = find(strcmpi(G.names,'tar'),1); ciR = find(strcmpi(G.names,'ref'),1);
        if ~isempty(ciT), tt(G.mask(:,ciT)) = "tar"; end
        if ~isempty(ciR), tt(G.mask(:,ciR)) = "ref"; end
    end

    rowsT = {};
    for gi = 1:nG
        gname = string(G.names{gi});
        for bi = 1:numel(bandList)
            bname = string(bandList{bi});
            spec  = string(specOfBand{bi});
            for wi = 1:numel(winDefs)
                L = MOut.coupling.group(gi).band(bi).win(wi).trial.peakLag_s;
                R = MOut.coupling.group(gi).band(bi).win(wi).trial.peakR;
                for tr = 1:nTr
                    if ~isfinite(L(tr)) || ~isfinite(R(tr)), continue, end
                    rowsT(end+1,:) = {string(sessname), string(modeName), gname, spec, bname, string(winDefs{wi}), ...
                        tr, tt(tr), L(tr), R(tr)}; %#ok<AGROW>
                end
            end
        end
    end
    if isempty(rowsT), rowsT = cell(0,10); end
    MOut.export.trialTable = cell2table(rowsT, 'VariableNames', ...
        {'session','piezoMode','group','spec','band','window','trialIdx','trialType','peakLag_s','peakR'});

    % =========================
    % Optional: behaviour-style figure for OB (Tar vs Ref) in a chosen window
    % =========================
    figB = [];
    if opts.make_behaviour_style
        figB = ra_plot_behaviour_style_ob(datapath, sessname, Baphy, G, stim_on_abs_s, ...
            V_L, t_abs_L, f_plot_L, V_M, t_abs_M, f_plot_M, hasMID, ...
            t_rel_L, bmask_t_L, ...
            winDefs, w1_all, w2_all, opts, ttl, modeName);
    end
    MOut.figureB_handle = figB;

    % =========================
    % finalize Out.mode(mi)
    % =========================
    Out.mode(mi).modeName        = modeName;
    Out.mode(mi).piezo_band_hz   = piezoBand;
    Out.mode(mi).piezo_smooth_s  = piezoSmooth;

    Out.mode(mi).low             = MOut.low;
    Out.mode(mi).mid             = MOut.mid;
    Out.mode(mi).coupling        = MOut.coupling;

    Out.mode(mi).figure1_handle  = MOut.figure1_handle;
    Out.mode(mi).figure2_handle  = MOut.figure2_handle;
    Out.mode(mi).figure3_handle  = MOut.figure3_handle;
    Out.mode(mi).figure4_handle  = MOut.figure4_handle;
    Out.mode(mi).figure5_handle  = MOut.figure5_handle;
    Out.mode(mi).figureB_handle  = MOut.figureB_handle;

    Out.mode(mi).phase           = MOut.phase;
    Out.mode(mi).export          = MOut.export;
end

end

% =========================================================================
% helpers
% =========================================================================

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

function [t_stimOff, t_arr, t_reward] = ra_group_event_lines(EvRel, idxUsed)
t_stimOff = nan; t_arr = nan; t_reward = nan;
if isempty(idxUsed), return, end
t_stimOff = median(EvRel.stimOff_rel(idxUsed),'omitnan');
t_arr     = median(EvRel.arr_rel(idxUsed),'omitnan');
t_reward  = median(EvRel.reward_rel(idxUsed),'omitnan');
end

function [Lb_mean, usedMask] = ra_eventlocked_mean_spec(V, t_abs, ev, t_rel, bmask_t)
usedMask = false(numel(ev),1);
X = nan(numel(ev), numel(t_rel), size(V,2));
for k = 1:numel(ev)
    tq = ev(k) + t_rel;
    if tq(1) < t_abs(1) || tq(end) > t_abs(end), continue, end
    X(k,:,:) = interp1(t_abs, V, tq, 'linear', NaN);
    usedMask(k) = any(isfinite(X(k,:,:)), 'all');
end
Xuse = X(usedMask,:,:);
M = squeeze(nanmean(Xuse,1)); % [t x f]
base = nanmean(M(bmask_t,:),1);
Lb_mean = M - base;
end

function [piezo_z, delta_z, theta_z] = ra_group_overlay_traces_low_onepiezo(V, t_abs, f_plot, ev, t_rel, bmask_t, t_piezo_abs, env_piezo, opts)
nEv = numel(ev);
PiezoMat = nan(nEv, numel(t_rel));
DelMat   = nan(nEv, numel(t_rel));
TheMat   = nan(nEv, numel(t_rel));

dmask = (f_plot >= opts.ob_band_delta(1)) & (f_plot <= opts.ob_band_delta(2));
tmask = (f_plot >= opts.ob_band_theta(1)) & (f_plot <= opts.ob_band_theta(2));

for k = 1:nEv
    tq = ev(k) + t_rel(:);
    PiezoMat(k,:) = interp1(t_piezo_abs, env_piezo, tq, 'linear', NaN);

    Pk = interp1(t_abs, V, tq, 'linear', NaN); % [t x f]
    basek = nanmean(Pk(bmask_t,:),1);
    Pkb = Pk - basek;

    DelMat(k,:) = nanmean(Pkb(:, dmask), 2);
    TheMat(k,:) = nanmean(Pkb(:, tmask), 2);
end

if opts.overlay_use_median
    p  = nanmedian(PiezoMat,1);
    d  = nanmedian(DelMat,1);
    th = nanmedian(TheMat,1);
else
    p  = nanmean(PiezoMat,1);
    d  = nanmean(DelMat,1);
    th = nanmean(TheMat,1);
end

piezo_z = ra_zscore_safe(p(:));
delta_z = ra_zscore_safe(d(:));
theta_z = ra_zscore_safe(th(:));
end

function [piezo_z, gamma_z, highgamma_z] = ra_group_overlay_traces_mid_onepiezo(V, t_abs, f_plot, ev, t_rel, bmask_t, t_piezo_abs, env_piezo, opts)
nEv = numel(ev);
PiezoMat = nan(nEv, numel(t_rel));
GamMat   = nan(nEv, numel(t_rel));
HGMat    = nan(nEv, numel(t_rel));

gmask = (f_plot >= opts.ob_band_gamma(1)) & (f_plot <= opts.ob_band_gamma(2));
hmask = (f_plot >= opts.ob_band_highgamma(1)) & (f_plot <= opts.ob_band_highgamma(2));

for k = 1:nEv
    tq = ev(k) + t_rel(:);
    PiezoMat(k,:) = interp1(t_piezo_abs, env_piezo, tq, 'linear', NaN);

    Pk = interp1(t_abs, V, tq, 'linear', NaN);
    basek = nanmean(Pk(bmask_t,:),1);
    Pkb = Pk - basek;

    GamMat(k,:) = nanmean(Pkb(:, gmask), 2);
    HGMat(k,:)  = nanmean(Pkb(:, hmask), 2);
end

if opts.overlay_use_median
    p  = nanmedian(PiezoMat,1);
    g  = nanmedian(GamMat,1);
    hg = nanmedian(HGMat,1);
else
    p  = nanmean(PiezoMat,1);
    g  = nanmean(GamMat,1);
    hg = nanmean(HGMat,1);
end

piezo_z     = ra_zscore_safe(p(:));
gamma_z     = ra_zscore_safe(g(:));
highgamma_z = ra_zscore_safe(hg(:));
end

function [w1, w2] = ra_get_window_edges_all(winDefs, EvRel)
n = numel(EvRel.stimOff_rel);
w1 = nan(n, numel(winDefs));
w2 = nan(n, numel(winDefs));
for wi = 1:numel(winDefs)
    wn = lower(string(winDefs{wi}));
    if wn == "stimon_to_stimoff"
        w1(:,wi) = 0;
        w2(:,wi) = EvRel.stimOff_rel;
    elseif wn == "stimoff_to_arr"
        w1(:,wi) = EvRel.stimOff_rel;
        w2(:,wi) = EvRel.arr_rel;
    elseif wn == "arr_to_stop"
        w1(:,wi) = EvRel.arr_rel;
        w2(:,wi) = EvRel.stop_rel;
    end
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

function z = ra_zscore_safe(x)
x = double(x(:));
m = mean(x,'omitnan');
s = std(x,'omitnan');
if ~isfinite(s) || s < 1e-9
    z = x - m;
else
    z = (x - m) ./ s;
end
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

function [XC, peakLag_samp, peakR] = ra_trialwise_xcorr_mat(A, B, wmask, maxLagSamp)
nTr = size(A,1);
XC = nan(nTr, 2*maxLagSamp+1);
peakLag_samp = nan(nTr,1);
peakR = nan(nTr,1);

for tr = 1:nTr
    x = A(tr,wmask);
    y = B(tr,wmask);
    ok = isfinite(x) & isfinite(y);
    if sum(ok) < 10, continue, end
    x = x(ok); y = y(ok);
    x = x - mean(x); y = y - mean(y);
    c = xcorr(x, y, maxLagSamp, 'coeff');
    XC(tr,:) = c;
    [m,im] = max(c);
    peakR(tr) = m;
    peakLag_samp(tr) = im;
end
peakLag_samp = peakLag_samp - (maxLagSamp+1); % samples
end

function phq = ra_interp_phase_circ(t, ph, tq)
u = exp(1i*ph(:));
urq = interp1(t, real(u), tq, 'linear', NaN);
uiq = interp1(t, imag(u), tq, 'linear', NaN);
phq = angle(urq + 1i*uiq);
end

function Y = ra_subtract_baseline_per_trial(X, bmask)
Y = nan(size(X));
for tr = 1:size(X,1)
    xb = X(tr,bmask);
    ok = isfinite(xb);
    if sum(ok) < 2, continue, end
    Y(tr,:) = X(tr,:) - mean(xb(ok));
end
end

function [phCenters, meanPow, nPerBin] = ra_phase_to_power(powMat, phMat, tGrid, win_s, nBins)
wm = (tGrid >= win_s(1)) & (tGrid <= win_s(2));
ph = phMat(:,wm); pw = powMat(:,wm);
ph = ph(:); pw = pw(:);
ok = isfinite(ph) & isfinite(pw);
ph = ph(ok); pw = pw(ok);

edges = linspace(-pi, pi, nBins+1);
phCenters = (edges(1:end-1) + edges(2:end))/2;
meanPow = nan(1,nBins);
nPerBin = zeros(1,nBins);

for b = 1:nBins
    m = (ph >= edges(b)) & (ph < edges(b+1));
    nPerBin(b) = sum(m);
    if nPerBin(b) > 0
        meanPow(b) = mean(pw(m));
    end
end
end

function [respSig, sniffSig] = ra_piezo_two_bands(Resp, fs, opts)
% slow respiration band
FilR = FilterLFP(Resp, opts.resp_band_hz, fs);
tR   = double(Range(FilR,'s')); tR = tR(:);
xR   = double(Data(FilR)); xR = xR(:);

respSig.t_abs = tR;
respSig.xF    = xR;
respSig.ph    = angle(hilbert(zscore(xR)));
respSig.env   = xR.^2;
if isfield(opts,'resp_smooth_s') && opts.resp_smooth_s > 0
    w = max(1, round(opts.resp_smooth_s * fs));
    respSig.env = movmean(respSig.env, w, 'omitnan');
end

% sniff band
FilS = FilterLFP(Resp, opts.sniff_band_hz, fs);
tS   = double(Range(FilS,'s')); tS = tS(:);
xS   = double(Data(FilS)); xS = xS(:);

sniffSig.t_abs = tS;
sniffSig.xF    = xS;
sniffSig.ph    = angle(hilbert(zscore(xS)));
sniffSig.env   = xS.^2;
if isfield(opts,'sniff_smooth_s') && opts.sniff_smooth_s > 0
    w = max(1, round(opts.sniff_smooth_s * fs));
    sniffSig.env = movmean(sniffSig.env, w, 'omitnan');
end
end

function figB = ra_plot_behaviour_style_ob(datapath, sessname, Baphy, G, stim_on_abs_s, ...
    V_L, t_abs_L, f_plot_L, V_M, t_abs_M, f_plot_M, hasMID, ...
    t_rel_L, bmask_t_L, winDefs, w1_all, w2_all, opts, ttl, modeName)

% Create a "behaviour_analysis-like" figure for OB band power classification (Tar vs Ref)
% Rows = bands (opts.behav_bands). Columns mimic your pupil figure:
% (1) band trace timecourse (Tar vs Ref), (2) histogram, (3) meanDiff, (4) CohenD, (5) AUROC text, (6) ROC

behW = lower(string(opts.behav_window));
% map user shortcut names
if behW == "stimoff_to_arrival", behW = "stimoff_to_arr"; end
if behW == "stimon_to_stimoff",  behW = "stimon_to_stimoff"; end

winIdx = find(lower(string(winDefs)) == behW, 1, 'first');
if isempty(winIdx)
    warning('behav_window not found: %s', opts.behav_window);
    figB = [];
    return
end

ciTar = find(strcmpi(G.names,'tar'),1);
ciRef = find(strcmpi(G.names,'ref'),1);
if isempty(ciTar) || isempty(ciRef)
    warning('Need G.names contain tar/ref for behaviour-style figure.');
    figB = [];
    return
end

bands = opts.behav_bands;
nB = numel(bands);

figB = figure('Color','w','Units','pixels','Position',[20 60 1900 850], ...
    'Name', sprintf('%s | OB behaviour-style | %s | %s', sessname, opts.behav_window, modeName));

sw = opts.behav_sw_s;
dt = median(diff(t_rel_L));
swSamp = max(1, round(sw / dt));

for bi = 1:nB
    bname = lower(string(bands{bi}));

    % bandRange + choose spec
    useMid = false;
    if bname == "delta"
        br = opts.ob_band_delta;
    elseif bname == "theta"
        br = opts.ob_band_theta;
    elseif bname == "gamma"
        br = opts.ob_band_gamma; useMid = true;
    elseif bname == "highgamma"
        br = opts.ob_band_highgamma; useMid = true;
    else
        continue
    end

    % choose V/t/f
    if useMid && hasMID
        V = V_M; t_abs = t_abs_M; f_plot = f_plot_M;
        bmask_grid = (t_rel_L >= opts.spec_baseline_s(1)) & (t_rel_L <= opts.spec_baseline_s(2));
    else
        V = V_L; t_abs = t_abs_L; f_plot = f_plot_L;
        bmask_grid = bmask_t_L;
    end

    % build band trace per trial (full t_rel_L grid)
    idxTar = find(G.mask(:,ciTar) & isfinite(stim_on_abs_s));
    idxRef = find(G.mask(:,ciRef) & isfinite(stim_on_abs_s));

    BTar = ra_bandtrace_from_spec_to_grid(V, t_abs, f_plot, stim_on_abs_s(idxTar), t_rel_L(:), bmask_grid, br);
    BRef = ra_bandtrace_from_spec_to_grid(V, t_abs, f_plot, stim_on_abs_s(idxRef), t_rel_L(:), bmask_grid, br);

    mTar = movmean(mean(BTar,1,'omitnan'), swSamp, 'omitnan');
    mRef = movmean(mean(BRef,1,'omitnan'), swSamp, 'omitnan');
    sTar = movmean(std(BTar,[],1,'omitnan'),  swSamp, 'omitnan') ./ sqrt(max(1,size(BTar,1)));
    sRef = movmean(std(BRef,[],1,'omitnan'),  swSamp, 'omitnan') ./ sqrt(max(1,size(BRef,1)));

    % scalar value per trial in chosen window (mean over time in that window)
    t1T = w1_all(idxTar, winIdx); t2T = w2_all(idxTar, winIdx);
    t1R = w1_all(idxRef, winIdx); t2R = w2_all(idxRef, winIdx);

    xTar = nan(numel(idxTar),1);
    for k = 1:numel(idxTar)
        m = (t_rel_L >= t1T(k)) & (t_rel_L <= t2T(k));
        if sum(m) < 5, continue, end
        xTar(k) = mean(BTar(k,m),'omitnan');
    end

    xRef = nan(numel(idxRef),1);
    for k = 1:numel(idxRef)
        m = (t_rel_L >= t1R(k)) & (t_rel_L <= t2R(k));
        if sum(m) < 5, continue, end
        xRef(k) = mean(BRef(k,m),'omitnan');
    end

    xTar = xTar(isfinite(xTar));
    xRef = xRef(isfinite(xRef));

    M = ra_effect_metrics_simple(xTar, xRef);

    % Row panels
    r0 = (bi-1)*6;

    % (1) timecourse
    subplot(nB,6,r0+1); hold on
    plot(t_rel_L, mRef, 'LineWidth',1.5);
    plot(t_rel_L, mTar, 'LineWidth',1.5);
    % lightweight shaded CI
    ra_shade(t_rel_L, mRef, sRef);
    ra_shade(t_rel_L, mTar, sTar);
    xline(0,'k-','LineWidth',1);
    yline(0,'k:');
    xlabel('Time from stimOn (s)');
    ylabel(char(bname));
    title(sprintf('OB %s (Tar=%d Ref=%d)', bname, numel(xTar), numel(xRef)), 'Interpreter','none');
    box off
    legend({'Ref','Tar'},'Box','off');
    hold off

    % (2) histogram
    subplot(nB,6,r0+2); hold on
    histogram(xRef, 30, 'Normalization','count');
    histogram(xTar, 30, 'Normalization','count');
    title(sprintf('AshmanD=%.2f', M.AshmanD), 'Interpreter','none');
    box off
    hold off

    % (3) meanDiff
    subplot(nB,6,r0+3);
    ra_spread_two(xRef, xTar);
    title(sprintf('meanDiff (Tar-Ref)=%.4g', M.meanDiff),'Interpreter','none');
    box off

    % (4) Cohen d
    subplot(nB,6,r0+4);
    ra_spread_two(xRef, xTar);
    title(sprintf('Cohen''s d = %.2f', M.CohensD),'Interpreter','none');
    box off

    % (5) AUROC numeric
    subplot(nB,6,r0+5);
    ra_spread_two(xRef, xTar);
    title(sprintf('AUROC = %.2f', M.AUROC),'Interpreter','none');
    box off

    % (6) ROC
    subplot(nB,6,r0+6); hold on
    plot(M.rocFPR, M.rocTPR, 'k-', 'LineWidth',1.5);
    plot([0 1],[0 1],'k:');
    xlabel('FPR'); ylabel('TPR');
    title(sprintf('AUC=%.2f', M.AUROC),'Interpreter','none');
    axis([0 1 0 1]); box off
    hold off
end

sgtitle(sprintf('%s | OB behaviour-style | %s | %s', ttl, opts.behav_window, upper(modeName)), 'Interpreter','none');
end

function M = ra_effect_metrics_simple(xTar, xRef)
% Returns AUROC, Cohen's d, meanDiff, AshmanD and an ROC curve (FPR/TPR).
xTar = double(xTar(:)); xRef = double(xRef(:));
xTar = xTar(isfinite(xTar)); xRef = xRef(isfinite(xRef));

M = struct('AUROC',nan,'CohensD',nan,'meanDiff',nan,'AshmanD',nan,'rocFPR',[],'rocTPR',[]);

if numel(xTar) < 3 || numel(xRef) < 3
    return
end

M.meanDiff = mean(xTar) - mean(xRef);

% Cohen d
m1 = mean(xTar); m0 = mean(xRef);
s1 = var(xTar,1); s0 = var(xRef,1);
sp = sqrt(0.5*(s1+s0));
if sp > 0
    M.CohensD = (m1 - m0) / sp;
else
    M.CohensD = 0;
end

% AshmanD
if s1 > 0 && s0 > 0
    M.AshmanD = sqrt(2) * abs(m1 - m0) / sqrt(s1 + s0);
else
    M.AshmanD = 0;
end

% AUROC via Mann–Whitney (rank) + ROC points
y = [zeros(numel(xRef),1); ones(numel(xTar),1)];
s = [xRef; xTar];

% ranks
[~,~,r] = unique(s);
r = double(r);

n0 = numel(xRef);
n1 = numel(xTar);
R1 = sum(r(n0+1:end));
U1 = R1 - n1*(n1+1)/2;
M.AUROC = U1 / (n0*n1);

% ROC curve using thresholds on score
thr = unique(s);
thr = sort(thr,'ascend');
FPR = nan(numel(thr)+2,1);
TPR = nan(numel(thr)+2,1);
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

% ensure monotonic by sorting on FPR
[FF,ord] = sort(FPR);
TT = TPR(ord);
M.rocFPR = FF;
M.rocTPR = TT;
end

function ra_spread_two(xRef, xTar)
% Simple spread + medians (no external deps).
xRef = xRef(:); xTar = xTar(:);
j1 = (rand(size(xRef))-0.5)*0.2;
j2 = (rand(size(xTar))-0.5)*0.2;
plot(1+j1, xRef, '.', 'MarkerSize',8); hold on
plot(2+j2, xTar, '.', 'MarkerSize',8);
m1 = median(xRef,'omitnan'); m2 = median(xTar,'omitnan');
plot([0.8 1.2],[m1 m1],'k-','LineWidth',2);
plot([1.8 2.2],[m2 m2],'k-','LineWidth',2);
set(gca,'XTick',[1 2],'XTickLabel',{'Ref','Tar'});
xlim([0.5 2.5]);
end

function ra_shade(t, m, s)
t = t(:); m = m(:); s = s(:);
ok = isfinite(t) & isfinite(m) & isfinite(s);
t = t(ok); m = m(ok); s = s(ok);
if numel(t) < 3, return, end
x = [t; flipud(t)];
y = [m-s; flipud(m+s)];
h = fill(x,y,'k');
set(h,'FaceAlpha',0.10,'EdgeColor','none');
end