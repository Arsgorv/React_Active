function Met = eventlocked_spec_grid_plot(Spectro, Baphy, event_abs_s, G, opts, ttl, alignName)
% eventlocked_spec_grid_plot
% Event-locked mean spectrogram + band-power traces + per-trial window metrics.
%
% Key properties:
% - Does NOT skip groups if n is small or zero. Metrics become NaN with n=0.
% - Uses EXACT 3 main windows: [stimOn stimOff], [stimOff arrival], [arrival stop]
% - Optional extra windows defined relative to alignment time 0 (typically stimOn).
% - Bands are supplied by opts.metricBands / opts.metricBandNames (LOW vs MID handled upstream).
%
% Aesthetics updates:
% - Row 1: event lines (stimOn/stimOff/arrival/reward), viridis colormap, shared caxis
% - Row 2: legend for band traces, same event lines
% - Rows 3+: use MakeSpreadAndBoxPlot3_SB distributions (instead of bars), shared ylim per band-row

if nargin < 7 || isempty(alignName), alignName = 'stimOn'; end
if nargin < 6, ttl = ''; end
if nargin < 5, opts = struct(); end

% Defaults
if ~isfield(opts,'spec_twin_s'),        opts.spec_twin_s = [-0.5 7.0]; end
if ~isfield(opts,'spec_baseline_s'),    opts.spec_baseline_s = [-0.3 0]; end
if ~isfield(opts,'spec_freq_xlim'),     opts.spec_freq_xlim = [0.5 10]; end
if ~isfield(opts,'spec_smooth_t'),      opts.spec_smooth_t = 0; end
if ~isfield(opts,'spec_smooth_f'),      opts.spec_smooth_f = 0; end
if ~isfield(opts,'min_events'),         opts.min_events = 5; end

% Bands MUST be provided (so LOW/MID is clean and no empty [] metrics)
if ~isfield(opts,'metricBands') || isempty(opts.metricBands)
    error('opts.metricBands must be provided (e.g., [0.5 4; 4 6] for LOW, [40 60; 60 80] for MID).');
end
if ~isfield(opts,'metricBandNames') || isempty(opts.metricBandNames)
    nb = size(opts.metricBands,1);
    tmp = cell(nb,1);
    for i = 1:nb, tmp{i} = sprintf('band%d',i); end
    opts.metricBandNames = tmp;
end

% Extra windows relative to alignment time 0
if ~isfield(opts,'relWinEdges_s'), opts.relWinEdges_s = []; end
if ~isfield(opts,'relWinNames'),   opts.relWinNames = {}; end

% ---- Extract (t,f,Praw) from Spectro robustly, then standardize ----
[t, f_all, Praw] = ra_extract_spectro_triplet(Spectro);
[t_abs_s, S_tf]  = standardize_spectro_time_and_orient(t, f_all, Praw);

if size(S_tf,1) ~= numel(t_abs_s) || size(S_tf,2) ~= numel(f_all)
    error('After standardize: size(P)=[%d %d], numel(t)=%d, numel(f)=%d', ...
        size(S_tf,1), size(S_tf,2), numel(t_abs_s), numel(f_all));
end

% Restrict freq
fmask = f_all >= opts.spec_freq_xlim(1) & f_all <= opts.spec_freq_xlim(2);
f_plot = f_all(fmask);
S_tf = S_tf(:, fmask); % [time x freq]

% Build time axis relative to event
dt_spec = median(diff(t_abs_s));
nRel = round((opts.spec_twin_s(2)-opts.spec_twin_s(1)) / max(dt_spec, eps)) + 1;
t_rel = linspace(opts.spec_twin_s(1), opts.spec_twin_s(2), nRel);

% Baseline mask in relative time
bmask_t = t_rel >= opts.spec_baseline_s(1) & t_rel <= opts.spec_baseline_s(2);

% ---- Relative event times for vertical lines (use REL fields; robust) ----
Ttr = Baphy.trial;

% If alignName is stimOn, alignment is rel_stim_start; we still compute relative-to-alignment via rel fields
stimOn_rel  = zeros(Baphy.n_trials,1);
stimOff_rel = nan(Baphy.n_trials,1);
arr_rel     = nan(Baphy.n_trials,1);
reward_rel  = nan(Baphy.n_trials,1);
stop_rel    = nan(Baphy.n_trials,1);

if isfield(Ttr,'rel_stim_start') && isfield(Ttr,'rel_stim_stop')
    stimOff_rel = Ttr.rel_stim_stop(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_stim_stop')
    stimOff_rel = Ttr.rel_stim_stop(:);
end
if isfield(Ttr,'spout_arrival_rel_s') && isfield(Ttr,'rel_stim_start')
    arr_rel = Ttr.spout_arrival_rel_s(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'spout_arrival_rel_s')
    arr_rel = Ttr.spout_arrival_rel_s(:);
end
if isfield(Ttr,'rel_reward') && isfield(Ttr,'rel_stim_start')
    reward_rel = Ttr.rel_reward(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_reward')
    reward_rel = Ttr.rel_reward(:);
end
if isfield(Ttr,'rel_trialend') && isfield(Ttr,'rel_stim_start')
    stop_rel = Ttr.rel_trialend(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_trialend')
    stop_rel = Ttr.rel_trialend(:);
end

% Define EXACT 3 main windows (per trial)
winDefs = { ...
    'stimOn_to_stimOff'; ...
    'stimOff_to_arr'; ...
    'arr_to_stop' ...
};

% Add optional relative windows (same for all trials, relative to t=0)
extraWinDefs = {};
if ~isempty(opts.relWinEdges_s)
    if size(opts.relWinEdges_s,2) ~= 2
        error('opts.relWinEdges_s must be [n x 2].');
    end
    nWextra = size(opts.relWinEdges_s,1);
    if isempty(opts.relWinNames)
        opts.relWinNames = cell(nWextra,1);
        for k = 1:nWextra
            opts.relWinNames{k} = sprintf('rel_%.3g_to_%.3g', opts.relWinEdges_s(k,1), opts.relWinEdges_s(k,2));
        end
    end
    if numel(opts.relWinNames) ~= nWextra
        error('opts.relWinNames must match number of rows in opts.relWinEdges_s.');
    end
    extraWinDefs = opts.relWinNames(:);
end
allWinNames = [winDefs; extraWinDefs];

% Prepare Met struct
Met = struct();
Met.groupNames      = G.names;
Met.metricBands     = opts.metricBands;
Met.metricBandNames = opts.metricBandNames;
Met.winDefs         = allWinNames;
Met.opts_snapshot   = opts;
Met.alignName       = alignName;
Met.tag             = '';
if isfield(opts,'tag'), Met.tag = opts.tag; end

% Figure layout
nGroups = numel(G.names);
nCols = nGroups;
nB = size(opts.metricBands,1);
nRows = 2 + nB; % spectrogram + band traces + one row per band distributions

fig = figure('Color','w','Name',ttl);
set(fig, 'Units','pixels', 'Position',[50 50 1800 950]);
set(fig, 'PaperPositionMode','auto');

% Global ranges
cax_all = [inf -inf];
ylim_all = [inf -inf];

% Precompute eventlocked matrices per group
groupData = cell(nGroups,1);
groupUsed = cell(nGroups,1);
groupIdx  = cell(nGroups,1);

for c = 1:nGroups
    idx_all = find(G.mask(:,c));
    ev_all  = event_abs_s(idx_all);
    valid = isfinite(ev_all);
    idx = idx_all(valid);
    ev  = ev_all(valid);

    groupIdx{c} = idx;

    if isempty(ev)
        groupData{c} = [];
        groupUsed{c} = false(0,1);
        continue
    end

    X = nan(numel(ev), numel(t_rel), numel(f_plot));
    used = false(numel(ev),1);

    for k = 1:numel(ev)
        tq = ev(k) + t_rel;
        if tq(1) < t_abs_s(1) || tq(end) > t_abs_s(end)
            continue
        end
        X(k,:,:) = interp1(t_abs_s, S_tf, tq, 'linear', NaN);
        used(k) = all(isfinite(X(k,1,:))) && all(isfinite(X(k,end,:)));
    end

    groupData{c} = X;
    groupUsed{c} = used;
end

% Storage to unify y-lims for distribution rows
axWin   = cell(nB, nGroups);
winVals = cell(nB, nGroups);

% Main plots + metrics fill
Met.group = struct([]);

for c = 1:nGroups
    gname = string(G.names{c});
    idx_used = groupIdx{c};
    X = groupData{c};
    used = groupUsed{c};

    if isempty(X)
        Xuse = nan(0, numel(t_rel), numel(f_plot));
        idx_used = zeros(0,1);
    else
        Xuse = X(used,:,:);
        idx_used = idx_used(used);
    end
    nEv = size(Xuse,1);

    % Group-specific median event lines (relative to stimOn=0)
    t_stimOn  = 0;
    t_stimOff = median(stimOff_rel(idx_used), 'omitnan');
    t_arr     = median(arr_rel(idx_used),     'omitnan');
    t_reward  = median(reward_rel(idx_used),  'omitnan');

    % Mean spectrogram (baseline-corrected log-power)
    M = squeeze(nanmean(Xuse,1)); % [time x freq]
    M(M<=0) = NaN;
    L = log10(M);
    base = nanmean(L(bmask_t,:),1);
    Lb = L - base;

    if opts.spec_smooth_t > 0
        Lb = movmean(Lb, max(1, round(opts.spec_smooth_t/median(diff(t_rel)))), 1, 'omitnan');
    end
    if opts.spec_smooth_f > 0
        Lb = movmean(Lb, max(1, round(opts.spec_smooth_f/median(diff(f_plot)))), 2, 'omitnan');
    end

    if nEv > 0
        cax_all(1) = min(cax_all(1), nanmin(Lb(:)));
        cax_all(2) = max(cax_all(2), nanmax(Lb(:)));
    end

    % ---- Row 1: Spectrogram
    ax1 = subplot(nRows, nCols, c);
    imagesc(t_rel, f_plot, Lb'); axis xy
    title(sprintf('%s (n=%d)', gname, nEv), 'Interpreter','none');
    ylabel('Hz'); xlabel(sprintf('t from %s (s)', alignName), 'Interpreter','none');
    box off
    % colormap viridis
    try
        colormap(ax1,'viridis');
    catch
        try
            colormap(ax1, viridis);
        catch
            % ignore
        end
    end
    % event lines
    xline(ax1, t_stimOn,'k-','LineWidth',1);
    if isfinite(t_stimOff), xline(ax1, t_stimOff,'k--','LineWidth',1); end
    if isfinite(t_arr),     xline(ax1, t_arr,    'k--','LineWidth',1); end
    if isfinite(t_reward),  xline(ax1, t_reward, 'k--','LineWidth',1); end

    % ---- Row 2: Band traces
    Met.group(c).name = char(gname);
    Met.group(c).nEvents = nEv;

    nW = numel(allWinNames);
    valMat = nan(nEv, nB, nW);
    slpMat = nan(nEv, nB, nW);

    hBand = gobjects(nB,1);
    bandLabel = cell(nB,1);

    ax2 = subplot(nRows, nCols, nCols + c);
    hold(ax2,'on');

    for b = 1:nB
        br = opts.metricBands(b,:);
        bname = string(opts.metricBandNames{b});
        bmask_f = f_plot >= br(1) & f_plot <= br(2);

        % band trace per trial (baseline-corrected in log10 domain)
        if nEv == 0
            bt_mean = nan(size(t_rel));
            bt_sem  = nan(size(t_rel));
        else
            Xb = nan(nEv, numel(t_rel));
            for k = 1:nEv
                Sk = squeeze(Xuse(k,:,:)); % [time x freq]
                Sk(Sk<=0) = NaN;
                Lk = log10(Sk);
                basek = nanmean(Lk(bmask_t,:),1);
                Lkb = Lk - basek;
                Xb(k,:) = nanmean(Lkb(:,bmask_f),2);
            end
            bt_mean = nanmean(Xb,1);
            nEff = sum(isfinite(Xb), 1);
            bt_sem = nanstd(Xb, [], 1) ./ max(1, sqrt(nEff));
        end

        hBand(b) = plot(ax2, t_rel, bt_mean, 'LineWidth', 1);
        bandLabel{b} = char(bname);

        if nEv > 0
            ylim_all(1) = min(ylim_all(1), nanmin(bt_mean - bt_sem));
            ylim_all(2) = max(ylim_all(2), nanmax(bt_mean + bt_sem));
        end

        % Fill Met band structs
        Met.group(c).band(b).name  = char(bname);
        Met.group(c).band(b).range = br;
        Met.group(c).band(b).win   = struct([]);

        % ---- Window metrics: EXACT 3 + optional relative windows
        for w = 1:nW
            wname = string(allWinNames{w});
            [w1, w2] = ra_get_window_edges(wname, stimOn_rel, stimOff_rel, arr_rel, stop_rel, opts);

            if nEv == 0
                v = nan(0,1);
                s = nan(0,1);
            else
                v = nan(nEv,1);
                s = nan(nEv,1);
                for k = 1:nEv
                    t1 = w1(idx_used(k));
                    t2 = w2(idx_used(k));
                    if ~isfinite(t1) || ~isfinite(t2) || t2 <= t1
                        continue
                    end
                    tm = (t_rel >= t1) & (t_rel <= t2);
                    if ~any(tm), continue, end

                    Sk = squeeze(Xuse(k,:,:)); % [time x freq]
                    Lk = log10(Sk);
                    basek = nanmean(Lk(bmask_t,:),1);
                    Lkb = Lk - basek;

                    y = nanmean(Lkb(:,bmask_f),2);
                    y = y(:);

                    v(k) = nanmean(y(tm));

                    tt = t_rel(tm); tt = tt(:);
                    yy = y(tm);     yy = yy(:);
                    if numel(tt) ~= numel(yy) || isempty(tt), continue, end
                    good = isfinite(tt) & isfinite(yy);
                    if sum(good) >= 3
                        p = polyfit(tt(good), yy(good), 1);
                        s(k) = p(1);
                    end
                end
            end

            valMat(:,b,w) = v;
            slpMat(:,b,w) = s;

            Met.group(c).band(b).win(w).name = char(wname);
            Met.group(c).band(b).win(w).mean = nanmean(v);
            Met.group(c).band(b).win(w).sem  = nanstd(v) ./ max(1, sqrt(sum(isfinite(v))));
            Met.group(c).band(b).win(w).n    = sum(isfinite(v));
            Met.group(c).band(b).win(w).slope_mean = nanmean(s);
            Met.group(c).band(b).win(w).slope_sem  = nanstd(s) ./ max(1, sqrt(sum(isfinite(s))));
            Met.group(c).band(b).win(w).slope_n    = sum(isfinite(s));
            Met.group(c).band(b).trace.t_rel = t_rel(:);
            Met.group(c).band(b).trace.mean = bt_mean(:);
            Met.group(c).band(b).trace.sem  = bt_sem(:);
        end

        % ---- Rows 3+: distributions using MakeSpreadAndBoxPlot3_SB (only the 3 main windows)
        axb = subplot(nRows, nCols, (2+b-1)*nCols + c);
        cla(axb); axes(axb); %#ok<LAXES>

        A = cell(1,3);
        for ww = 1:3
            A{ww} = valMat(:,b,ww);
        end
        Leg = {'stimOn\rightarrowstimOff','stimOff\rightarrowarr','arr\rightarrowstop'};

        try
            MakeSpreadAndBoxPlot3_SB(A, {}, [1 2 3], Leg, 'newfig',0, 'paired',0, 'showpoints',1, 'showsigstar','none');
        catch
            text(0.1,0.5,'MakeSpreadAndBoxPlot3_SB not found on path','Units','normalized');
            set(gca,'XTick',[1 2 3],'XTickLabel',Leg);
        end
        if c==1
            ylabel(sprintf('%s mean', char(bname)), 'Interpreter','none');
        end
        yline(0,'k:'); box off
        set(gca,'XTickLabelRotation',20);

        axWin{b,c}   = gca;
        winVals{b,c} = [A{1}(:); A{2}(:); A{3}(:)];
    end

    % Band-trace axes cosmetics + legend + event lines (once per group)
    yline(ax2, 0,'k:');
    xline(ax2, t_stimOn,'k-','LineWidth',1);
    if isfinite(t_stimOff), xline(ax2, t_stimOff,'k--','LineWidth',1); end
    if isfinite(t_arr),     xline(ax2, t_arr,    'k--','LineWidth',1); end
    if isfinite(t_reward),  xline(ax2, t_reward, 'k--','LineWidth',1); end
    hold(ax2,'off');
    box(ax2,'off');
    if c==1, ylabel(ax2,'log10 power (BL-corr)'); end
    if c==1, title(ax2, sprintf('Band traces (%s)', alignName), 'Interpreter','none'); end

    goodH = isgraphics(hBand);
    if any(goodH)
        legend(ax2, hBand(goodH), bandLabel(goodH), 'Location','northwest', 'Box','off');
    end

    % Save per-trial matrices for export
    Met.group(c).trialIdx_used = idx_used(:);
    Met.group(c).trialIdxUsed  = idx_used(:); % alias
    Met.group(c).valMat = valMat;
    Met.group(c).slpMat = slpMat;
end

% Global cosmetics: shared caxis + ylims for band traces
for c = 1:nGroups
    ax = subplot(nRows, nCols, c);
    if isfinite(cax_all(1)) && isfinite(cax_all(2)) && cax_all(2) > cax_all(1)
        caxis(ax, cax_all);
    end
end
for c = 1:nGroups
    ax = subplot(nRows, nCols, nCols + c);
    if isfinite(ylim_all(1)) && isfinite(ylim_all(2)) && ylim_all(2) > ylim_all(1)
        ylim(ax, ylim_all);
    end
end

% Shared y-lims for distribution rows (per band)
for b = 1:nB
    allv = [];
    for c = 1:nGroups
        v = winVals{b,c};
        if ~isempty(v)
            allv = [allv; v(isfinite(v))]; %#ok<AGROW>
        end
    end
    if isempty(allv), continue, end
    ymin = min(allv); ymax = max(allv);
    pad = 0.05 * max(1e-6, ymax-ymin);
    yl = [min(ymin,0)-pad, max(ymax,0)+pad];
    for c = 1:nGroups
        axb = axWin{b,c};
        if ~isempty(axb) && isgraphics(axb)
            ylim(axb, yl);
        end
    end
end

% Figure-wide font and axes styling
set(findall(fig,'-property','FontSize'),'FontSize',10);
set(findall(fig,'type','axes'),'TickDir','out','Box','off');

sgtitle(ttl,'Interpreter','none');

end

% ---------------- helper: window edges per trial ----------------
function [w1, w2] = ra_get_window_edges(wname, stimOn_rel, stimOff_rel, arr_rel, stop_rel, opts)
w1 = nan(size(stimOn_rel));
w2 = nan(size(stimOn_rel));

wn = lower(string(wname));

if wn == "stimon_to_stimoff"
    w1 = stimOn_rel;
    w2 = stimOff_rel;
    return
end
if wn == "stimoff_to_arr"
    w1 = stimOff_rel;
    w2 = arr_rel;
    return
end
if wn == "arr_to_stop"
    w1 = arr_rel;
    w2 = stop_rel;
    return
end

% optional relative windows
if isfield(opts,'relWinNames') && isfield(opts,'relWinEdges_s') && ~isempty(opts.relWinNames)
    nm = lower(string(opts.relWinNames(:)));
    k = find(nm == wn, 1, 'first');
    if ~isempty(k)
        e = opts.relWinEdges_s(k,:);
        w1(:) = e(1);
        w2(:) = e(2);
        return
    end
end

end

function [t, f, Praw] = ra_extract_spectro_triplet(Spectro)
% Supports your convention: Spectro = {P, t, f} (but also works if permuted)

t = []; f = []; Praw = [];

if isstruct(Spectro)
    if isfield(Spectro,'t'), t = Spectro.t; end
    if isempty(t) && isfield(Spectro,'time'), t = Spectro.time; end

    if isfield(Spectro,'f'), f = Spectro.f; end
    if isempty(f) && isfield(Spectro,'freq'), f = Spectro.freq; end

    if isfield(Spectro,'P'), Praw = Spectro.P; end
    if isempty(Praw) && isfield(Spectro,'S'), Praw = Spectro.S; end

    if isempty(t) || isempty(f) || isempty(Praw)
        error('Spectro struct: could not find fields t/f/P (or S).');
    end

    t = double(t(:));
    f = double(f(:));
    Praw = double(Praw);
    return
end

if iscell(Spectro)
    if numel(Spectro) < 3
        error('Spectro cell: expected at least 3 elements.');
    end

    A = Spectro{1}; B = Spectro{2}; C = Spectro{3};
    if ~(isnumeric(A) && isnumeric(B) && isnumeric(C))
        error('Spectro cell: first 3 elements must be numeric.');
    end

    A = double(A); B = double(B); C = double(C);

    isMat = [~isvector(A), ~isvector(B), ~isvector(C)];
    if sum(isMat) ~= 1
        error('Spectro cell: need exactly one matrix + two vectors.');
    end

    if isMat(1)
        P = A; v1 = B(:); v2 = C(:);
    elseif isMat(2)
        P = B; v1 = A(:); v2 = C(:);
    else
        P = C; v1 = A(:); v2 = B(:);
    end

    n1 = numel(v1); n2 = numel(v2);
    [nR,nC] = size(P);

    ok = (n1==nR && n2==nC) || (n1==nC && n2==nR);
    if ~ok
        error('Spectro cell: dims mismatch: size(P)=[%d %d], numel(v1)=%d, numel(v2)=%d', nR,nC,n1,n2);
    end

    if n1 == nR
        t = v1; f = v2; Praw = P;
    else
        t = v2; f = v1; Praw = P;
    end

    t = double(t(:));
    f = double(f(:));
    Praw = double(Praw);
    return
end

error('Unknown Spectro type.');
end
