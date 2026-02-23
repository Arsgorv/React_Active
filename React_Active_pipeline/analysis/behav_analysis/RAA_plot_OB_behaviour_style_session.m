% Per-session "behaviour-style" OB figure (Tar vs Ref) for band power in a chosen window.
% Assumes you already have:
%   - datapath, Baphy, stim_on_abs_s
%   - G with G.names containing 'tar' and 'ref' and G.mask [nTrials x nGroups]
%   - SpectroLow already loaded (or pass the one you use in ob_event_analysis)
%
% Usage:
%   fig = RAA_plot_OB_behaviour_style_session(datapath, SpectroLow, Baphy, stim_on_abs_s, G, opts);

function fig = RAA_plot_OB_behaviour_style_session(datapath, SpectroLow, Baphy, stim_on_abs_s, G, opts)

if nargin < 6 || isempty(opts), opts = struct(); end

% defaults
if ~isfield(opts,'ttl'),               opts.ttl = ''; end
if ~isfield(opts,'spec_twin_whole_s'), opts.spec_twin_whole_s = [-0.5 7.0]; end
if ~isfield(opts,'spec_baseline_s'),   opts.spec_baseline_s   = [-0.5 -0.1]; end
if ~isfield(opts,'spec_freq_xlim_low'),opts.spec_freq_xlim_low = [0.5 12]; end
if ~isfield(opts,'spec_freq_xlim_mid'),opts.spec_freq_xlim_mid = [20 100]; end
if ~isfield(opts,'ob_band_delta'),     opts.ob_band_delta = [0.5 4]; end
if ~isfield(opts,'ob_band_theta'),     opts.ob_band_theta = [4 12]; end
if ~isfield(opts,'ob_band_gamma'),     opts.ob_band_gamma = [40 60]; end
if ~isfield(opts,'ob_band_highgamma'), opts.ob_band_highgamma = [60 80]; end
if ~isfield(opts,'windowName'),        opts.windowName = 'stimoff_to_arr'; end % 'stimon_to_stimoff'|'stimoff_to_arr'|'arr_to_stop'
if ~isfield(opts,'bands'),             opts.bands = {'delta','theta','gamma','highgamma'}; end
if ~isfield(opts,'sw_s'),              opts.sw_s = 0.10; end

[~, sessname] = fileparts(datapath);

ciTar = find(strcmpi(G.names,'tar'),1);
ciRef = find(strcmpi(G.names,'ref'),1);
if isempty(ciTar) || isempty(ciRef)
    error('G.names must contain ''tar'' and ''ref''');
end

EvRel = ra_get_trial_rel_times(Baphy);
[winDefs, w1_all, w2_all] = ra_get_window_edges_from_EvRel(EvRel);
wi = find(strcmpi(winDefs, opts.windowName), 1);
if isempty(wi), error('Unknown windowName: %s', opts.windowName); end

% --- parse LOW spectrogram
[tL, fL, PLraw] = ra_extract_spectro_triplet(SpectroLow);
[t_abs_L, P_tf_L] = ra_standardize_spectro_time_and_orient(tL, fL, PLraw);
P_tf_L = log10(max(double(P_tf_L), eps));

dtL = median(diff(t_abs_L));
t_rel = (opts.spec_twin_whole_s(1):dtL:opts.spec_twin_whole_s(2))';
bmask = (t_rel >= opts.spec_baseline_s(1)) & (t_rel <= opts.spec_baseline_s(2));

fmaskL = (fL >= opts.spec_freq_xlim_low(1)) & (fL <= opts.spec_freq_xlim_low(2));
f_plot_L = fL(fmaskL);
V_L = P_tf_L(:, fmaskL);

% --- load MID if exists
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

% --- pick Tar/Ref trials
idxTar = find(G.mask(:,ciTar) & isfinite(stim_on_abs_s));
idxRef = find(G.mask(:,ciRef) & isfinite(stim_on_abs_s));

% --- figure
nB = numel(opts.bands);
fig = figure('Color','w','Units','pixels','Position',[20 60 1900 250*nB], ...
    'Name', sprintf('%s | OB behaviour-style | %s', sessname, opts.windowName));

dt = median(diff(t_rel));
swSamp = max(1, round(opts.sw_s / dt));

for bi = 1:nB
    bname = lower(string(opts.bands{bi}));

    [useMid, br] = ra_band_range_from_name(bname, opts);
    if useMid && ~hasMID
        continue
    end

    if useMid
        V = V_M; t_abs = t_abs_M; f_plot = f_plot_M;
        bmask_grid = (t_rel >= opts.spec_baseline_s(1)) & (t_rel <= opts.spec_baseline_s(2));
    else
        V = V_L; t_abs = t_abs_L; f_plot = f_plot_L;
        bmask_grid = bmask;
    end

    % full timecourse band-power matrices (baseline-sub per trial)
    BTar = ra_bandtrace_from_spec_to_grid(V, t_abs, f_plot, stim_on_abs_s(idxTar), t_rel, bmask_grid, br);
    BRef = ra_bandtrace_from_spec_to_grid(V, t_abs, f_plot, stim_on_abs_s(idxRef), t_rel, bmask_grid, br);

    % mean traces + SEM (smoothed)
    mTar = movmean(mean(BTar,1,'omitnan'), swSamp, 'omitnan');
    mRef = movmean(mean(BRef,1,'omitnan'), swSamp, 'omitnan');
    sTar = movmean(std(BTar,[],1,'omitnan'), swSamp, 'omitnan') ./ sqrt(max(1,size(BTar,1)));
    sRef = movmean(std(BRef,[],1,'omitnan'), swSamp, 'omitnan') ./ sqrt(max(1,size(BRef,1)));

    % scalar per trial in the chosen window (trial-specific edges)
    xTar = ra_trial_scalar_in_window(BTar, t_rel, w1_all(idxTar,wi), w2_all(idxTar,wi));
    xRef = ra_trial_scalar_in_window(BRef, t_rel, w1_all(idxRef,wi), w2_all(idxRef,wi));

    M = ra_effect_metrics_simple(xTar, xRef);

    r0 = (bi-1)*6;

    % (1) traces
    subplot(nB,6,r0+1); hold on
    plot(t_rel, mRef, 'LineWidth',1.5);
    plot(t_rel, mTar, 'LineWidth',1.5);
    ra_shade(t_rel, mRef, sRef);
    ra_shade(t_rel, mTar, sTar);
    xline(0,'k-','LineWidth',1);
    yline(0,'k:');
    xlabel('Time from stimOn (s)');
    ylabel(sprintf('%s pow (BL-sub)', bname));
    title(sprintf('%s (Tar=%d Ref=%d)', bname, numel(xTar), numel(xRef)), 'Interpreter','none');
    legend({'Ref','Tar'},'Box','off');
    box off
    hold off

    % (2) histogram
    subplot(nB,6,r0+2); hold on
    histogram(xRef, 30);
    histogram(xTar, 30);
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

sgtitle(sprintf('%s | OB behaviour-style | %s', opts.ttl, opts.windowName), 'Interpreter','none');

end

% ================= helpers =================

function [useMid, br] = ra_band_range_from_name(bname, opts)
useMid = false;
if bname == "delta"
    br = opts.ob_band_delta;
elseif bname == "theta"
    br = opts.ob_band_theta;
elseif bname == "gamma"
    br = opts.ob_band_gamma; useMid = true;
elseif bname == "highgamma" || bname == "high_gamma" || bname == "higamma"
    br = opts.ob_band_highgamma; useMid = true;
else
    error('Unknown band: %s', bname);
end
end

function x = ra_trial_scalar_in_window(B, t_rel, t1, t2)
% B: nTrials x nT. t1/t2: nTrials x 1 (trial-specific).
x = nan(size(B,1),1);
for k = 1:size(B,1)
    if ~isfinite(t1(k)) || ~isfinite(t2(k)) || t2(k) <= t1(k), continue, end
    m = (t_rel >= t1(k)) & (t_rel <= t2(k));
    if sum(m) < 5, continue, end
    x(k) = mean(B(k,m), 'omitnan');
end
x = x(isfinite(x));
end

function [winDefs, w1, w2] = ra_get_window_edges_from_EvRel(EvRel)
winDefs = {'stimon_to_stimoff','stimoff_to_arr','arr_to_stop'};
n = numel(EvRel.stimOff_rel);
w1 = nan(n,3); w2 = nan(n,3);
w1(:,1) = 0;               w2(:,1) = EvRel.stimOff_rel;
w1(:,2) = EvRel.stimOff_rel; w2(:,2) = EvRel.arr_rel;
w1(:,3) = EvRel.arr_rel;     w2(:,3) = EvRel.stop_rel;
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

[~,~,r] = unique(s);
r = double(r);

n0 = numel(xRef);
n1 = numel(xTar);
R1 = sum(r(n0+1:end));
U1 = R1 - n1*(n1+1)/2;
M.AUROC = U1 / (n0*n1);

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
[FF,ord] = sort(FPR);
TT = TPR(ord);
M.rocFPR = FF;
M.rocTPR = TT;
end

function ra_spread_two(xRef, xTar)
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
