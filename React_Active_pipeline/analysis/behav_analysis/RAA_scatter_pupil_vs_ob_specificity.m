function RAA_scatter_pupil_vs_ob_specificity(RESob, BRES, opts)
% RAA_scatter_pupil_vs_ob_specificity
% (1) Adds regression line (+ optional CI band)
% (2) Adds partial correlation controlling for session index
% (3) 4-panel specificity: delta/theta/gamma/highgamma vs pupil
%
% Inputs:
%   RESob.Sess : table with per-session metrics for OB, produced by RAA_collect_ob_trialmetrics_across_sessions
%               Required columns: session, band, window, AUROC, nTar, nRef
%   BRES       : behaviour across-session results. Must contain per-session AUROC for pupil in matching window.
%               Supported formats:
%                 - BRES.Sess (table) with columns: session, bodypart, marker, subset, window, AUROC
%                 - BRES.Tall (table) same idea
%                 - BRES.T (table) same idea
%   opts:
%     .windowOB      = 'stimoff_to_arrival'
%     .windowPupil   = 'stimoff_to_arrival'
%     .subset        = 'regular'
%     .pupilBodypart = 'Pupil-EyeCam'
%     .pupilMarker   = 'pupil_area_007'
%     .bands         = {'delta','theta','gamma','highgamma'}
%     .useSpearman   = false (Pearson default)
%     .minSessN      = 10
%     .minTrials     = 5 (min nTar and nRef per session in OB)
%     .showCI        = false (set true if you want CI band; uses bootstrapping)
%     .nBoot         = 1000
%
% Notes:
% - "partial correlation controlling for session index" uses session order after aligning sessions between BRES and RESob.
% - Session index here is *within the provided sessions list order* (or sorted by session string if ambiguous).

if nargin < 3, opts = struct(); end
if ~isfield(opts,'windowOB'),      opts.windowOB = 'stimoff_to_arrival'; end
if ~isfield(opts,'windowPupil'),   opts.windowPupil = 'stimoff_to_arrival'; end
if ~isfield(opts,'subset'),        opts.subset = 'regular'; end
if ~isfield(opts,'pupilBodypart'), opts.pupilBodypart = 'Pupil-EyeCam'; end
if ~isfield(opts,'pupilMarker'),   opts.pupilMarker = 'pupil_area_007'; end
if ~isfield(opts,'bands'),         opts.bands = {'delta','theta','gamma','highgamma'}; end
if ~isfield(opts,'useSpearman'),   opts.useSpearman = false; end
if ~isfield(opts,'minSessN'),      opts.minSessN = 10; end
if ~isfield(opts,'minTrials'),     opts.minTrials = 5; end
if ~isfield(opts,'showCI'),        opts.showCI = false; end
if ~isfield(opts,'nBoot'),         opts.nBoot = 1000; end

% -------------------- pull pupil AUROC per session --------------------
Tpup = [];
if isstruct(BRES)
    if isfield(BRES,'Sess') && istable(BRES.Sess), Tpup = BRES.Sess; end
    if isempty(Tpup) && isfield(BRES,'Tall') && istable(BRES.Tall), Tpup = BRES.Tall; end
    if isempty(Tpup) && isfield(BRES,'T') && istable(BRES.T), Tpup = BRES.T; end
elseif istable(BRES)
    Tpup = BRES;
end
if isempty(Tpup)
    error('BRES does not contain a recognizable table (Sess/Tall/T).');
end

% Normalize string columns if present
if ismember('session', Tpup.Properties.VariableNames), Tpup.session = string(Tpup.session); end
if ismember('bodypart', Tpup.Properties.VariableNames), Tpup.bodypart = string(Tpup.bodypart); end
if ismember('marker', Tpup.Properties.VariableNames), Tpup.marker = string(Tpup.marker); end
if ismember('subset', Tpup.Properties.VariableNames), Tpup.subset = string(Tpup.subset); end
if ismember('window', Tpup.Properties.VariableNames), Tpup.window = string(Tpup.window); end

needCols = {'session','AUROC'};
for i = 1:numel(needCols)
    if ~ismember(needCols{i}, Tpup.Properties.VariableNames)
        error('Pupil table missing required column: %s', needCols{i});
    end
end

% Filter pupil metric row(s)
mP = true(height(Tpup),1);
if ismember('bodypart', Tpup.Properties.VariableNames)
    mP = mP & (lower(Tpup.bodypart) == lower(string(opts.pupilBodypart)));
end
if ismember('marker', Tpup.Properties.VariableNames)
    mP = mP & (lower(Tpup.marker) == lower(string(opts.pupilMarker)));
end
if ismember('subset', Tpup.Properties.VariableNames)
    mP = mP & (lower(Tpup.subset) == lower(string(opts.subset)));
end
if ismember('window', Tpup.Properties.VariableNames)
    mP = mP & (lower(Tpup.window) == lower(string(opts.windowPupil)));
end

Tp = Tpup(mP, {'session','AUROC'});
Tp.Properties.VariableNames = {'session','pupilAUROC'};
Tp.session = string(Tp.session);

% If duplicates per session remain, collapse by median
% Ensure Tp is not empty
if isempty(Tp) || height(Tp)==0
    error('Tp is empty after filtering. Check subset/window selection for pupil table.');
end

% Ensure AUROC is numeric
Tp.pupilAUROC = double(Tp.pupilAUROC);

% Recompute groups on the exact Tp currently used
[Gp, sessKey] = findgroups(Tp.session);

% splitapply requires full positive integer group vector
Gp = full(double(Gp));

pAU = splitapply(@(x) median(x,'omitnan'), Tp.pupilAUROC, Gp);
Tp = table(sessKey, pAU, 'VariableNames', {'session','pupilAUROC'});

% -------------------- pull OB AUROC per session x band --------------------
S = RESob.Sess;
if ~istable(S), error('RESob.Sess must be a table.'); end

% Normalize strings
if ismember('session', S.Properties.VariableNames), S.session = string(S.session); end
if ismember('band', S.Properties.VariableNames),    S.band = string(S.band); end
if ismember('window', S.Properties.VariableNames),  S.window = string(S.window); end

req = {'session','band','window','AUROC'};
for i = 1:numel(req)
    if ~ismember(req{i}, S.Properties.VariableNames)
        error('RESob.Sess missing required column: %s', req{i});
    end
end

% Ensure nTar/nRef exist for filtering (optional but recommended)
hasNR = ismember('nTar', S.Properties.VariableNames) && ismember('nRef', S.Properties.VariableNames);

% Filter window
mW = (lower(S.window) == lower(string(opts.windowOB)));
Sw = S(mW,:);

% -------------------- 4-panel plot --------------------
bands = string(opts.bands(:)');
nb = numel(bands);

figure('Color','w','Units','pixels','Position',[60 60 1200 900], 'Name', 'Pupil_vs_OB_specificity');

for bi = 1:nb
    bname = bands(bi);

    % filter band
    mb = (lower(Sw.band) == lower(bname));
    Sb = Sw(mb, {'session','AUROC'});
    Sb.Properties.VariableNames = {'session','obAUROC'};
    Sb.session = string(Sb.session);

    if hasNR
        SbFull = Sw(mb,:);
        keep = (SbFull.nTar >= opts.minTrials) & (SbFull.nRef >= opts.minTrials);
        Sb = Sb(keep,:);
    end

    % align sessions with pupil
    [T, sessIdx] = ra_align_two_tables_by_session(Sb, Tp);

    x = T.obAUROC;
    y = T.pupilAUROC;

    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    si = sessIdx(ok); % session index covariate for partial corr

    subplot(2,2,bi); hold on

    % scatter with session-index shading
    % use si (aligned+filtered session index covariate)
    s0 = min(si); s1 = max(si);
    den = max(1, (s1 - s0));
    u = (si - s0) / den; % 0..1 within aligned sessions
    
    base = ra_band_base_color(bname);  % <-- FIX: bname, not bandName
    
    for ii = 1:numel(x)
        if ~isfinite(x(ii)) || ~isfinite(y(ii)) || ~isfinite(u(ii)), continue, end
        ci = ra_color_by_session(base, u(ii));
        plot(x(ii), y(ii), 'o', ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', ci, ...
            'MarkerEdgeColor', ci);
    end
    % reference lines
    xline(0.5,'k:'); yline(0.5,'k:');
    xlim([0.35 0.7])
    ylim([0.35 0.7])
    axis square

    % correlation
    if numel(x) >= 3
        if opts.useSpearman
            [r,p] = corr(x, y, 'Type','Spearman', 'Rows','complete');
            corrName = 'Spearman';
        else
            [r,p] = corr(x, y, 'Type','Pearson', 'Rows','complete');
            corrName = 'Pearson';
        end
    else
        r = nan; p = nan; corrName = ternary(opts.useSpearman,'Spearman','Pearson');
    end

    % regression line (least squares)
    if numel(x) >= 2
        [b0,b1] = ra_fit_line_ols(x,y); % y = b0 + b1*x
        xx = linspace(min(x), max(x), 100);
        yy = b0 + b1*xx;
        plot(xx, yy, 'k-', 'LineWidth', 1.8);

        if opts.showCI
            [lo,hi] = ra_boot_ci_line(x,y,xx,opts.nBoot);
            patch([xx fliplr(xx)], [lo fliplr(hi)], [0 0 0], 'FaceAlpha',0.08, 'EdgeColor','none');
        end
    end

    % partial correlation controlling for session index
    rp = nan; pp = nan;
    if numel(x) >= 5
        if opts.useSpearman
            % Spearman partial corr: rank-transform x,y then partial corr (Pearson) on ranks
            xr = tiedrank(x);
            yr = tiedrank(y);
            [rp, pp] = partialcorr(xr, yr, si, 'Rows','complete', 'Type','Pearson');
            pcName = 'partial(Spearman ranks | sess#)';
        else
            [rp, pp] = partialcorr(x, y, si, 'Rows','complete', 'Type','Pearson');
            pcName = 'partial(Pearson | sess#)';
        end
    else
        pcName = 'partial(| sess#)';
    end

    xlabel(sprintf('OB %s AUROC (%s)', char(bname), opts.windowOB), 'Interpreter','none');
    ylabel(sprintf('Pupil AUROC (%s)', opts.windowPupil), 'Interpreter','none');

    title(sprintf('%s | n=%d', char(bname), numel(x)), 'Interpreter','none');

    % Put stats in text box
    txt = sprintf('%s r=%.2f, p=%.3g\n%s r=%.2f, p=%.3g', corrName, r, p, pcName, rp, pp);
    yl = ylim; xl = xlim;
    text(xl(1)+0.02*(xl(2)-xl(1)), yl(2)-0.06*(yl(2)-yl(1)), txt, 'FontSize', 10, 'VerticalAlignment','top');

    box off
    hold off
end

sgtitle(sprintf('Pupil vs OB AUROC | subset=%s | window=%s', opts.subset, opts.windowOB), 'Interpreter','none');

end

% ======================= helpers =======================

function [T, sessIdx] = ra_align_two_tables_by_session(Tob, Tpup)
% Tob: columns {session, obAUROC}
% Tpup: columns {session, pupilAUROC}
Tob.session = string(Tob.session);
Tpup.session = string(Tpup.session);

% Use intersection; preserve "training order" by sorting session strings
sess = intersect(Tob.session, Tpup.session);
sess = sort(sess); % safe default; replace with your sessions list order if you have it

% Session index covariate
sessIdx = (1:numel(sess))';

% build aligned table
x = nan(numel(sess),1);
y = nan(numel(sess),1);

for i = 1:numel(sess)
    s = sess(i);
    a = Tob.obAUROC(Tob.session==s);
    b = Tpup.pupilAUROC(Tpup.session==s);
    x(i) = median(a,'omitnan');
    y(i) = median(b,'omitnan');
end

T = table(sess, x, y, 'VariableNames', {'session','obAUROC','pupilAUROC'});

end

function [b0,b1] = ra_fit_line_ols(x,y)
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 2
    b0 = nan; b1 = nan; return
end
X = [ones(numel(x),1) x];
b = X \ y;
b0 = b(1); b1 = b(2);
end

function [lo,hi] = ra_boot_ci_line(x,y,xx,nBoot)
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);

n = numel(x);
Yhat = nan(nBoot, numel(xx));
for b = 1:nBoot
    ii = randi(n, n, 1);
    xb = x(ii); yb = y(ii);
    [b0,b1] = ra_fit_line_ols(xb,yb);
    Yhat(b,:) = b0 + b1*xx;
end
lo = prctile(Yhat, 2.5, 1);
hi = prctile(Yhat, 97.5, 1);
end

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end

function c = ra_color_by_session(baseRGB, u)
% u in [0..1], early=0 -> light, late=1 -> dark
u = max(0, min(1, u));
a = 1 - u;               % early: a=1 (white), late: a=0 (base)
c = baseRGB*(1-a) + [1 1 1]*a;
c = max(0, min(1, c));
end

function base = ra_band_base_color(bname)
bn = lower(string(bname));
switch bn
    case "delta"
        base = [0.10 0.35 0.90]; % blue
    case "theta"
        base = [0.45 0.15 0.70]; % purple
    case "gamma"
        base = [0.10 0.60 0.20]; % green
    case "highgamma"
        base = [0.80 0.10 0.10]; % red
    otherwise
        base = [0.20 0.20 0.20]; % fallback
end
end