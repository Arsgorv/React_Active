function RAA_scatter_pupil_vs_ob_theta_AUROC(RESob, BRES, opts)
% Scatter: each point is a session.
% X = OB theta AUROC (Tar vs Ref)
% Y = pupil AUROC (Tar vs Ref)

if nargin < 3, opts = struct(); end
if ~isfield(opts,'band'),   opts.band = 'theta'; end
if ~isfield(opts,'window'), opts.window = 'stimoff_to_arrival'; end
if ~isfield(opts,'minN'),   opts.minN = 10; end     % min trials per class for both metrics
if ~isfield(opts,'useSpearman'), opts.useSpearman = true; end

S = RESob.Sess;
if isempty(S) || isempty(BRES.T)
    error('Empty RESob.Sess or BRES.T');
end

% filter OB theta/window
Ss = S(strcmpi(string(S.band), string(opts.band)) & strcmpi(string(S.window), string(opts.window)), :);
Ss.session = string(Ss.session);

Tb = BRES.T(strcmpi(string(BRES.T.window), string(opts.window)), :);
Tb.session = string(Tb.session);

% join
T = outerjoin(Ss, Tb, 'Keys','session', 'MergeKeys', true);

% enforce minN on both sides (OB and pupil)
ok = isfinite(T.AUROC_Ss) & isfinite(T.AUROC_Tb);
if ismember('nTar', T.Properties.VariableNames) && ismember('nRef', T.Properties.VariableNames)
    ok = ok & (T.nTar_Ss >= opts.minN) & (T.nRef_Ss >= opts.minN) & (T.nTar_Tb >= opts.minN) & (T.nRef_Tb >= opts.minN);
end
T = T(ok,:);

x = T.AUROC_Ss;  % OB theta AUROC
y = T.AUROC_Tb;  % pupil AUROC

figure('Color','w','Units','pixels','Position',[200 120 650 520], ...
    'Name', sprintf('PupilAUROC_vs_OB_%s_AUROC_%s', opts.band, opts.window));
plot(x, y, 'k.', 'MarkerSize', 14); hold on
xline(0.5,'k:'); yline(0.5,'k:');
xlabel(sprintf('OB %s AUROC (%s)', opts.band, opts.window), 'Interpreter','none');
ylabel(sprintf('Pupil AUROC (%s)', opts.window), 'Interpreter','none');
axis square
xlim([0 1]); ylim([0 1]);
box off

% correlation
if numel(x) >= 3
    if opts.useSpearman
        [rho,p] = corr(x, y, 'Type','Spearman', 'Rows','complete');
        txt = sprintf('Spearman rho=%.2f, p=%.3g, n=%d', rho, p, numel(x));
    else
        [r,p] = corr(x, y, 'Type','Pearson', 'Rows','complete');
        txt = sprintf('Pearson r=%.2f, p=%.3g, n=%d', r, p, numel(x));
    end
else
    txt = sprintf('n=%d (too few for corr)', numel(x));
end
title(txt, 'Interpreter','none');

% optional: label points
if isfield(opts,'labelPoints') && opts.labelPoints
    for i = 1:height(T)
        text(x(i)+0.01, y(i), char(T.session(i)), 'FontSize', 7, 'Interpreter','none');
    end
end

hold off
end
