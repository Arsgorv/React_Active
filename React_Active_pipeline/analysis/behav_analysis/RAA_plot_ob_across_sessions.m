function RAA_plot_ob_across_sessions(RES, opts)
% R2018b-compatible plotting:
% (f) heatmaps: session rows, windows columns, per band and group
% (g) AUROC dynamics across sessions (tar vs ref)

if nargin < 2, opts = struct(); end

H = RES.Heat;
S = RES.Sess;

if isempty(S)
    warning('No session summary to plot.');
    return
end

bands = H.bands;
sessions = H.sessions;
winOrder = H.winOrder;

gNames = {'all','tar','ref','nosound','nomotor'};

% ---- heatmaps (one figure per band, 5 subplots)
for bi = 1:numel(bands)
    b = bands(bi);
    figure('Color','w','Units','pixels','Position',[50 50 1700 650], 'Name', sprintf('Heat_%s', char(b)));

    for gi = 1:numel(gNames)
        subplot(1,5,gi);
        M = H.(char(b)).(gNames{gi});
        imagesc(M);
        axis xy
        title(sprintf('%s | %s', char(b), gNames{gi}), 'Interpreter','none');
        set(gca,'XTick',1:numel(winOrder),'XTickLabel',cellstr(winOrder));
        xtickangle(25);
        set(gca,'YTick',1:numel(sessions),'YTickLabel',cellstr(sessions));
        colorbar
        box off
    end
end

% ---- AUROC dynamics (tar vs ref) per band x window
bandsU = unique(S.band);
winsU  = unique(S.window);

figure('Color','w','Units','pixels','Position',[50 50 1700 900], 'Name', 'TarRef_AUROC_dynamics');

nR = numel(bandsU);
nC = numel(winsU);

for bi = 1:nR
    for wi = 1:nC
        subplot(nR, nC, (bi-1)*nC + wi);
        Ss = S(S.band==bandsU(bi) & S.window==winsU(wi),:);

        [~,ord] = ismember(Ss.session, RES.SessInfo.session);
        [~,ix] = sort(ord);
        Ss = Ss(ix,:);

        plot(Ss.AUROC, 'LineWidth', 1); hold on
        yline(0.5,'k:');
        ylim([0 1]);
        xlabel('session #'); ylabel('AUROC');
        title(sprintf('%s | %s', char(bandsU(bi)), char(winsU(wi))), 'Interpreter','none');
        box off
    end
end

% ---- Dynamics across sessions (AUROC / CohenD / meanDiff)
metrics = {'AUROC','CohenD','meanDiff'};
bandsU = unique(S.band);
winsU  = unique(S.window);

for mi = 1:numel(metrics)
    metName = metrics{mi};

    figure('Color','w','Units','pixels','Position',[50 50 1700 900], ...
        'Name', sprintf('TarRef_%s_dynamics', metName));

    nR = numel(bandsU);
    nC = numel(winsU);

    for bi = 1:nR
        for wi = 1:nC
            subplot(nR, nC, (bi-1)*nC + wi);

            Ss = S(S.band==bandsU(bi) & S.window==winsU(wi),:);
            if isempty(Ss), axis off; continue, end

            % order by sessions list
            [~,ord] = ismember(Ss.session, RES.SessInfo.session);
            [~,ix] = sort(ord);
            Ss = Ss(ix,:);

            y = Ss.(metName);
            plot(y, 'LineWidth', 1); hold on

            if strcmpi(metName,'AUROC')
                yline(0.5,'k:');
                ylim([0 1]);
            else
                yline(0,'k:');
            end

            xlabel('session #');
            ylabel(metName);
            title(sprintf('%s | %s', char(bandsU(bi)), char(winsU(wi))), 'Interpreter','none');
            box off
        end
    end
end


end
