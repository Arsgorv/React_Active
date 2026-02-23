function plot_ob_tar_ref_boxplots(RES_trialmetrics, bandName, winName)
% Uses RES_trialmetrics.Tall from RAA_collect_ob_trialmetrics_across_sessions

T = RES_trialmetrics.Tall;
bandName = string(bandName);
winName  = string(winName);

sessions = unique(T.session,'stable');

figure('Color','w','Units','pixels','Position',[50 50 1600 700]);
nS = numel(sessions);

for si = 1:nS
    sess = sessions(si);
    Ts = T(T.session==sess & T.band==bandName & T.window==winName,:);
    if isempty(Ts), continue, end

    xTar = Ts.meanPower(Ts.isTar);
    xRef = Ts.meanPower(Ts.isRef);

    subplot(ceil(nS/6), 6, si);
    A = {xTar(:), xRef(:)};
    Leg = {'Tar','Ref'};
    try
        MakeSpreadAndBoxPlot3_SB(A, {}, [1 2], Leg, 'newfig',0, 'paired',0, 'showpoints',1, 'showsigstar','none');
    catch
        plot(ones(size(xTar)), xTar, '.'); hold on
        plot(2*ones(size(xRef)), xRef, '.'); hold off
        set(gca,'XTick',[1 2],'XTickLabel',Leg);
    end
    title(char(sess),'Interpreter','none');
    yline(0,'k:');
    box off
end

sgtitle(sprintf('OB %s | %s | Tar vs Ref (per session)', char(bandName), char(winName)), 'Interpreter','none');
end
