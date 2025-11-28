function pupil_gamma_performance_corr(datapath, sm_win_behav, sm_win_ephys)
% Plot z-scored pupil & OB-gamma (raw + ultra-slow band), plus sliding hit-rate;
% then compute a 3×3 correlation matrix across trials (gamma_slow, pupil_slow, hit_MA).
%
% win = moving-average window in TRIALS for hit-rate.
%

%% Options
fig_visibility = 'off';

figdir = fullfile(datapath,'Figures');
if ~exist(figdir,'dir'), mkdir(figdir); end

%% Load
load(fullfile(datapath, 'stim', 'trial_structure.mat'))
load(fullfile(datapath, 'video', 'DLC_data.mat'), 'pupil_area_007', 'pupil_center_007_mvt')

if exist(fullfile(datapath,'SleepScoring_OBGamma.mat'),'file')
    S = load(fullfile(datapath, 'SleepScoring_OBGamma.mat'),'BrainPower');
else
    S = struct;
end
if isfield(S,'BrainPower') && ~isempty(S.BrainPower) && numel(S.BrainPower.Power)>=1
    SmoothGamma = S.BrainPower.Power{1};
else
    sm_w = .1; % seconds
    load([datapath filesep 'ChannelsToAnalyse/Bulb_deep.mat'])              %#ok<LOAD>
    LFP_load = load(strcat([datapath,'/LFPData/LFP',num2str(channel),'.mat'])); %#ok<NODEF>
    LFP = LFP_load.LFP;                                                     %#ok<NODEF>
    FilGamma = FilterLFP(LFP,[40 60],1024);
    tEnveloppeGamma = tsd(Range(LFP), abs(hilbert(Data(FilGamma))) );
    SmoothGamma = tsd(Range(tEnveloppeGamma), runmean(Data(tEnveloppeGamma), ...
        ceil(sm_w/median(diff(Range(tEnveloppeGamma,'s'))))));
    close; clear LFP_load LFP FilGamma tEnveloppeGamma
    BrainPower.Power{1} = SmoothGamma; %#ok<NASGU>
    save(fullfile(datapath, 'SleepScoring_OBGamma.mat'), 'BrainPower', '-append');
end

%% Resample gamma to pupil timeline
resamp = @(tsdObj,refTsd) tsd( Range(refTsd), ...
    interp1( Range(tsdObj), Data(tsdObj), Range(refTsd), 'linear','extrap') );

Gamma_rs = resamp(SmoothGamma, pupil_area_007);
t = Range(pupil_area_007);

%% Z-score raw
Gamma_z_cont = tsd(t, zscore(runmean(Data(Gamma_rs), max(1,sm_win_ephys))));

hit_raw = double(trial_structure.hit_d(:));
tTrials_ts = trial_structure.trial_onset.Target;
tTrials = Range(tTrials_ts);
hitMA = runmean(hit_raw, max(1,sm_win_behav));
Hit_z = tsd(tTrials, zscore(hitMA));

%% Sliding hit-rate over TRIALS
tRewards = tTrials + trial_structure.reward_onset.Target*1e4;

a = tTrials(:);
dur = tRewards(:) - a + 2*1e4;
dur = dur(dur>0 & isfinite(dur));
b = tTrials(:) + median(dur);

trial_epoch = intervalSet(trial_structure.trial_onset.Target,ts(Range(trial_structure.trial_onset.Target)+dur));
% pupil_wtrial = Restrict(pupil_area_007, trial_epoch);
% gamma_wtrial = Restrict(Gamma_rs, trial_epoch);

% Pre-trial baseline window (?2?0 s)
a_pre = a - 2*1e4;  b_pre = a;

%% Run the full analysis twice: AREA and POS
out_all = struct();

for MOD = 1:2
    if MOD==1
        modeName   = 'AREA';
        pupil_raw  = tsd(t, Data(pupil_area_007));
        pupil_lab  = 'pupil area';
        suffix     = '_AREA';
    else
        modeName   = 'MVT';
        % radial displacement from median center (robust to camera drift)
        C  = Data(pupil_center_007_mvt);
        if size(C,2)>=2
            medC = [median(C(:,1),'omitnan') median(C(:,2),'omitnan')];
            mvtAmp = sqrt((C(:,1)-medC(1)).^2 + (C(:,2)-medC(2)).^2);
        else
            mvtAmp = C; % if center is 1D, just use it
        end
        pupil_raw = tsd(t, mvtAmp); 
        pupil_lab = 'pupil movement';
        suffix    = '_MVT';
    end
    % continuous Z-scored pupil
    Pupil_z_cont = tsd(t, zscore(runmean(Data(pupil_raw), max(1,sm_win_ephys))));
    
    %% Within-trial & pre-trial features (means from continuous Z traces)
    % within-trial
    tp = Range(Pupil_z_cont);  xp = Data(Pupil_z_cont);
    tg = Range(Gamma_z_cont);  xg = Data(Gamma_z_cont);
    
    pupil_mean_trial = nan(numel(a),1);
    gamma_mean_trial = nan(numel(a),1);
    for i = 1:numel(a)
        wP = tp>=a(i) & tp<=b(i);
        if any(wP), pupil_mean_trial(i) = mean(xp(wP),'omitnan'); end
        wG = tg>=a(i) & tg<=b(i);
        if any(wG), gamma_mean_trial(i) = mean(xg(wG),'omitnan'); end
    end
    
    % pre-trial baseline (?2?0 s)
    pupil_pre = nan(numel(a),1);
    gamma_pre = nan(numel(a),1);
    for i = 1:numel(a)
        wP = tp>=a_pre(i) & tp<=b_pre(i);
        if any(wP), pupil_pre(i) = mean(xp(wP),'omitnan'); end
        wG = tg>=a_pre(i) & tg<=b_pre(i);
        if any(wG), gamma_pre(i) = mean(xg(wG),'omitnan'); end
    end
    
    % trial-MA on within-trial features (long-scale vigilance)
    tEndTrials = a + median(dur);
    Pupil_smoothed = tsd(tEndTrials, runmean(pupil_mean_trial, max(1,sm_win_behav)));
    Gamma_smoothed = tsd(tEndTrials, runmean(gamma_mean_trial, max(1,sm_win_behav)));
    
    %% PLOT: continuous + trial-MA
    f = figure('visible', fig_visibility, 'Color','w','Position',[100 80 1500 520]);
    subplot(211),  hold on
    plot(Range(Gamma_z_cont, 'min'), Data(Gamma_z_cont), 'Color',[0.5 0.5 1], 'DisplayName','OB gamma');
    plot(Range(Pupil_z_cont, 'min'), Data(Pupil_z_cont), 'k','DisplayName',pupil_lab);
    plot(Range(Hit_z,'min'), Data(Hit_z), 'LineWidth',2, 'DisplayName',sprintf('hit-rate MA %d',sm_win_behav));
    xlabel('Time (min)'); ylabel('z-score'); box off; xlim([ t(1)/1e4/60 t(end)/1e4/60]); try, makepretty; end
    title([upper(modeName) ' | Pupil & OB gamma with sliding hit-rate (trials)'])
    legend('Location','best');
        
    subplot(212),  hold on
    plot(Range(Pupil_smoothed, 'min'), Data(Pupil_smoothed), '.-k','DisplayName',[pupil_lab ' (within-trial MA)']);
    plot(Range(Gamma_smoothed, 'min'), Data(Gamma_smoothed), '.-','Color',[0.5 0.5 1], 'DisplayName','OB gamma (within-trial MA)');
    plot(Range(Hit_z,'min'), Data(Hit_z), '.-','DisplayName',sprintf('hit-rate MA %d',sm_win_behav));
    xlabel('Time (min)'); ylabel('z-score'); box off; xlim([ t(1)/1e4/60 t(end)/1e4/60]); try, makepretty; end
    title(['smoothed with ' num2str(sm_win_behav) ' trials'])
    legend('Location','best');
        
    saveas(f, fullfile(figdir,['session_timeseries' suffix '.png']));
    try, saveas(f, fullfile(figdir,['session_timeseries' suffix '.svg'])); end
    
    %% Trial-wise matrices: within-trial & MA 
    Xtrial = [zscore(pupil_mean_trial), zscore(gamma_mean_trial), zscore(hitMA(:))];
    n = min([size(Xtrial,1), numel(hit_raw)]);
    X = Xtrial(1:n,:);
    m1 = all(isfinite(X),2);
    X = X(m1,:);
    labels = {pupil_lab,'OB gamma power','hit rate'};
    
    [R,P] = corr(X, 'rows','pairwise', 'type','Pearson');
    
    f2 = figure('visible', fig_visibility, 'Color','w','Position',[200 120 450 420]);
    try, imagesc(R, [-1 1]); colormap('redblue'); catch, imagesc(R, [0 1]); colormap(parula); end
    axis square; colorbar; caxis([-1 1]);
    set(gca,'XTick',1:3,'XTickLabel',labels,'YTick',1:3,'YTickLabel',labels,'XTickLabelRotation',20);
    title(['Pair-wise correlations (within-trial means) | ' modeName]); box on
    for i=1:3, for j=1:3
            val = R(i,j);
            text(j,i,sprintf('%.2f', val),'HorizontalAlignment','center','FontWeight','bold', ...
                'Color', tern(abs(val)>0.6,'w','k'));
        end, end
    saveas(f2, fullfile(figdir,['session_corr_within' suffix '.png']));
    try, saveas(f2, fullfile(figdir,['session_corr_within' suffix '.svg'])); end
    
    % Trial-MA correlations (long-scale vigilance)
    Xtrial_smoothed = [Data(Pupil_smoothed), Data(Gamma_smoothed), hitMA(:)];
    X_sm = Xtrial_smoothed(1:n,:);
    m2 = all(isfinite(X_sm),2);
    X_sm = X_sm(m2,:);
    [R_sm,P_sm] = corr(X_sm, 'rows','pairwise', 'type','Pearson');
    
    f3 = figure('visible', fig_visibility, 'Color','w','Position',[200 120 450 420]);
    try, imagesc(R_sm, [-1 1]); colormap('redblue'); catch, imagesc(R_sm, [0 1]); colormap(parula); end
    axis square; colorbar; caxis([-1 1]);
    set(gca,'XTick',1:3,'XTickLabel',labels,'YTick',1:3,'YTickLabel',labels,'XTickLabelRotation',20);
    title(['Smoothed pair-wise correlations (within-trial MA) | ' modeName]); box on
    for i=1:3, for j=1:3
            val = R_sm(i,j);
            text(j,i,sprintf('%.2f', val),'HorizontalAlignment','center','FontWeight','bold', ...
                'Color', tern(abs(val)>0.6,'w','k'));
        end, end
    saveas(f3, fullfile(figdir,['session_corr_MA' suffix '.png']));
    try, saveas(f3, fullfile(figdir,['session_corr_MA' suffix '.svg'])); end
    
    %% Partial correlations (within-trial & MA)
    % within-trial partials
    pr_gamma_hit_given_pupil = partialcorr(X(:,2), X(:,3), X(:,1), 'rows','pairwise');
    pr_pupil_hit_given_gamma = partialcorr(X(:,1), X(:,3), X(:,2), 'rows','pairwise');
    % MA partials
    if ~isempty(X_sm)
        pr_gamma_hit_given_pupil_MA = partialcorr(X_sm(:,2), X_sm(:,3), X_sm(:,1), 'rows','pairwise');
        pr_pupil_hit_given_gamma_MA = partialcorr(X_sm(:,1), X_sm(:,3), X_sm(:,2), 'rows','pairwise');
    else
        pr_gamma_hit_given_pupil_MA = NaN; pr_pupil_hit_given_gamma_MA = NaN;
    end
    
    %% Pre-trial baseline correlations (and partials)
    Xpre = [zscore(pupil_pre), zscore(gamma_pre), zscore(hitMA(:))];
    Xpre = Xpre(1:n,:); mp = all(isfinite(Xpre),2); Xpre = Xpre(mp,:);
    [R_pre,P_pre] = corr(Xpre, 'rows','pairwise', 'type','Pearson');
    pr_gamma_hit_given_pupil_pre = partialcorr(Xpre(:,2), Xpre(:,3), Xpre(:,1), 'rows','pairwise');
    pr_pupil_hit_given_gamma_pre = partialcorr(Xpre(:,1), Xpre(:,3), Xpre(:,2), 'rows','pairwise');
    
    f_pre = figure('visible', fig_visibility, 'Color','w','Position',[200 120 450 420]);
    try, imagesc(R_pre, [-1 1]); colormap('redblue'); catch, imagesc(R_pre, [0 1]); colormap(parula); end
    axis square; colorbar; caxis([-1 1]);
    set(gca,'XTick',1:3,'XTickLabel',{[pupil_lab ' (pre)'],'OB gamma (pre)','hit rate'},'YTick',1:3,'YTickLabel',{[pupil_lab ' (pre)'],'OB gamma (pre)','hit rate'},'XTickLabelRotation',20);
    title(['Pre-trial baseline correlations | ' modeName]); box on
    for i=1:3, for j=1:3
            val = R_pre(i,j);
            text(j,i,sprintf('%.2f', val),'HorizontalAlignment','center','FontWeight','bold', ...
                'Color', tern(abs(val)>0.6,'w','k'));
        end, end
    saveas(f_pre, fullfile(figdir,['session_corr_PRE' suffix '.png']));
    try, saveas(f_pre, fullfile(figdir,['session_corr_PRE' suffix '.svg'])); end
    
    %% Logistic regression + 5-fold CV AUC (within-trial features)
    y = hit_raw(1:n);  y = y(m1);
    feat = zscore(X(:,1:2));
    
    has_two_classes = numel(unique(y(~isnan(y)))) == 2;
    var_ok = all(std(feat,0,1) > eps);
    K = 5;
    
    AUCp = nan(K,1); AUCg = nan(K,1); AUCb = nan(K,1);
    BETA = struct('pupil',[],'gamma',[],'both',[]);
    
    if has_two_classes && var_ok && numel(y) >= max(10, 2*K)
        try, C = cvpartition(y,'KFold',K); catch, C = cvpartition(length(y),'KFold',K); end
        for k = 1:K
            idxTr = training(C,k); idxTe = test(C,k);
            if numel(unique(y(idxTr))) < 2 || numel(unique(y(idxTe))) < 2, continue, end
            
            Xtr = feat(idxTr,:); ytr = y(idxTr);
            Xte = feat(idxTe,:); yte = y(idxTe);
            
            keepCols = std(Xtr,0,1) > eps;
            
            if keepCols(1)
                bp = glmfit(Xtr(:,1), ytr, 'binomial','link','logit');
                sp = glmval(bp, Xte(:,1), 'logit');
                AUCp(k) = local_auc(yte, sp);
            end
            if keepCols(2)
                bg = glmfit(Xtr(:,2), ytr, 'binomial','link','logit');
                sg = glmval(bg, Xte(:,2), 'logit');
                AUCg(k) = local_auc(yte, sg);
            end
            if all(keepCols) && abs(corr(Xtr(:,1),Xtr(:,2),'rows','complete'))<0.999
                bb = glmfit(Xtr, ytr, 'binomial','link','logit');
                sb = glmval(bb, Xte, 'logit');
                AUCb(k) = local_auc(yte, sb);
            end
        end
        
        if std(feat(:,1))>eps, BETA.pupil = glmfit(feat(:,1), y, 'binomial','link','logit'); end
        if std(feat(:,2))>eps, BETA.gamma = glmfit(feat(:,2), y, 'binomial','link','logit'); end
        if all(std(feat,0,1)>eps) && abs(corr(feat(:,1),feat(:,2),'rows','complete'))<0.999
            BETA.both  = glmfit(feat, y, 'binomial','link','logit');
        end
    end
    
    % AUC box
    f4 = figure('visible', fig_visibility, 'Color','w','Position',[140 100 260 320]); hold on
    hasMS = exist('MakeSpreadAndBoxPlot3_SB','file')==2;
    Acell = {AUCp(isfinite(AUCp)), AUCg(isfinite(AUCg)), AUCb(isfinite(AUCb))};
    if hasMS && any(cellfun(@(v)~isempty(v), Acell))
        MakeSpreadAndBoxPlot3_SB(Acell, {}, [], {'pupil','gamma','both'}, 'ShowSigstar','none', 'paired',0, 'showpoints',1);
    else
        Y = [AUCp, AUCg, AUCb];
        for i=1:3
            boxchart(i*ones(size(Y,1),1), Y(:,i)); hold on
            scatter(i+(rand(size(Y,1),1)-0.5)*0.08, Y(:,i), 12, 'k','filled','MarkerFaceAlpha',0.35);
        end
        set(gca,'XTick',1:3,'XTickLabel',{'pupil','gamma','both'});
    end
    ylabel('CV AUC'); ylim([0.5 1]); title(sprintf('Logistic regression (K=%d) | %s',K,modeName)); box off; grid on
    saveas(f4, fullfile(figdir,['session_auc_box' suffix '.png']));
    try, saveas(f4, fullfile(figdir,['session_auc_box' suffix '.svg'])); end
    
    %% stats vis
    % pupil_mean_trial, gamma_mean_trial, hitMA are already built above
    zp = zscore(pupil_mean_trial(:));
    zg = zscore(gamma_mean_trial(:));
    zh = zscore(hitMA(:));
    
    edges = -3:0.25:3;
    fD = figure('visible', fig_visibility, 'Color','w','Position',[120 120 700 450]); hold on
    
    % histograms (PDF-normalized) + light transparency
    histogram(zp(isfinite(zp)), edges, 'Normalization','pdf', 'FaceColor','k',  'FaceAlpha',0.15, 'EdgeColor','none');
    histogram(zg(isfinite(zg)), edges, 'Normalization','pdf', 'FaceColor',[0.5 0.5 1], 'FaceAlpha',0.20, 'EdgeColor','none');
    histogram(zh(isfinite(zh)), edges, 'Normalization','pdf', 'FaceColor',[1 0.5 0], 'FaceAlpha',0.20, 'EdgeColor','none');
    
    % smooth KDE overlays when available
    try
        [fx,xx] = ksdensity(zp(isfinite(zp))); plot(xx,fx,'k','LineWidth',1.5);
        [fx,xx] = ksdensity(zg(isfinite(zg))); plot(xx,fx,'Color',[0.3 0.3 1],'LineWidth',1.5);
        [fx,xx] = ksdensity(zh(isfinite(zh))); plot(xx,fx,'Color',[1 0.3 0],'LineWidth',1.5);
    end
    
    % medians
    mzp = median(zp,'omitnan'); xline(mzp, '--k',  'LineWidth',1);
    mzg = median(zg,'omitnan'); xline(mzg, '--',   'Color',[0.3 0.3 1], 'LineWidth',1);
    mzh = median(zh,'omitnan'); xline(mzh, '--',   'Color',[1 0.3 0],   'LineWidth',1);
    
    xlabel('Z (within-trial means)'); ylabel('density'); grid on; box off
    title(['Distributions | ' modeName])
    legend({pupil_lab,'OB gamma','hit-rate MA', ...
        [pupil_lab ' KDE'], 'OB gamma KDE', 'hit-rate KDE'}, 'Location','best')
    try, makepretty; end
    saveas(fD, fullfile(figdir, ['session_distributions' suffix '.png']));
    try, saveas(fD, fullfile(figdir, ['session_distributions' suffix '.svg'])); end
    
    %% Two zoom panels: low vs high engagement (~15 trials), with trial patches
    % z-score across the whole session (not per window)
    zH = zscore(hitMA(:));
    zG = zscore(gamma_mean_trial(:));
    % robust sliding means using convolution, ignoring NaNs
    winTrials = 10;
    kbox = ones(winTrials,1);
    zH0 = zH; zH0(~isfinite(zH0)) = 0; vH = double(isfinite(zH));
    zG0 = zG; zG0(~isfinite(zG0)) = 0; vG = double(isfinite(zG));
    sumH = conv(zH0, kbox, 'valid'); nH = conv(vH,kbox, 'valid');
    sumG = conv(zG0, kbox, 'valid'); nG = conv(vG,kbox, 'valid');
    muHit = sumH ./ max(1,nH);
    muGam = sumG ./ max(1,nG);
    okWin = (nH >= ceil(0.8*winTrials)) & (nG >= ceil(0.8*winTrials)); % require coverage
    score = 0.5*muHit + 0.5*muGam;
    score(~okWin) = NaN;
    [~,kHigh] = max(score,[],'omitnan');
    [~,kLow]  = min(score,[],'omitnan');
    if isempty(kHigh) || isnan(kHigh), kHigh = 1; end
    if isempty(kLow)  || isnan(kLow),  kLow  = max(1, numel(score)-winTrials+1); end
    
    fz = figure('visible', fig_visibility, 'Color','w','Position',[80 80 800 420]);
    yl = [-3 3.5];  % same y for both panels (Z units)
%     yl = [-2 2];  % same y for both panels (Z units)

    tpad = 1*1e4;
    
    subplot(1,2,1); hold on
%     yl = ylim;                                     % use actual y-lims
    t0 = a(kLow)-tpad; t1 = b(kLow+winTrials-1)+tpad;

    for k = 1:numel(a)
        patch([a(k) a(k) b(k) b(k)]/60/1e4, ...
            [yl(1) yl(2) yl(2) yl(1)], ...
            [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.3);
    end
    uistack(findobj(gca,'Type','line'),'top');     % keep traces on top
    
    plot(Range(pupil_raw,'min'), Data(Pupil_z_cont), 'k');
    plot(Range(pupil_raw,'min'), Data(Gamma_z_cont), 'Color',[0.5 0.5 1]);
    %     plot(Range(Pupil_smoothed,'min'), Data(Pupil_smoothed), 'k');
    %     plot(Range(Pupil_smoothed,'min'), Data(Gamma_smoothed), 'Color',[0.5 0.5 1]);
    plot(Range(Hit_z,'min'), Data(Hit_z), 'LineWidth',1.5);
    yline(0,'--r', 'LineWidth', 1)
    xlim([t0 t1]/60/1e4); ylim(yl); box off
    title(['LOW engagement | ' modeName]); xlabel('Time (min)'); ylabel('z-score'); makepretty
    
    subplot(1,2,2); hold on
    t0 = a(kHigh)-tpad; t1 = b(kHigh+winTrials-1)+tpad;
    for k=kHigh:(kHigh+winTrials-1)
        patch([a(k) a(k) b(k) b(k)]/60/1e4, [yl(1) yl(2) yl(2) yl(1)], [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.3);
    end
    plot(Range(pupil_raw,'min'), Data(Pupil_z_cont), 'k');
    plot(Range(pupil_raw,'min'), Data(Gamma_z_cont), 'Color',[0.5 0.5 1]);
%     plot(Range(Pupil_smoothed,'min'), Data(Pupil_smoothed), 'k');
%     plot(Range(Pupil_smoothed,'min'), Data(Gamma_smoothed), 'Color',[0.5 0.5 1]);
    plot(Range(Hit_z,'min'), Data(Hit_z), 'LineWidth',1.5);
    yline(0,'--r', 'LineWidth', 1)
    xlim([t0 t1]/60/1e4); ylim(yl); box off
    title(['HIGH engagement | ' modeName]); xlabel('Time (min)'); makepretty
    saveas(fz, fullfile(figdir,['session_zoom_lowHigh' suffix '.png']));
    try, saveas(fz, fullfile(figdir,['session_zoom_lowHigh' suffix '.svg'])); end

    %% Save stats for this modality
    out = struct;
    out.params = struct('sm_win_behav',sm_win_behav,'sm_win_ephys',sm_win_ephys,'Kfold',K);
    out.traces = struct('Pupil_z',Pupil_z_cont,'Gamma_z',Gamma_z_cont,'Hit_z',Hit_z);
    hit_sub = hit_raw(1:n);
    hit_sub = hit_sub(m1);
    out.trial_tbl = table(tTrials(m1), X(:,1), X(:,2), X(:,3), hit_sub, ...  
        'VariableNames', {'t_onset','pupil_withinZ','gamma_withinZ','hit_MA_Z','hit_raw'});
        out.R = R; out.P = P;  out.R_sm = R_sm; out.P_sm = P_sm; out.labels = labels;
    out.R_pre = R_pre; out.P_pre = P_pre;
    out.partial = struct( ...
        'within',  [pr_pupil_hit_given_gamma, pr_gamma_hit_given_pupil], ...
        'MA',      [pr_pupil_hit_given_gamma_MA, pr_gamma_hit_given_pupil_MA], ...
        'pre',     [pr_pupil_hit_given_gamma_pre, pr_gamma_hit_given_pupil_pre] );
    out.cvAUC = struct('pupil',AUCp,'gamma',AUCg,'both',AUCb, ...
        'mean',[mean(AUCp,'omitnan'), mean(AUCg,'omitnan'), mean(AUCb,'omitnan')]);
    out.beta = BETA;
    
    save(fullfile(figdir, ['session_stats_' lower(modeName) '.mat']), 'out');
    out_all.(lower(modeName)) = out;
end
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end

function A = local_auc(y, score)
y = y(:); score = score(:);
if exist('perfcurve','file')==2
    [~,~,~,A] = perfcurve(y, score, 1);
else
    pos = score(y==1); neg = score(y==0);
    if isempty(pos) || isempty(neg), A = NaN; return; end
    r = tiedrank([pos;neg]);
    Rpos = sum(r(1:numel(pos)));
    n1 = numel(pos); n0 = numel(neg);
    U = Rpos - n1*(n1+1)/2;
    A = U / (n1*n0);
end
end

% x = zscore(Data(pupil_area_007));
% pupil_mean_trial = nan(numel(a),1);
% for i = 1:numel(a)
%     w = Range(pupil_area_007)>=a(i) & Range(pupil_area_007)<=b(i);
%     if any(w), pupil_mean_trial(i) = mean(x(w),'omitnan'); end
% end
% 
% x1 = zscore(Data(Gamma_rs));
% gamma_mean_trial = nan(numel(a),1);
% for i = 1:numel(a)
%     w = Range(Gamma_rs)>=a(i) & Range(Gamma_rs)<=b(i);
%     if any(w), gamma_mean_trial(i) = mean(x1(w),'omitnan'); end
% end
% 
% % t_start = Start(trial_epoch);
% % t_end = End(trial_epoch);
% % pupil_data = Data(pupil_wtrial);
% % gamma_data = Data(gamma_wtrial);
% % for tr_n = 1:n
% %     mean_pupil_trial = pupil_data(t_start(tr_n):t_end(tr_n));
% %     mean_gamma_trial = gamma_data(t_start(tr_n):t_end(tr_n));
% % end
% 
% Gamma_smoothed = tsd(End(trial_epoch), runmean(gamma_mean_trial, sm_win_behav));
% Pupil_smoothed = tsd(End(trial_epoch), runmean(pupil_mean_trial, sm_win_behav));
