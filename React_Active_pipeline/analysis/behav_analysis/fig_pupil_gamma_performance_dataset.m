function S = fig_pupil_gamma_performance_dataset(sessionPaths,signalType)
% Aggregate pupil/gamma/performance results across sessions by LOADING the
% per-session 'out' struct saved at <session>/Figures/session_stats.mat.
%
% Outputs figures:
%   - Mean imagesc correlation matrices (within-trial, trial-MA)
%   - Spread/box plots for pairwise r's (within-trial)
%   - Spread/box plots for CV AUC (pupil, gamma, both)
%   - Spread/box plots for GLM slopes from the combined model (beta_pupil, beta_gamma)
% Saves figures to <root>/Figures/dataset as PNG + SVG and a dataset .mat.

%% Options
if nargin<2 || isempty(signalType), signalType = 'area'; end
signalType = lower(signalType);
assert(ismember(signalType,{'area','pos'}),'signalType must be ''area'' or ''pos''');

assert(iscell(sessionPaths) && ~isempty(sessionPaths), 'sessionPaths must be a non-empty cell array of session folders.');
root   = fileparts(sessionPaths{1});
outdir = fullfile(root,'Figures',['dataset_' signalType]);
if ~exist(outdir,'dir'), mkdir(outdir); end

%%
nS = numel(sessionPaths);
R_within = nan(3,3,nS);
R_MA = nan(3,3,nS);
R_pre = nan(3,3,nS);
r_pupil_hit = nan(nS,1);
r_gamma_hit = nan(nS,1);
auc_p = nan(nS,1); auc_g = nan(nS,1); auc_b = nan(nS,1);
beta_p = nan(nS,1); beta_g = nan(nS,1);
pc_within_pupil_given_gamma = nan(nS,1);
pc_within_gamma_given_pupil = nan(nS,1);
pc_pre_pupil_given_gamma = nan(nS,1);
pc_pre_gamma_given_pupil = nan(nS,1);

kept = false(nS,1);

for s = 1:nS
    try
        statsFile = fullfile(sessionPaths{s},'Figures',['session_stats_' signalType '.mat']);
        tmp = load(statsFile,'out'); out = tmp.out;

        if isfield(out,'R'),    R_within(:,:,s) = out.R;    end
        if isfield(out,'R_sm'), R_MA(:,:,s)     = out.R_sm; end
        if isfield(out,'R_pre'),R_pre(:,:,s)    = out.R_pre;end

        r_pupil_hit(s) = safe_pick(out,'R',1,3);
        r_gamma_hit(s) = safe_pick(out,'R',2,3);

        if isfield(out,'cvAUC')
            if isfield(out.cvAUC,'pupil'), auc_p(s) = mean(out.cvAUC.pupil,'omitnan'); end
            if isfield(out.cvAUC,'gamma'), auc_g(s) = mean(out.cvAUC.gamma,'omitnan'); end
            if isfield(out.cvAUC,'both'),  auc_b(s) = mean(out.cvAUC.both, 'omitnan'); end
        end

        if isfield(out,'beta') && isfield(out.beta,'both') && numel(out.beta.both)>=3
            beta_p(s) = out.beta.both(2);
            beta_g(s) = out.beta.both(3);
        end

        if isfield(out,'partial') && numel(out.partial.within)==2
            pc_within_pupil_given_gamma(s) = out.partial.within(1);
            pc_within_gamma_given_pupil(s) = out.partial.within(2);
        end
        if isfield(out,'partial') && numel(out.partial.pre)==2
            pc_pre_pupil_given_gamma(s) = out.partial.pre(1);
            pc_pre_gamma_given_pupil(s) = out.partial.pre(2);
        end

        kept(s) = true;
    catch ME
        warning('Session %d skipped: %s', s, ME.message);
    end
end

% Keep only sessions that loaded
fn = @(v) v(kept,:);
R_within = R_within(:,:,kept); R_MA = R_MA(:,:,kept); R_pre = R_pre(:,:,kept);
r_pupil_hit = r_pupil_hit(kept); r_gamma_hit = r_gamma_hit(kept);
auc_p = auc_p(kept); auc_g = auc_g(kept); auc_b = auc_b(kept);
beta_p = beta_p(kept); beta_g = beta_g(kept);
pc_within_pupil_given_gamma = pc_within_pupil_given_gamma(kept);
pc_within_gamma_given_pupil = pc_within_gamma_given_pupil(kept);
pc_pre_pupil_given_gamma = pc_pre_pupil_given_gamma(kept);
pc_pre_gamma_given_pupil = pc_pre_gamma_given_pupil(kept);
sessionPaths = sessionPaths(kept); nS = sum(kept);

%%
% [p_r,~,stats_r]   = signrank(r_gamma_hit, r_pupil_hit);      % ? vs pupil
% [p_auc,~,stats_a] = signrank(auc_b, auc_g);                  % both vs ?
% [p_b,~,stats_b]   = signrank(beta_g);                        % ?_? > 0 ?
% 
% fprintf('Across sessions: r ?–hit > pupil–hit? p=%.3g (signed-rank z=%.2f)\n', p_r, stats_r.zval);
% fprintf('AUC both > ?? p=%.3g (z=%.2f)\n', p_auc, stats_a.zval);
% fprintf('?_gamma > 0 ? p=%.3g (z=%.2f)\n', p_b, stats_b.zval);

%% Symmetric Fisher means
Rmean_within = fisher_mean(R_within);
Rmean_MA = fisher_mean(R_MA);
Rmean_pre = fisher_mean(R_pre);

% imagesc: within-trial
f1 = figure('Color','w','Position',[180 100 460 420]);
try, imagesc(Rmean_within, [-1 1]); colormap('redblue'); catch, imagesc(Rmean_within, [-1 1]); colormap(parula); end
axis square; colorbar; caxis([-1 1]);
set(gca,'XTick',1:3,'XTickLabel',{'pupil','gamma','hit'},'YTick',1:3,'YTickLabel',{'pupil','gamma','hit'});
title('Mean correlation across sessions (within-trial means)'); box on
overlay_vals(Rmean_within);
saveas(f1, fullfile(outdir,'dataset_corr_within.png'));  try, saveas(f1, fullfile(outdir,'dataset_corr_within.svg')); end

% imagesc: trial-MA
if any(isfinite(Rmean_MA),'all')
    f2 = figure('Color','w','Position',[180 100 460 420]);
    try, imagesc(Rmean_MA, [-1 1]); colormap('redblue'); catch, imagesc(Rmean_MA, [-1 1]); colormap(parula); end
    axis square; colorbar; caxis([-1 1]);
    set(gca,'XTick',1:3,'XTickLabel',{'pupil','gamma','hit'},'YTick',1:3,'YTickLabel',{'pupil','gamma','hit'});
    title('Mean correlation across sessions (trial-MA features)'); box on
    overlay_vals(Rmean_MA);
    saveas(f2, fullfile(outdir,'dataset_corr_MA.png'));  try, saveas(f2, fullfile(outdir,'dataset_corr_MA.svg')); end
end

% imagesc: pre-trial baseline
if any(isfinite(Rmean_pre),'all')
    f2b = figure('Color','w','Position',[180 100 460 420]);
    try, imagesc(Rmean_pre, [-1 1]); colormap('redblue'); catch, imagesc(Rmean_pre, [-1 1]); colormap(parula); end
    axis square; colorbar; caxis([-1 1]);
    set(gca,'XTick',1:3,'XTickLabel',{'pupil(pre)','gamma(pre)','hit'},'YTick',1:3,'YTickLabel',{'pupil(pre)','gamma(pre)','hit'});
    title(['Mean correlation across sessions (pre-trial) | ' upper(signalType)]); box on
    overlay_vals(Rmean_pre);
    saveas(f2b, fullfile(outdir,'dataset_corr_pre.png'));  try, saveas(f2b, fullfile(outdir,'dataset_corr_pre.svg')); end
end

%% Box/Spread of pairwise r’s (within-trial)
hasMS = exist('MakeSpreadAndBoxPlot3_SB','file')==2;
f3 = figure('Color','w','Position',[160 100 260 320]); hold on
Acell = {r_pupil_hit(isfinite(r_pupil_hit)), r_gamma_hit(isfinite(r_gamma_hit))};
if hasMS && nnz(cellfun(@(v)~isempty(v), Acell))>=1
    MakeSpreadAndBoxPlot3_SB(Acell, {}, [], {'pupil-hit','gamma-hit'},'ShowSigstar','none','paired',1,'showpoints',1);
else
    Y = [r_pupil_hit, r_gamma_hit]; for i=1:2, boxchart(i*ones(nS,1),Y(:,i)); hold on; end
    set(gca,'XTick',1:2,'XTickLabel',{'pupil-hit','gamma-hit'});
end
hold on ;yline(0,'--r', 'LineWidth', 1)
ylabel('r'); ylim([-1 1]); title(['Pairwise r | ' upper(signalType)]); box off; grid on
saveas(f3, fullfile(outdir,'dataset_r_within.png'));  try, saveas(f3, fullfile(outdir,'dataset_r_within.svg')); end

% AUC
f4 = figure('Color','w','Position',[160 100 260 320]); hold on
Acell = {auc_p(isfinite(auc_p)), auc_g(isfinite(auc_g)), auc_b(isfinite(auc_b))};
if hasMS && nnz(cellfun(@(v)~isempty(v), Acell))>=1
    MakeSpreadAndBoxPlot3_SB(Acell, {}, [], {'pupil','gamma','both'}, 'ShowSigstar','none','paired',0,'showpoints',1);
else
    Y = [auc_p, auc_g, auc_b]; for i=1:3, boxchart(i*ones(nS,1),Y(:,i)); hold on; end
    set(gca,'XTick',1:3,'XTickLabel',{'pupil','gamma','both'});
end
ylabel('CV AUC'); ylim([0.5 1]); title(['Logistic regression (5-fold) | ' upper(signalType)]); box off; grid on
saveas(f4, fullfile(outdir,'dataset_auc.png'));  try, saveas(f4, fullfile(outdir,'dataset_auc.svg')); end

% Betas
f5 = figure('Color','w','Position',[160 100 260 320]); hold on
Acell = {beta_p(isfinite(beta_p)), beta_g(isfinite(beta_g))};
if hasMS && nnz(cellfun(@(v)~isempty(v), Acell))>=1
    MakeSpreadAndBoxPlot3_SB(Acell, {}, [], {'\beta_{pupil}','\beta_{\gamma}'}, 'ShowSigstar','none','paired',0,'showpoints',1);
else
    Y = [beta_p, beta_g]; for i=1:2, boxchart(i*ones(nS,1),Y(:,i)); hold on; end
    set(gca,'XTick',1:2,'XTickLabel',{'beta pupil','beta gamma'});
end
ylabel('\beta (log-odds per +1 SD)'); title(['Combined model coefficients | ' upper(signalType)]); box off; grid on
disp('? is log-odds change for a +1 SD increase in the predictor (others held fixed); exp(?)=odds ratio.');
saveas(f5, fullfile(outdir,'dataset_glm_betas.png'));  try, saveas(f5, fullfile(outdir,'dataset_glm_betas.svg')); end

% ---- NEW: Odds ratio (exp(beta)) histograms with small bins ----
bp = beta_p(isfinite(beta_p));  bg = beta_g(isfinite(beta_g));
ORp = exp(bp); ORg = exp(bg);

f5b = figure('Color','w','Position',[160 100 560 320]); hold on
edges_or = 0.5:0.1:8;   % finer, linear bins
histogram(ORp, edges_or, 'Normalization','pdf', 'FaceColor',[0.6 0.6 0.6], 'FaceAlpha',0.45, 'DisplayName','OR_{pupil}');
histogram(ORg, edges_or, 'Normalization','pdf', 'FaceColor',[0.3 0.3 0.3], 'FaceAlpha',0.45, 'DisplayName','OR_{\gamma}');
xline(1,'--r','DisplayName','OR = 1'); legend('Location','best');
grid on; box off
xlabel('Odds ratio = exp(\beta)'); ylabel('density');
title(['Odds ratios (per +1 SD) | ' upper(signalType)]);

med_bp = median(bp,'omitnan'); ci_bp = prctile(bp,[2.5 97.5]);
med_bg = median(bg,'omitnan'); ci_bg = prctile(bg,[2.5 97.5]);
fprintf('GLM ? (pupil): median=%.3f 95%%CI[%.3f, %.3f] ? OR median=%.2f\n', med_bp, ci_bp(1), ci_bp(2), exp(med_bp));
fprintf('GLM ? (gamma): median=%.3f 95%%CI[%.3f, %.3f] ? OR median=%.2f\n', med_bg, ci_bg(1), ci_bg(2), exp(med_bg));

saveas(f5b, fullfile(outdir,'dataset_glm_odds_ratio.png'));
try, saveas(f5b, fullfile(outdir,'dataset_glm_odds_ratio.svg')); end

% Partial correlations (within-trial)
f6 = figure('Color','w','Position',[160 100 260 320]); hold on
Acell = {pc_within_pupil_given_gamma(isfinite(pc_within_pupil_given_gamma)), ...
         pc_within_gamma_given_pupil(isfinite(pc_within_gamma_given_pupil))};
if hasMS && nnz(cellfun(@(v)~isempty(v), Acell))>=1
    MakeSpreadAndBoxPlot3_SB(Acell, {}, [], {'r(hit, pupil | gamma)','r(hit, gamma | pupil)'}, 'ShowSigstar','none','paired',1,'showpoints',1);
else
    Y = [pc_within_pupil_given_gamma, pc_within_gamma_given_pupil]; for i=1:2, boxchart(i*ones(nS,1),Y(:,i)); hold on; end
    set(gca,'XTick',1:2,'XTickLabel',{'r(hit,pupil|gamma)','r(hit,gamma|pupil)'});
end
hold on ;yline(0,'--r', 'LineWidth', 1)
ylabel('partial r'); ylim([-1 1]); title(['Partial (within-trial) | ' upper(signalType)]); box off; grid on
saveas(f6, fullfile(outdir,'dataset_partial_within.png'));  try, saveas(f6, fullfile(outdir,'dataset_partial_within.svg')); end

% Partial correlations (pre-trial)
f7 = figure('Color','w','Position',[160 100 260 320]); hold on
Acell = {pc_pre_pupil_given_gamma(isfinite(pc_pre_pupil_given_gamma)), ...
         pc_pre_gamma_given_pupil(isfinite(pc_pre_gamma_given_pupil))};
if hasMS && nnz(cellfun(@(v)~isempty(v), Acell))>=1
    MakeSpreadAndBoxPlot3_SB(Acell, {}, [], {'r(hit,pupil pre|gamma pre)','r(hit,gamma pre|pupil pre)'}, 'ShowSigstar','none','paired',1,'showpoints',1);
else
    Y = [pc_pre_pupil_given_gamma, pc_pre_gamma_given_pupil]; for i=1:2, boxchart(i*ones(nS,1),Y(:,i)); hold on; end
    set(gca,'XTick',1:2,'XTickLabel',{'r(hit,pupil{pre}|gamma{pre})','r(hit,gamma{pre}|pupil{pre})'});
end
hold on ;yline(0,'--r', 'LineWidth', 1)
ylabel('partial r'); ylim([-1 1]); title(['Partial (pre-trial) | ' upper(signalType)]); box off; grid on
saveas(f7, fullfile(outdir,'dataset_partial_pre.png'));  try, saveas(f7, fullfile(outdir,'dataset_partial_pre.svg')); end

%% Dataset-level distributions (pairwise r, partial r, ?AUC) -- smaller bins
deltaAUC = auc_b - auc_g;

fD = figure('Color','w','Position',[160 100 980 640]);

% finer, consistent bins
edges_r   = -1:0.05:1;       % correlations
edges_pr  = -1:0.05:1;       % partial correlations
edges_dA  = -0.20:0.01:0.20; % ?AUC

subplot(2,2,1); hold on
histogram(r_pupil_hit(isfinite(r_pupil_hit)), edges_r, 'Normalization','pdf', 'FaceAlpha',0.30, 'FaceColor',[0 0 0],   'DisplayName','pupil?hit');
histogram(r_gamma_hit(isfinite(r_gamma_hit)), edges_r, 'Normalization','pdf', 'FaceAlpha',0.30, 'FaceColor',[0.4 0.4 1], 'DisplayName','\gamma?hit');
xline(median(r_pupil_hit,'omitnan'),'--k','DisplayName','med(pupil)');
xline(median(r_gamma_hit,'omitnan'),'--','Color',[0.3 0.3 1],'DisplayName','med(\gamma)');
ylabel('density'); xlabel('r'); title(['Pairwise r | ' upper(signalType)]);
grid on; box off; legend('Location','best')

subplot(2,2,2); hold on
histogram(pc_within_gamma_given_pupil(isfinite(pc_within_gamma_given_pupil)), edges_pr, 'Normalization','pdf', 'FaceAlpha',0.30, 'FaceColor',[0.4 0.4 1], 'DisplayName','r(hit,\gamma | pupil)');
histogram(pc_within_pupil_given_gamma(isfinite(pc_within_pupil_given_gamma)), edges_pr, 'Normalization','pdf', 'FaceAlpha',0.30, 'FaceColor',[0 0 0],   'DisplayName','r(hit,pupil | \gamma)');
xline(0,'--r','DisplayName','0');
xline(median(pc_within_gamma_given_pupil,'omitnan'),'--','Color',[0.3 0.3 1],'DisplayName','med(\gamma|pupil)');
xline(median(pc_within_pupil_given_gamma,'omitnan'),'--k','DisplayName','med(pupil|\gamma)');
ylabel('density'); xlabel('partial r'); title('Partial (within-trial)');
grid on; box off; legend('Location','best')

subplot(2,2,3); hold on
histogram(pc_pre_gamma_given_pupil(isfinite(pc_pre_gamma_given_pupil)), edges_pr, 'Normalization','pdf', 'FaceAlpha',0.30, 'FaceColor',[0.4 0.4 1], 'DisplayName','r(hit,\gamma_{pre}|pupil_{pre})');
histogram(pc_pre_pupil_given_gamma(isfinite(pc_pre_pupil_given_gamma)), edges_pr, 'Normalization','pdf', 'FaceAlpha',0.30, 'FaceColor',[0 0 0],   'DisplayName','r(hit,pupil_{pre}|\gamma_{pre})');
xline(0,'--r','DisplayName','0');
xline(median(pc_pre_gamma_given_pupil,'omitnan'),'--','Color',[0.3 0.3 1],'DisplayName','med(\gamma_{pre}|pupil_{pre})');
xline(median(pc_pre_pupil_given_gamma,'omitnan'),'--k','DisplayName','med(pupil_{pre}|\gamma_{pre})');
ylabel('density'); xlabel('partial r'); title('Partial (pre-trial)');
grid on; box off; legend('Location','best')

subplot(2,2,4); hold on
histogram(deltaAUC(isfinite(deltaAUC)), edges_dA, 'Normalization','pdf', 'FaceAlpha',0.30, 'FaceColor',[0.2 0.2 0.2], 'DisplayName','\DeltaAUC');
xline(0,'--r','DisplayName','0');
xline(median(deltaAUC,'omitnan'),'--k','DisplayName','median');
ylabel('density'); xlabel('\DeltaAUC = AUC(both) - AUC(\gamma)'); title('\DeltaAUC across sessions');
grid on; box off; legend('Location','best')

saveas(fD, fullfile(outdir, 'dataset_stats_distributions.png'));
try, saveas(fD, fullfile(outdir, 'dataset_stats_distributions.svg')); end


%% Save dataset table/mats
OR_p = exp(beta_p); OR_g = exp(beta_g);

T = table(sessionPaths(:), ...
          r_pupil_hit, r_gamma_hit, ...
          auc_p, auc_g, auc_b, ...
          beta_p, beta_g, OR_p, OR_g, ...
          pc_within_pupil_given_gamma, pc_within_gamma_given_pupil, ...
          pc_pre_pupil_given_gamma,    pc_pre_gamma_given_pupil, ...
          'VariableNames', {'session', ...
                            'r_pupil_hit','r_gamma_hit', ...
                            'auc_pupil','auc_gamma','auc_both', ...
                            'beta_pupil','beta_gamma','OR_pupil','OR_gamma', ...
                            'pc_within_pupil_given_gamma','pc_within_gamma_given_pupil', ...
                            'pc_pre_pupil_given_gamma','pc_pre_gamma_given_pupil'});
                        
save(fullfile(outdir,'dataset_table.mat'),'T','R_within','R_MA','R_pre','Rmean_within','Rmean_MA','Rmean_pre');

S = struct('table',T,'Rmean_within',Rmean_within,'Rmean_MA',Rmean_MA,'Rmean_pre',Rmean_pre,'outdir',outdir);

end

%% --------- helpers ----------
function val = safe_pick(out, fieldname, i, j)
val = NaN;
if isfield(out, fieldname)
    M = out.(fieldname);
    if all(size(M) >= [i j]) && isfinite(M(i,j)), val = M(i,j); end
end
end

function Rm = fisher_mean(Rstack)
% Rstack: 3x3xN
Rm = nan(3,3);
if isempty(Rstack), return; end
for i = 1:3
    for j = 1:3
        if i==j
            Rm(i,j) = 1;
        else
            rij = squeeze(Rstack(i,j,:));
            rij = max(min(rij,0.999999),-0.999999);  % avoid Inf in atanh
            z   = atanh(rij);
            Rm(i,j) = tanh(nanmean(z));    % mean in z-space
        end
    end
end
% enforce symmetry just in case of unequal NaNs
Rm = (Rm + Rm.')/2;
end

function overlay_vals(M)
for i=1:size(M,1)
    for j=1:size(M,2)
        val = M(i,j);
        txt = sprintf('%.2f',val);
        if abs(val)>0.6, tc='w'; else, tc='k'; end
        text(j,i,txt,'HorizontalAlignment','center','FontWeight','bold','Color',tc);
    end
end
end
