function behaviour_analysis_dataset(sessionPaths, SW, bodypart, marker)

% function behaviour_analysis_dataset(sessionPaths, SW, AW, bodypart, marker)
% behaviour_analysis_dataset  Aggregate and visualise marker metrics across sessions
%
%   behaviour_analysis_dataset(rootDir, bodypart, marker)
%
% Inputs
%   rootDir  – directory that contains all session folders (one level down)
%   bodypart – e.g. 'Pupil'   (exactly the same name used in behaviour_analysis)
%   marker   – e.g. 'pupil_area_004'
%
% Example
%   across_session_analysis('Z:\Arsenii\React_Active\training\Tvorozhok', 'Pupil_EyeCam', 'pupil_area_007')

%% Parameters of your choice to label figures:
AW = [3.7 4.1];

%% 1. locate every metrics file ------------------------------------------------
files = struct('folder',{},'name',{});

for p = 1:numel(sessionPaths)
    metricsDir = fullfile(sessionPaths{p}, 'analysis', 'behaviour', bodypart);
    hit        = dir(fullfile(metricsDir, '*all_metrics.mat'));
    if isempty(hit)
        warning('No *all_metrics.mat file found in %s – skipping', metricsDir);
        continue
    end
    % keep only the first hit (normally exactly one)
    files(end+1).folder = hit(1).folder; %#ok<AGROW>
    files(end  ).name   = hit(1).name;
end

if isempty(files)
    error('No metrics files found in the provided session paths.');
end

files = files(:);

%% 2. load each table into a big table ---------------------------------------
bigT = table();
sessions = strings(numel(files),1);

for k = 1:numel(files)
    L = load(fullfile(files(k).folder, files(k).name), 'T');
    if ~isfield(L,'T') || isempty(L.T),  continue, end
    T = L.T;

    % ---- FIX: correct way to check for a marker column in a table
    if ~any(strcmp(T.Properties.VariableNames, marker))
        warning('Marker "%s" not found in %s – skipping', marker, fullfile(files(k).folder, files(k).name));
        continue
    end

    % --- flatten metrics (stay close to original) -----------------------
    S = T.(marker);                    % may be a struct or a cell containing a struct
    if iscell(S), S = S{1}; end        % unwrap cell to struct if needed
    if istable(S), S = table2struct(S); end
    if ~isstruct(S) || numel(S) ~= 1
        warning('Marker "%s" in %s is not a 1x1 struct – skipping', marker, files(k).name);
        continue
    end

    metricsWanted = {'AshmanD','AUROC','CohensD','meanDiff','meanRef','meanTar','stdRef','stdTar'};
    % make sure all requested fields exist; fill missing with NaN (robust)
    miss = setdiff(metricsWanted, fieldnames(S));
    for m = 1:numel(miss), S.(miss{m}) = NaN; end
    % drop everything else to avoid vector fields exploding struct2table
    S = rmfield(S, setdiff(fieldnames(S), metricsWanted));

    metricsTbl = struct2table(S);
    T = [ removevars(T, marker) , metricsTbl ];  %#ok<AGROW>

    % --- add meta-data (kept identical in spirit) -----------------------
    [~,sessName] = fileparts(sessionPaths{k});
    T.Session    = repmat(string(sessName), height(T), 1);
    T.BodyPart   = repmat(string(bodypart), height(T), 1);
    T.Marker     = repmat(string(marker),   height(T), 1);

    % --- accumulate ------------------------------------------------------
    bigT = [bigT ;
        T(:,{'Session','AshmanD','AUROC','CohensD','meanDiff','meanRef','meanTar','stdRef','stdTar'})];
    sessions(k) = sessName;
end

if isempty(bigT)
    error('Metrics loaded but empty—check column names in the per-session files.');
end

sessions  = sessions(sessions~="");      % remove skipped ones
% keep only sessions that actually appear in bigT to sync x-axis
sessOrder = unique(bigT.Session,'stable');
nSess     = numel(sessOrder);

% --- NEW: compute one row per session for plotting (means) ------------------
mNamesAll = {'AshmanD','AUROC','CohensD','meanDiff','meanRef','meanTar','stdRef','stdTar'};
sessT = groupsummary(bigT,'Session','mean',mNamesAll);
% rename mean_* back to original names
for i = 1:numel(mNamesAll)
    old = "mean_" + mNamesAll{i};
    if any(strcmp(sessT.Properties.VariableNames, old))
        sessT.(mNamesAll{i}) = sessT.(old);
        sessT.(old) = [];
    end
end
% reorder rows to the stable session order used on x-axis
[~,ord] = ismember(sessOrder, sessT.Session);
sessT = sessT(ord,:);
sessList = sessOrder;  % used by the plotting helper below

%% 3. figure 1  : metrics vs session -----------------------------------------
f{1} = figure('Color','w','Visible','on', 'Position', [0 0 1920 1080]);
ax = gobjects(6,1);
for n = 1:6
    ax(n) = subplot(3,2,n);     % classical placement
    set(ax(n),'Box','off');     % cosmetic: no surrounding box
end
x = (1:nSess)';

% Helper for trend line that ignores NaNs and annotates slope & R2
    function plot_with_trend(xv,yv,clr)
        idx = isfinite(xv) & isfinite(yv);
        plot(xv(idx), yv(idx), 'o-','LineWidth',1.4,'Color',clr,'MarkerFaceColor',clr*0+0.05); hold on
        if sum(idx) >= 2
            p = polyfit(xv(idx), yv(idx), 1);
            yHat = polyval(p, xv(idx));
            plot(xv(idx), yHat, '--','LineWidth',1.4,'Color',clr);
            SSres = sum((yv(idx) - yHat).^2);
            SStot = sum((yv(idx) - mean(yv(idx))).^2);
            Rsq   = 1 - SSres / max(SStot, eps);
            txt = sprintf('slope = %.3g /sess\nR^2 = %.2f', p(1), Rsq);
            text(0.04, 0.90, txt, 'Units','normalized','FontSize',8,'Color',clr);
        end
        xlim([0.7 nSess+0.3]); xticks(x); xticklabels(sessList); xtickangle(45); grid on
    end

% Colors
cAsh   = [0 0 0];
cAuroc = [.25 .5 .9];
cCoh   = [.3 .3 .3];
cRef   = [0.20 0.60 1.00];
cTar   = [1.00 0.40 0.40];
cDiff  = [.2 .7 .2];

% (1) AshmanD
axes(ax(1)); plot_with_trend(x, sessT.AshmanD, cAsh); title('Ashman D');

% (2) AUROC
axes(ax(2)); plot_with_trend(x, sessT.AUROC, cAuroc); title('AUROC');

% (3) Cohen''s d
axes(ax(3)); plot_with_trend(x, sessT.CohensD, cCoh); title('Cohen''s d');

% (4) meanDiff (Tar-Ref in your per-session metric)
axes(ax(4)); plot_with_trend(x, sessT.meanDiff, cDiff); title('meanDiff (Tar - Ref)');

% (5) Combined: meanRef vs meanTar on ONE plot (optional error bars if stds exist)
axes(ax(5)); hold on
hasStd = ismember('stdRef', sessT.Properties.VariableNames) && ismember('stdTar', sessT.Properties.VariableNames);

if hasStd
    er1 = errorbar(x, sessT.meanRef, sessT.stdRef, 'o-','LineWidth',1.2,'Color',cRef,'MarkerFaceColor',cRef*0+0.05);
    er2 = errorbar(x, sessT.meanTar, sessT.stdTar, 'o-','LineWidth',1.2,'Color',cTar,'MarkerFaceColor',cTar*0+0.05);
else
    er1 = plot(x, sessT.meanRef, 'o-','LineWidth',1.4,'Color',cRef,'MarkerFaceColor',cRef*0+0.05);
    er2 = plot(x, sessT.meanTar, 'o-','LineWidth',1.4,'Color',cTar,'MarkerFaceColor',cTar*0+0.05);
end
% Trends for each
plot_with_trend(x, sessT.meanRef, cRef);
plot_with_trend(x, sessT.meanTar, cTar);
title('Session means: Ref vs Tar (one axis)');
legend([er1 er2], {'Ref','Tar'}, 'Location','best'); legend boxoff

% (6) Explicit Tar-Ref per session (derived)
axes(ax(6));
tarMinusRef = sessT.meanTar - sessT.meanRef;
plot_with_trend(x, tarMinusRef, cDiff);
yline(0,'k:'); title('Tar - Ref (derived)');

sgtitle(sprintf('%s – %s | metrics across sessions  (SW=%.3g, AW=[%.1f %.1f])', ...
    bodypart, marker, SW, AW(1), AW(2)));

%% 4. figure 2  : session-mean Ref vs Tar -------------------------------------
refSess = sessT.meanRef;
tarSess = sessT.meanTar;

f{2} = figure('Color','w', 'Position', [0 0 350 500]);
MakeSpreadAndBoxPlot3_SB({refSess, tarSess}, ...
    {[0.2 0.6 1], [1 0.4 0.4]}, ...
    [1 2], {'Ref','Tar'}, ...
    'paired',1, 'optiontest','ttest', ...
    'newfig',0,'showpoints',1);
ylabel(sprintf('Session mean %s (anticipatory window)', marker))
title('Across-sess (paired t-test)')

% paired t-test on session means
[~,p,ci,stats] = ttest(tarSess, refSess);
fprintf('\nPaired t-test for %s Tar vs Ref across %d sessions:  t = %.3f  p = %.4g\n', ...
    marker, nSess, stats.tstat, p);
fprintf('95%% CI of difference: [%.3f  %.3f]\n', ci(1), ci(2));

% Annotate stats on the figure too
ax2 = gca;
txt = sprintf('t(%d)=%.2f  p=%.3g\nCI: [%.3f, %.3f]', stats.df, stats.tstat, p, ci(1), ci(2));
text(ax2, 0.05, 0.95, txt, 'Units','normalized','VerticalAlignment','top','FontSize',9);

%% save figures
outDir = fullfile(fileparts(sessionPaths{1}), 'DS_figures', 'behaviour');
if ~exist(outDir,'dir'), mkdir(outDir), end
saveas(f{1},fullfile(outDir,[bodypart '_' marker '_sw_' num2str(SW) '_aw_' num2str(round(AW(1))) '_' num2str(round(AW(2))) '_trend.png']));
saveas(f{2},fullfile(outDir,[bodypart '_' marker '_sw_' num2str(SW) '_aw_' num2str(round(AW(1))) '_' num2str(round(AW(2))) '_ttest.png']));
end
