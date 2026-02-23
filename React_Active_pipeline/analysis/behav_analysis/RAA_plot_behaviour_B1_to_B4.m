function RAA_plot_behaviour_B1_to_B4(BigT, SessT, markerName, outRoot)
% RAA_plot_behaviour_B1_to_B4
% Uses ONLY one selected window for plotting (default: 'stimoff_to_arrival').
% B3 AUROC threshold = 0.6.

if nargin < 4 || isempty(outRoot)
    outRoot = fullfile(pwd, 'DS_figures', 'behaviour');
end
if ~exist(outRoot,'dir'), mkdir(outRoot); end

% -------- selection knobs ----------
winSel  = 'stimoff_to_arrival';
winMode = 'exact';   % 'exact' or 'contains'
thrAUC  = 0.55;
% ---------------------------------

mk = char(string(markerName));

% Ensure required columns exist and are comparable (cellstr)
BigT = ra_ensure_plot_columns(BigT);
SessT = ra_ensure_sess_columns(SessT);

% markerName can be either:
%   - an actual marker (e.g. pupil_area_007)
%   - a bodypart label (e.g. Pupil-EyeCam)
% In all cases: plot ONE marker at a time, never average across markers, never use __MEAN__.

mask_marker  = strcmpi_safe(BigT.Marker, mk);
mask_bodypart = strcmpi_safe(BigT.BodyPart, mk);

if any(mask_marker)
    % user passed an explicit marker
    markers_to_plot = {mk};
elseif any(mask_bodypart)
    % user passed a bodypart: plot all markers within it (excluding __MEAN__)
    Tb = BigT(mask_bodypart, :);
    mlist = unique(ra_to_cellstr(Tb.Marker), 'stable');
    mlist = mlist(~strcmpi(mlist,'__MEAN__'));  % exclude __MEAN__ always
    markers_to_plot = mlist(:)';
else
    warning('No rows found for marker/bodypart %s', mk);
    return
end

if isempty(markers_to_plot)
    warning('No valid markers to plot for %s (after excluding __MEAN__).', mk);
    return
end

for im = 1:numel(markers_to_plot)
    mk_i = markers_to_plot{im};

    Ti = BigT(strcmpi_safe(BigT.Marker, mk_i), :);
    if isempty(Ti) || height(Ti)==0
        continue
    end

    animals = unique(Ti.Animal,'stable');
    for ia = 1:numel(animals)
        plot_set_for_one_animal(Ti, SessT, mk_i, animals{ia}, outRoot, winSel, winMode, thrAUC);
    end
    plot_set_for_one_animal(Ti, SessT, mk_i, 'pooled', outRoot, winSel, winMode, thrAUC);
end

end

% ======================================================================
function plot_set_for_one_animal(T, SessT, mk, animal, outRoot, winSel, winMode, thrAUC)

FS = 11; LW = 1.8; MS = 5; LWthin = 1.1;
cAuroc = [.25 .5 .9];
cDiff  = [.2 .7 .2];
cClassic = [0 0 0];
cOperant = [.25 .25 .25];

subsets = {'regular','nosound','nomotor'};

if strcmpi(animal,'pooled')
    Ta = T;
    Sa = SessT;
    tag = 'POOLED';
else
    Ta = T(strcmpi_safe(T.Animal, animal), :);
    Sa = SessT(strcmpi_safe(SessT.Animal, animal), :);
    tag = animal;
end
if isempty(Ta) || height(Ta)==0, return; end
if isempty(Sa) || height(Sa)==0, return; end
mk_file = regexprep(mk,'[^a-zA-Z0-9_]+','_');

% ---- select window for plotting ----
Ta = filter_by_window(Ta, winSel, winMode);
if isempty(Ta) || height(Ta)==0
    warning('No rows left after window selection: %s (%s) for %s', winSel, winMode, tag);
    return
end

[Sa, sessNames, x_global, blockChangeX, blockLabels, motorSwitchX] = make_session_axis(Sa, tag);
nSess = numel(sessNames);

% Robust block-aligned x-axis:
blockIdxMap = make_block_index_map(Sa);

% Collapse duplicates per session+subset+window (mean over duplicates)
TaC = collapse_by_session_subset_window(Ta);

% ---------------- B1a (global) ----------------
ylA = ylims_across_subsets(TaC, sessNames, subsets, 'AUROC', [0.35 0.75]);
ylD = ylims_across_subsets(TaC, sessNames, subsets, 'meanDiff', []);

f1 = figure('Color','w','Position',[50 50 1300 900], 'visible', 'off');

for r = 1:3
    subset = subsets{r};

    ax = subplot(3,2,(r-1)*2+1);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);

    if strcmpi(tag,'POOLED')
        uA = unique(TaC.Animal,'stable');
        for ii = 1:numel(uA)
            Ts = TaC(strcmpi_safe(TaC.Subset,subset) & strcmpi_safe(TaC.Animal,uA{ii}), :);
            y = series_by_session(Ts, sessNames, 'AUROC');
            plot(ax, x_global, y, 'o-','LineWidth',LW,'MarkerSize',MS,'DisplayName',uA{ii});
            add_trend(ax, x_global, y, LWthin);
        end
        legend(ax,'Location','best'); legend(ax,'boxoff');
    else
        Ts = TaC(strcmpi_safe(TaC.Subset,subset), :);
        y = series_by_session(Ts, sessNames, 'AUROC');
        plot(ax, x_global, y, 'o-','Color',cAuroc,'LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',cAuroc);
        add_trend(ax, x_global, y, LWthin);
    end

    yline(ax,0.5,'k:','LineWidth',LWthin);
    add_block_lines(ax, blockChangeX, blockLabels);
    xlim(ax,[0.7 nSess+0.3]);
    ylim(ax, ylA);
    xticks(ax, nice_xticks(nSess));
    ylabel(ax,'AUROC');
    title(ax, sprintf('%s | AUROC', subset), 'Interpreter','none','FontWeight','bold');
    grid(ax,'on');

    ax = subplot(3,2,(r-1)*2+2);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);

    if strcmpi(tag,'POOLED')
        uA = unique(TaC.Animal,'stable');
        for ii = 1:numel(uA)
            Ts = TaC(strcmpi_safe(TaC.Subset,subset) & strcmpi_safe(TaC.Animal,uA{ii}), :);
            y = series_by_session(Ts, sessNames, 'meanDiff');
            plot(ax, x_global, y, 'o-','LineWidth',LW,'MarkerSize',MS,'DisplayName',uA{ii});
            add_trend(ax, x_global, y, LWthin);
        end
        legend(ax,'Location','best'); legend(ax,'boxoff');
    else
        Ts = TaC(strcmpi_safe(TaC.Subset,subset), :);
        y = series_by_session(Ts, sessNames, 'meanDiff');
        plot(ax, x_global, y, 'o-','Color',cDiff,'LineWidth',LW,'MarkerSize',MS,'MarkerFaceColor',cDiff);
        add_trend(ax, x_global, y, LWthin);
    end

    yline(ax,0,'k:','LineWidth',LWthin);
    add_block_lines(ax, blockChangeX, blockLabels);
    xlim(ax,[0.7 nSess+0.3]);
    if ~isempty(ylD), ylim(ax, ylD); end
    xticks(ax, nice_xticks(nSess));
    ylabel(ax,'meanDiff');
    title(ax, sprintf('%s | meanDiff (Tar-Ref)', subset), 'Interpreter','none','FontWeight','bold');
    grid(ax,'on');
end

xlabel(subplot(3,2,5),'Global session #');
xlabel(subplot(3,2,6),'Global session #');
sgtitle(sprintf('B1a %s | %s | window=%s', tag, mk, winSel), 'Interpreter','none','FontWeight','bold');

saveas(f1, fullfile(outRoot, sprintf('%s_%s_B1a_global_%s.png', tag, mk_file, winSel)));
close(f1);


% ---------------- B1b (block-aligned, regular) ----------------
Tr = TaC(strcmpi_safe(TaC.Subset,'regular'), :);
if ~isempty(Tr) && height(Tr)>0
    f2 = figure('Color','w','Position',[50 50 600 700], 'visible', 'off');
    blocks = unique(Tr.SoundPair,'stable');
    C = gradient_multicolor(numel(blocks));

    ax = subplot(2,1,1);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
    for b = 1:numel(blocks)
        Tb = Tr(strcmpi_safe(Tr.SoundPair, blocks{b}), :);
        [x,y] = xy_sorted(Tb, blockIdxMap, 'AUROC');
        plot(ax, x, y, 'o-','Color',C(b,:),'LineWidth',LW,'MarkerSize',MS,'DisplayName',blocks{b});
    end
    yline(ax,0.5,'k:','LineWidth',LWthin);
    ylabel(ax,'AUROC');
    title(ax, sprintf('B1b %s | AUROC aligned by sound pair (regular)', tag), 'Interpreter','none','FontWeight','bold');
    grid(ax,'on');
    legend(ax,'Location','best'); legend(ax,'boxoff');

    ax = subplot(2,1,2);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
    for b = 1:numel(blocks)
        Tb = Tr(strcmpi_safe(Tr.SoundPair, blocks{b}), :);
        [x,y] = xy_sorted(Tb, blockIdxMap, 'meanDiff');
        plot(ax, x, y, 'o-','Color',C(b,:),'LineWidth',LW,'MarkerSize',MS,'DisplayName',blocks{b});
    end
    yline(ax,0,'k:','LineWidth',LWthin);
    xlabel(ax,'Session # within sound pair'); ylabel(ax,'meanDiff');
    grid(ax,'on');

    sgtitle(sprintf('B1b %s | %s | window=%s', tag, mk, winSel), 'Interpreter','none','FontWeight','bold');
    saveas(f2, fullfile(outRoot, sprintf('%s_%s_B1b_block_aligned_%s.png', tag, mk, winSel)));
    close(f2);
end

% ---------------- B2 (first vs expert + names) ----------------
if ~isempty(Tr) && height(Tr)>0
    blocks = unique(Tr.SoundPair,'stable');
    C = gradient_multicolor(numel(blocks));
    firstA = nan(numel(blocks),1); expertA = nan(numel(blocks),1);
    firstD = nan(numel(blocks),1); expertD = nan(numel(blocks),1);
    firstName = cell(numel(blocks),1);
    expertNames = cell(numel(blocks),1);

    for b = 1:numel(blocks)
        Tb = Tr(strcmpi_safe(Tr.SoundPair, blocks{b}), :);
        [firstA(b), expertA(b), firstName{b}, expertNames{b}] = first_vs_expert_with_names(Tb, blockIdxMap, 'AUROC');
        [firstD(b), expertD(b)] = first_vs_expert(Tb, blockIdxMap, 'meanDiff');
    end

    fprintf('\nB2 definitions (%s | %s | window=%s):\n', tag, mk, winSel);
    for b = 1:numel(blocks)
        fprintf('  Block %s: FIRST=%s ; EXPERT=mean(%s)\n', blocks{b}, firstName{b}, strjoin(expertNames{b}, ', '));
    end

    f3 = figure('Color','w','Position',[50 50 650 650], 'visible', 'off');

    ax = subplot(2,2,1);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
    for b = 1:numel(blocks)
        plot(ax, [1 2], [firstA(b) expertA(b)], 'o-','Color',C(b,:),'LineWidth',LW,'MarkerSize',MS);
    end
    xlim(ax,[0.7 2.3]); xticks(ax,[1 2]); xticklabels(ax,{'First','Expert'});
    yline(ax,0.5,'k:','LineWidth',LWthin);
    title(ax,'B2 | AUROC: per sound pair','FontWeight','bold');
    ylabel(ax,'AUROC'); grid(ax,'on');

    subplot(2,2,2);
    MakeSpreadAndBoxPlot3_SB({firstA, expertA}, {[.25 .5 .9], [.25 .5 .9]}, [1 2], {'First','Expert'}, ...
        'paired',1,'newfig',0,'showpoints',1);
    set(gca,'Box','off','FontSize',FS);
    yline(0.5,'k:','LineWidth',LWthin);
    title('B2 | AUROC: paired across sound pairs','FontWeight','bold');
    ylabel('AUROC'); grid on

    ax = subplot(2,2,3);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
    for b = 1:numel(blocks)
        plot(ax, [1 2], [firstD(b) expertD(b)], 'o-','Color',C(b,:),'LineWidth',LW,'MarkerSize',MS);
    end
    xlim(ax,[0.7 2.3]); xticks(ax,[1 2]); xticklabels(ax,{'First','Expert'});
    yline(ax,0,'k:','LineWidth',LWthin);
    title(ax,'B2 | meanDiff: per sound pair','FontWeight','bold');
    ylabel(ax,'meanDiff'); grid(ax,'on');

    subplot(2,2,4);
    MakeSpreadAndBoxPlot3_SB({firstD, expertD}, {[.2 .7 .2], [.2 .7 .2]}, [1 2], {'First','Expert'}, ...
        'paired',1,'newfig',0,'showpoints',1);
    set(gca,'Box','off','FontSize',FS);
    yline(0,'k:','LineWidth',LWthin);
    title('B2 | meanDiff: paired across sound pairs','FontWeight','bold');
    ylabel('meanDiff'); grid on

    sgtitle(sprintf('B2 %s | %s | window=%s', tag, mk, winSel), 'Interpreter','none','FontWeight','bold');
    saveas(f3, fullfile(outRoot, sprintf('%s_%s_B2_first_vs_expert_%s.png', tag, mk, winSel)));
    close(f3);
end

% ---------------- B3 (speed across blocks) ----------------
if ~isempty(Tr) && height(Tr)>0
    blocks = unique(Tr.SoundPair,'stable');
    sessTo = nan(numel(blocks),1);
    slope  = nan(numel(blocks),1);

    for b = 1:numel(blocks)
        Tb = Tr(strcmpi_safe(Tr.SoundPair, blocks{b}), :);
        sessTo(b) = sessions_to_threshold(Tb, blockIdxMap, thrAUC);
        slope(b)  = early_slope(Tb, blockIdxMap, 5);
    end

    f4 = figure('Color','w','Position',[50 50 650 520], 'visible', 'off');

    ax = subplot(1,2,1);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
    bar(ax, sessTo);
    ylabel(ax, sprintf('Sessions to AUROC >= %.2f (within block)', thrAUC));
    set(ax,'XTick',1:numel(blocks),'XTickLabel',blocks); xtickangle(ax,45);
    grid(ax,'on');
    title(ax,'B3 | sessions-to-threshold','FontWeight','bold');

    ax = subplot(1,2,2);
    hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
    plot(ax, 1:numel(blocks), slope, 'o-','LineWidth',LW,'MarkerSize',MS);
    set(ax,'XTick',1:numel(blocks),'XTickLabel',blocks); xtickangle(ax,45);
    ylabel(ax,'Early slope (AUROC / session), first 5 valid pts');
    yline(ax,0,'k:','LineWidth',LWthin);
    grid(ax,'on');
    title(ax,'B3 | early learning slope','FontWeight','bold');

    sgtitle(sprintf('B3 %s | %s | window=%s', tag, mk, winSel), 'Interpreter','none','FontWeight','bold');
    saveas(f4, fullfile(outRoot, sprintf('%s_%s_B3_block_speed_%s.png', tag, mk, winSel)));
    close(f4);
end

% ---------------- B4 (classic vs operand) ----------------
if strcmpi(tag,'POOLED') || ~strcmpi(tag,'Mochi')
    return
end

Tr = TaC(strcmpi_safe(TaC.Subset,'regular'), :);
Tr = enforce_operant_sessions(Tr);

% Determine switch as the first operant session in chronological order.
Sa_m = Sa; % already sorted by make_session_axis
Sa_m = enforce_operant_sessions(Sa_m);
iSwitch = find(strcmpi_safe(Sa_m.MotorCondition,'operant'), 1, 'first');
if isempty(iSwitch), return; end

if iSwitch > 1
    classicSess = Sa_m.SessionName(1:iSwitch-1);
else
    classicSess = {};
end
operantSess = Sa_m.SessionName(iSwitch:end);

Tc = Tr(strcmpi_safe(Tr.MotorCondition,'classic') & ismember(string(Tr.SessionName), string(classicSess)), :);
To = Tr(strcmpi_safe(Tr.MotorCondition,'operant') & ismember(string(Tr.SessionName), string(operantSess)), :);
if isempty(Tc) || isempty(To) || height(Tc)==0 || height(To)==0, return; end

Tc = sortrows(Tc, {'SessionName'}); 
To = sortrows(To, {'SessionName'});
N = min(height(Tc), height(To));
TcM = Tc(end-N+1:end,:);
ToM = To(1:N,:);

f5 = figure('Color','w','Position',[50 50 900 650], 'visible', 'off');

ax = subplot(2,2,1);
hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
xC = (-N:-1)'; yC = TcM.AUROC;
xO = (1:N)';   yO = ToM.AUROC;
plot(ax, xC, yC, 'o-','Color',cClassic,'LineWidth',LW,'MarkerSize',MS,'DisplayName','classic (last N)');
plot(ax, xO, yO, 'o-','Color',cOperant,'LineWidth',LW,'MarkerSize',MS,'DisplayName','operant (first N)');
yline(ax,0.5,'k:','LineWidth',LWthin);
xline(ax,0,'k-','LineWidth',1.2);
xlabel(ax,'Sessions relative to switch'); ylabel(ax,'AUROC');
title(ax,'B4 | AUROC around switch (matched N)','FontWeight','bold');
grid(ax,'on'); legend(ax,'Location','best'); legend(ax,'boxoff');

subplot(2,2,2);
MakeSpreadAndBoxPlot3_SB({TcM.AUROC, ToM.AUROC}, {cClassic, cOperant}, [1 2], {'classic','operant'}, ...
    'paired',1,'newfig',0,'showpoints',1);
set(gca,'Box','off','FontSize',FS);
yline(0.5,'k:','LineWidth',LWthin);
title('B4 | AUROC matched (paired)','FontWeight','bold');
ylabel('AUROC'); grid on

ax = subplot(2,2,3);
hold(ax,'on'); set(ax,'Box','off','FontSize',FS);
yC = TcM.meanDiff; yO = ToM.meanDiff;
plot(ax, xC, yC, 'o-','Color',cClassic,'LineWidth',LW,'MarkerSize',MS);
plot(ax, xO, yO, 'o-','Color',cOperant,'LineWidth',LW,'MarkerSize',MS);
yline(ax,0,'k:','LineWidth',LWthin);
xline(ax,0,'k-','LineWidth',1.2);
xlabel(ax,'Sessions relative to switch'); ylabel(ax,'meanDiff');
title(ax,'B4 | meanDiff around switch (matched N)','FontWeight','bold');
grid(ax,'on');

subplot(2,2,4);
MakeSpreadAndBoxPlot3_SB({TcM.meanDiff, ToM.meanDiff}, {cClassic, cOperant}, [1 2], {'classic','operant'}, ...
    'paired',1,'newfig',0,'showpoints',1);
set(gca,'Box','off','FontSize',FS);
yline(0,'k:','LineWidth',LWthin);
title('B4 | meanDiff matched (paired)','FontWeight','bold');
ylabel('meanDiff'); grid on

sgtitle(sprintf('B4 Mochi | %s | window=%s | classic vs operant (matched)', mk, winSel), ...
    'Interpreter','none','FontWeight','bold');
saveas(f5, fullfile(outRoot, sprintf('Mochi_%s_B4_motor_condition_matched_%s.png', mk, winSel)));
close(f5);

end

% ======================================================================
% Helpers
% ======================================================================
function Tout = filter_by_window(T, winSel, winMode)
w = lower(char(winSel));
Tout = T;
if ~ismember('Window', T.Properties.VariableNames), return; end

W = lower(string(T.Window));
if strcmpi(winMode,'exact')
    mask = (W == string(w));
else
    mask = contains(W, string(w));
end
Tout = T(mask,:);
end

function TaC = collapse_by_session_subset_window(Ta)

keys = strcat(string(Ta.Animal), "|", string(Ta.SessionName), "|", string(Ta.Subset), "|", string(Ta.Window));
[uk,~,ic] = unique(keys,'stable');

TaC = table();
TaC.Animal = cell(numel(uk),1);
TaC.SessionName = cell(numel(uk),1);
TaC.Subset = cell(numel(uk),1);
TaC.Window = cell(numel(uk),1);
TaC.SoundPair = cell(numel(uk),1);
TaC.BlockSessionIndex = nan(numel(uk),1);
TaC.MotorCondition = cell(numel(uk),1);

metrics = {'AUROC','meanDiff','meanTar','meanRef','stdTar','stdRef','CohensD','AshmanD','nA','nB','nTar','nRef'};

for mm = 1:numel(metrics)
    if ~ismember(metrics{mm}, Ta.Properties.VariableNames)
        Ta.(metrics{mm}) = nan(height(Ta),1);
    end
    TaC.(metrics{mm}) = nan(numel(uk),1);
end

for ii = 1:numel(uk)
    idx = (ic==ii);
    k0 = find(idx,1,'first');

    TaC.Animal{ii} = Ta.Animal{k0};
    TaC.SessionName{ii} = Ta.SessionName{k0};
    TaC.Subset{ii} = Ta.Subset{k0};
    TaC.Window{ii} = Ta.Window{k0};
    TaC.SoundPair{ii} = Ta.SoundPair{k0};
    TaC.BlockSessionIndex(ii) = min(Ta.BlockSessionIndex(idx), [], 'omitnan');
    TaC.MotorCondition{ii} = Ta.MotorCondition{k0};

    for mm = 1:numel(metrics)
        v = Ta.(metrics{mm})(idx);
        TaC.(metrics{mm})(ii) = mean(v, 'omitnan');
    end

    % Backfill nA/nB from nTar/nRef (per row) if needed
    if ismember('nA', TaC.Properties.VariableNames) && ~isfinite(TaC.nA(ii))
        if ismember('nTar', TaC.Properties.VariableNames)
            TaC.nA(ii) = TaC.nTar(ii);
        end
    end
    if ismember('nB', TaC.Properties.VariableNames) && ~isfinite(TaC.nB(ii))
        if ismember('nRef', TaC.Properties.VariableNames)
            TaC.nB(ii) = TaC.nRef(ii);
        end
    end
end

end

function [Sa, sessNames, x, blockChangeX, blockLabels, motorSwitchX] = make_session_axis(Sa, tag)

% Ensure DateNum exists (derive from SessionName if missing)
if ~ismember('DateNum', Sa.Properties.VariableNames)
    Sa.DateNum = ra_extract_datenum_from_sessionname(Sa.SessionName);
else
    % If present but empty/NaN, still try to backfill
    dn = Sa.DateNum;
    if ~isnumeric(dn), dn = nan(height(Sa),1); end
    need = ~isfinite(dn);
    if any(need)
        dn2 = ra_extract_datenum_from_sessionname(Sa.SessionName);
        dn(need) = dn2(need);
        Sa.DateNum = dn;
    end
end

% Sort chronologically
Sa = sortrows(Sa, {'DateNum','SessionName'});

sessNames = Sa.SessionName;
nSess = numel(sessNames);
x = (1:nSess)';

blockChangeX = [];
blockLabels  = {};
motorSwitchX = [];

if strcmpi(tag,'POOLED')
    return
end

if ismember('SoundPair', Sa.Properties.VariableNames)
    sp = Sa.SoundPair;
    for i = 2:nSess
        if ~strcmpi(sp{i}, sp{i-1})
            blockChangeX(end+1) = i - 0.5; %#ok<AGROW>
            blockLabels{end+1}  = sp{i}; %#ok<AGROW>
        end
    end
end

if ismember('MotorCondition', Sa.Properties.VariableNames)
    mc = Sa.MotorCondition;
    for i = 2:nSess
        if ~strcmpi(mc{i}, mc{i-1})
            motorSwitchX(end+1) = i - 0.5; %#ok<AGROW>
        end
    end
end
end

function dn = ra_extract_datenum_from_sessionname(sessionNames)
% Extract yyyymmdd from session name like "...20251202..." and convert to datenum.
sn = ra_to_cellstr(sessionNames);
dn = nan(numel(sn),1);

for i = 1:numel(sn)
    s = sn{i};
    if isempty(s), continue; end

    tok = regexp(s, '(?<y>\d{4})(?<m>\d{2})(?<d>\d{2})', 'names', 'once');
    if isempty(tok), continue; end

    y = str2double(tok.y);
    m = str2double(tok.m);
    d = str2double(tok.d);

    if isfinite(y) && isfinite(m) && isfinite(d)
        dn(i) = datenum(y,m,d);
    end
end
end


function y = series_by_session(Ts, sessNames, metricName)
n = numel(sessNames);
y = nan(n,1);
if isempty(Ts) || height(Ts)==0, return; end
if ~ismember(metricName, Ts.Properties.VariableNames), return; end
for i = 1:n
    idx = strcmpi_safe(Ts.SessionName, sessNames{i});
    if any(idx)
        v = Ts.(metricName)(idx);
        y(i) = mean(v, 'omitnan');
    end
end
end

function yl = ylims_across_subsets(TaC, sessNames, subsets, metricName, defaultYL)
vals = [];
for r = 1:numel(subsets)
    Ts = TaC(strcmpi_safe(TaC.Subset,subsets{r}), :);
    y = series_by_session(Ts, sessNames, metricName);
    vals = [vals; y(:)]; %#ok<AGROW>
end
vals = vals(isfinite(vals));
if isempty(vals)
    yl = defaultYL;
    return
end
mn = min(vals); mx = max(vals);
if mn==mx, mn = mn-0.05; mx = mx+0.05; end
pad = 0.08*(mx-mn);
yl = [mn-pad mx+pad];
if ~isempty(defaultYL)
    yl(1) = min(yl(1), defaultYL(1));
    yl(2) = max(yl(2), defaultYL(2));
end
end

function add_block_lines(ax, xLines, labels)
if isempty(xLines), return; end
yl = ylim(ax);
for i = 1:numel(xLines)
    line(ax, [xLines(i) xLines(i)], yl, 'Color',[0 0 0], 'LineStyle',':', 'LineWidth',1.2);
    if i <= numel(labels)
        text(xLines(i)+0.1, yl(2), labels{i}, 'FontSize',9, 'VerticalAlignment','top', 'Color',[0 0 0]);
    end
end
ylim(ax, yl);
end

function add_trend(ax, x, y, lw)
idx = isfinite(x) & isfinite(y);
if sum(idx) < 2, return; end
p = polyfit(x(idx), y(idx), 1);
yh = polyval(p, x(idx));
plot(ax, x(idx), yh, '-', 'Color',[0 0 0], 'LineWidth',lw);
end

function xt = nice_xticks(nSess)
if nSess <= 12
    xt = 1:nSess;
elseif nSess <= 25
    xt = 1:2:nSess;
elseif nSess <= 40
    xt = 1:3:nSess;
else
    xt = 1:5:nSess;
end
end

function [x,y] = xy_sorted(T, blockIdxMap, yName)
% Sort rows by within-block session index (computed from SessT), not by any
% potentially polluted BlockSessionIndex in the metrics table.
x = get_block_x(T, blockIdxMap);
y = T.(yName);
[~,ord] = sort(x);
x = x(ord);
y = y(ord);
end

function [firstVal, expertVal] = first_vs_expert(T, blockIdxMap, metricName)
y = T.(metricName);
x = get_block_x(T, blockIdxMap);
[~,ord] = sort(x);
y = y(ord);
idx = find(isfinite(y));
if isempty(idx)
    firstVal = NaN; expertVal = NaN; return
end
firstVal = y(idx(1));
k = min(3, numel(idx));
expertVal = mean(y(idx(end-k+1:end)), 'omitnan');
end

function [firstVal, expertVal, firstName, expertNames] = first_vs_expert_with_names(T, blockIdxMap, metricName)
y = T.(metricName);
x = get_block_x(T, blockIdxMap);
nm = T.SessionName;
[~,ord] = sort(x);
y = y(ord); nm = nm(ord);
idx = find(isfinite(y));
if isempty(idx)
    firstVal = NaN; expertVal = NaN;
    firstName = 'NA';
    expertNames = {};
    return
end
firstVal = y(idx(1));
firstName = nm{idx(1)};
k = min(3, numel(idx));
use = idx(end-k+1:end);
expertNames = nm(use);
expertVal = mean(y(use), 'omitnan');
end

function n = sessions_to_threshold(T, blockIdxMap, thr)
y = T.AUROC;
x = get_block_x(T, blockIdxMap);
[~,ord] = sort(x);
y = y(ord);
n = NaN;
idx = find(isfinite(y));
for i = 1:numel(idx)
    if y(idx(i)) >= thr
        n = idx(i);
        return
    end
end
end

function s = early_slope(T, blockIdxMap, K)
y = T.AUROC;
x = get_block_x(T, blockIdxMap);
[~,ord] = sort(x);
y = y(ord); x = x(ord);
idx = find(isfinite(y));
if numel(idx) < 2
    s = NaN; return
end
idx = idx(1:min(K, numel(idx)));
p = polyfit(x(idx), y(idx), 1);
s = p(1);
end

function x = get_block_x(T, blockIdxMap)
% Map each SessionName -> session index within its SoundPair block.
% If missing in the map, fall back to any numeric BlockSessionIndex.
n = height(T);
x = nan(n,1);
if isempty(T) || n==0, return; end

sess = string(T.SessionName);
for i = 1:n
    key = char(sess(i));
    if isKey(blockIdxMap, key)
        x(i) = blockIdxMap(key);
    end
end

if any(~isfinite(x)) && ismember('BlockSessionIndex', T.Properties.VariableNames)
    xb = T.BlockSessionIndex;
    if isnumeric(xb)
        x(~isfinite(x)) = xb(~isfinite(x));
    end
end
end

function C = gradient_multicolor(n)
if n <= 1
    C = [0.3 0.5 0.9];
    return
end
anchors = [
    0.90 0.93 0.98
    0.70 0.85 0.95
    0.45 0.75 0.80
    0.30 0.60 0.60
    0.40 0.45 0.75
    0.25 0.25 0.55
];
x_anchor = linspace(0,1,size(anchors,1));
x_query  = linspace(0,1,n);
C = interp1(x_anchor, anchors, x_query, 'pchip');
C = max(0,min(1,C));
end

function M = make_block_index_map(Sa)
% Build SessionName -> (session number within SoundPair) map.
% Sa must already be sorted chronologically (as returned by make_session_axis).
M = containers.Map('KeyType','char','ValueType','double');
if isempty(Sa) || height(Sa)==0, return; end
if ~ismember('SoundPair', Sa.Properties.VariableNames)
    % fallback: global index
    for i = 1:height(Sa)
        M(char(string(Sa.SessionName{i}))) = i;
    end
    return
end

sp = string(Sa.SoundPair);
sess = string(Sa.SessionName);
cnt = containers.Map('KeyType','char','ValueType','double');
for i = 1:height(Sa)
    ksp = char(sp(i));
    if ~isKey(cnt, ksp), cnt(ksp) = 0; end
    cnt(ksp) = cnt(ksp) + 1;
    M(char(sess(i))) = cnt(ksp);
end
end

function T = enforce_operant_sessions(T)

operantDates = {
    '20251128_n'
    '20251202_m'
    '20251202_n'
    '20251204_n'
    '20251209_m'
    '20251210_n'
    '20251211_m'
    '20251211_n'
    '20251212_n'
    '20251216_n'
    '20251217_n'
    '20251218_n'
    '20251219_m'
    '20251222_m'
    '20251224_n'
};

if ~ismember('MotorCondition', T.Properties.VariableNames)
    T.MotorCondition = repmat({''}, height(T), 1);
end
T.MotorCondition = ra_to_cellstr(T.MotorCondition);

% spelling fix
mask_operand = strcmpi_safe(T.MotorCondition, 'operand');
T.MotorCondition(mask_operand) = {'operant'};

if ~ismember('SessionName', T.Properties.VariableNames), return; end
sn = ra_to_cellstr(T.SessionName);

for i = 1:height(T)
    if isempty(sn{i}), continue; end
    for d = 1:numel(operantDates)
        if contains(sn{i}, operantDates{d})
            T.MotorCondition{i} = 'operant';
            break
        end
    end
end
end

function T = ra_ensure_plot_columns(T)
T = ra_alias(T, 'Marker',       {'Marker','marker'});
T = ra_alias(T, 'Animal',       {'Animal','animal'});
T = ra_alias(T, 'Subset',       {'Subset','subset'});
T = ra_alias(T, 'Window',       {'Window','window'});
T = ra_alias(T, 'BodyPart',     {'BodyPart','bodypart'});
T = ra_alias(T, 'SessionName',  {'SessionName','session'});
T = ra_alias(T, 'SoundPair',    {'SoundPair','soundPair','pair'});
T = ra_alias(T, 'MotorCondition', {'MotorCondition','motorCondition'});
if ~ismember('nA', T.Properties.VariableNames) && ismember('nTar', T.Properties.VariableNames)
    T.nA = T.nTar;
end
if ~ismember('nB', T.Properties.VariableNames) && ismember('nRef', T.Properties.VariableNames)
    T.nB = T.nRef;
end

if ~ismember('BlockSessionIndex', T.Properties.VariableNames)
    if ismember('SoundPair', T.Properties.VariableNames) && ismember('SessionName', T.Properties.VariableNames)
        T.BlockSessionIndex = ra_make_block_session_index(T.SoundPair, T.SessionName);
    else
        T.BlockSessionIndex = nan(height(T),1);
    end
end

txtCols = {'Marker','Animal','Subset','Window','SessionName','SoundPair','MotorCondition','BodyPart'};
for i = 1:numel(txtCols)
    c = txtCols{i};
    if ismember(c, T.Properties.VariableNames)
        T.(c) = ra_to_cellstr(T.(c));
    end
end
end

function SessT = ra_ensure_sess_columns(SessT)
% minimal normalization for SessT used here
if isempty(SessT) || height(SessT)==0, return; end

need = {'Animal','SessionName','SoundPair','MotorCondition'};
for i = 1:numel(need)
    if ~ismember(need{i}, SessT.Properties.VariableNames)
        SessT.(need{i}) = repmat({''}, height(SessT), 1);
    end
end

SessT.Animal = ra_to_cellstr(SessT.Animal);
SessT.SessionName = ra_to_cellstr(SessT.SessionName);
SessT.SoundPair = ra_to_cellstr(SessT.SoundPair);
SessT.MotorCondition = ra_to_cellstr(SessT.MotorCondition);

% Ensure DateNum exists (derived from SessionName)
if ~ismember('DateNum', SessT.Properties.VariableNames)
    SessT.DateNum = ra_extract_datenum_from_sessionname(SessT.SessionName);
else
    dn = SessT.DateNum;
    if ~isnumeric(dn), dn = nan(height(SessT),1); end
    need = ~isfinite(dn);
    if any(need)
        dn2 = ra_extract_datenum_from_sessionname(SessT.SessionName);
        dn(need) = dn2(need);
        SessT.DateNum = dn;
    end
end
end

function T = ra_alias(T, outName, candidates)
if ismember(outName, T.Properties.VariableNames)
    return
end
for k = 1:numel(candidates)
    nm = candidates{k};
    if ismember(nm, T.Properties.VariableNames)
        T.(outName) = T.(nm);
        return
    end
end
T.(outName) = repmat({''}, height(T), 1);
end

function cs = ra_to_cellstr(x)
if iscell(x), cs = x; return; end
if isstring(x), cs = cellstr(x); return; end
if iscategorical(x), cs = cellstr(x); return; end
if ischar(x), cs = repmat({x}, size(x,1), 1); return; end
cs = repmat({''}, numel(x), 1);
end

function idx = ra_make_block_session_index(soundPair, sessionName)
n = numel(soundPair);
idx = nan(n,1);
sp = ra_to_cellstr(soundPair);
u = unique(sp,'stable');
for b = 1:numel(u)
    m = strcmp(sp, u{b});
    ii = find(m);
    idx(ii) = (1:numel(ii))';
end
end

function tf = strcmpi_safe(a, b)
% a: cellstr (or convertible), b: char/cellstr
aa = ra_to_cellstr(a);
if iscell(b), bb = b; else, bb = {char(string(b))}; end
tf = false(numel(aa),1);
for i = 1:numel(aa)
    tf(i) = any(strcmpi(aa{i}, bb));
end
end
