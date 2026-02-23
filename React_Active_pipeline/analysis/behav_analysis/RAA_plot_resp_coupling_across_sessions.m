function RAA_plot_resp_coupling_across_sessions(RES, opts)
% Plot respiration/sniff–OB coupling across sessions.
% Uses RES.Ssess if present (session summary); otherwise derives it from RES.Ttrial.
%
% Figures:
% 1) meanR across sessions (band x window grid; separate figure per piezoMode if available)
% 2) medianLag_s across sessions (same layout)
%
% Expected columns in RES.Ssess:
% session, band, window, meanR, medianLag_s, nUsed
% Optional:
% piezoMode (resp/sniff), group, spec

if nargin < 2, opts = struct(); end
if ~isfield(opts,'minN'), opts.minN = 20; end
if ~isfield(opts,'alignName'), opts.alignName = 'stimOn'; end
if ~isfield(opts,'showN'), opts.showN = true; end
if ~isfield(opts,'plotIndividualWindows'), opts.plotIndividualWindows = true; end

S = RES.Ssess;
if isempty(S) && ~isempty(RES.Ttrial)
    S = ra_compute_resp_session_summary_from_trials_fixed(RES.Ttrial);
end
if isempty(S)
    warning('No coupling summary available.');
    return
end

% normalize strings
if ismember('session', S.Properties.VariableNames), S.session = string(S.session); end
if ismember('band', S.Properties.VariableNames),    S.band    = string(S.band); end
if ismember('window', S.Properties.VariableNames),  S.window  = string(S.window); end
if ismember('piezoMode', S.Properties.VariableNames), S.piezoMode = string(S.piezoMode); end

% session order from SessInfo (if exists)
sessOrder = unique(S.session,'stable');
if isfield(RES,'SessInfo') && ~isempty(RES.SessInfo) && ismember('session', RES.SessInfo.Properties.VariableNames)
    so = string(RES.SessInfo.session);
    sessOrder = so;
end

% piezo modes
if ismember('piezoMode', S.Properties.VariableNames)
    modes = unique(S.piezoMode,'stable');
else
    modes = "resp";
    S.piezoMode = repmat(modes, height(S), 1);
end

bandsU = unique(S.band,'stable');
winsU  = unique(S.window,'stable');

% enforce canonical window ordering if your names match
canonWins = ["stimOn_to_stimOff","stimOff_to_arr","arr_to_stop"];
if all(ismember(canonWins, winsU))
    winsU = canonWins;
end

% helper to reorder per-session rows
sess_index = @(x) arrayfun(@(z) find(sessOrder==z,1,'first'), x);

for mi = 1:numel(modes)
    modeName = modes(mi);
    Sm = S(strcmpi(S.piezoMode, modeName), :);
    if isempty(Sm), continue, end

    % ---- meanR dynamics
    figure('Color','w','Units','pixels','Position',[50 50 1750 900], ...
        'Name', sprintf('RespCoupling_meanR_%s_%s', char(modeName), opts.alignName));

    nR = numel(bandsU);
    nC = numel(winsU);

    for bi = 1:nR
        for wi = 1:nC
            ax = subplot(nR, nC, (bi-1)*nC + wi);

            T = Sm(strcmpi(Sm.band, bandsU(bi)) & strcmpi(Sm.window, winsU(wi)), :);
            if isempty(T), axis(ax,'off'); continue, end

            % order by session
            ord = sess_index(T.session);
            [~,ix] = sort(ord);
            T = T(ix,:);

            y = T.meanR;
            nUsed = nan(height(T),1);
            if ismember('nUsed', T.Properties.VariableNames), nUsed = T.nUsed; end

            % mask low-N
            good = isfinite(y);
            if isfinite(opts.minN) && ismember('nUsed',T.Properties.VariableNames)
                good = good & (nUsed >= opts.minN);
            end

            plot(ax, find(good), y(good), 'k-', 'LineWidth', 1); hold(ax,'on');
            plot(ax, find(good), y(good), 'k.', 'MarkerSize', 10);

            yline(ax, 0, 'k:');
            xlabel(ax,'session #'); ylabel(ax,'meanR');
            title(ax, sprintf('%s | %s | %s', char(modeName), char(bandsU(bi)), char(winsU(wi))), 'Interpreter','none');
            box(ax,'off');

            if opts.showN && ismember('nUsed',T.Properties.VariableNames)
                % annotate median N
                txt = sprintf('median n=%d', round(median(nUsed(good),'omitnan')));
                xlimv = xlim(ax); ylimv = ylim(ax);
                text(ax, xlimv(1)+0.02*range(xlimv), ylimv(2)-0.08*range(ylimv), txt, 'FontSize', 8, 'Interpreter','none');
            end
            hold(ax,'off');
        end
    end

    % ---- medianLag dynamics
    figure('Color','w','Units','pixels','Position',[50 50 1750 900], ...
        'Name', sprintf('RespCoupling_medianLag_%s_%s', char(modeName), opts.alignName));

    for bi = 1:nR
        for wi = 1:nC
            ax = subplot(nR, nC, (bi-1)*nC + wi);

            T = Sm(strcmpi(Sm.band, bandsU(bi)) & strcmpi(Sm.window, winsU(wi)), :);
            if isempty(T), axis(ax,'off'); continue, end

            ord = sess_index(T.session);
            [~,ix] = sort(ord);
            T = T(ix,:);

            y = T.medianLag_s;
            nUsed = nan(height(T),1);
            if ismember('nUsed', T.Properties.VariableNames), nUsed = T.nUsed; end

            good = isfinite(y);
            if isfinite(opts.minN) && ismember('nUsed',T.Properties.VariableNames)
                good = good & (nUsed >= opts.minN);
            end

            plot(ax, find(good), y(good), 'k-', 'LineWidth', 1); hold(ax,'on');
            plot(ax, find(good), y(good), 'k.', 'MarkerSize', 10);

            yline(ax, 0, 'k:');
            xlabel(ax,'session #'); ylabel(ax,'medianLag (s)');
            title(ax, sprintf('%s | %s | %s', char(modeName), char(bandsU(bi)), char(winsU(wi))), 'Interpreter','none');
            box(ax,'off');

            if opts.showN && ismember('nUsed',T.Properties.VariableNames)
                txt = sprintf('median n=%d', round(median(nUsed(good),'omitnan')));
                xlimv = xlim(ax); ylimv = ylim(ax);
                text(ax, xlimv(1)+0.02*range(xlimv), ylimv(2)-0.08*range(ylimv), txt, 'FontSize', 8, 'Interpreter','none');
            end
            hold(ax,'off');
        end
    end
end

end

% --- fixed version of your helper (your earlier code had T vs Ttrial mixups)
function S = ra_compute_resp_session_summary_from_trials_fixed(Ttrial)
T = Ttrial;

T.session = string(T.session);
if ismember('piezoMode', T.Properties.VariableNames), T.piezoMode = string(T.piezoMode); else, T.piezoMode = repmat("resp", height(T),1); end
if ismember('band', T.Properties.VariableNames),     T.band = string(T.band); end
if ismember('window', T.Properties.VariableNames),   T.window = string(T.window); end

% group by session x piezoMode x band x window
[G, sSess, sMode, sBand, sWin] = findgroups(T.session, T.piezoMode, T.band, T.window);

S = table();
S.session   = sSess;
S.piezoMode = sMode;
S.band      = sBand;
S.window    = sWin;

S.nUsed = splitapply(@(x) sum(isfinite(x)), T.peakR, G);
S.medianLag_s = splitapply(@(x) median(x,'omitnan'), T.peakLag_s, G);
S.iqrLag_s = splitapply(@(x) iqr(x(isfinite(x))), T.peakLag_s, G);
S.meanR  = splitapply(@(x) mean(x,'omitnan'), T.peakR, G);
S.medianR= splitapply(@(x) median(x,'omitnan'), T.peakR, G);
end