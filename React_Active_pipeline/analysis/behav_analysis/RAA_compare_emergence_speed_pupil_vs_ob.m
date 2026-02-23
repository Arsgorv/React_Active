function OUT = RAA_compare_emergence_speed_pupil_vs_ob(RESob, BRES, opts)
% Compare emergence speed of pupil AUROC vs OB AUROC across sessions.
%
% Inputs:
% RESob.Sess: table with session, band, window, AUROC, nTar, nRef
% BRES.T:     table with session, AUROC, nTar, nRef (already filtered to one marker/subset/window)
%
% opts:
%   .bandList   = {'delta','theta','gamma','highgamma'}  % OB bands to test
%   .window     = 'stimoff_to_arrival'
%   .thr        = 0.60
%   .kConsec    = 2
%   .minN       = 10
%   .plot       = true

if nargin < 3, opts = struct(); end
if ~isfield(opts,'bandList'), opts.bandList = {'delta','theta','gamma','highgamma'}; end
if ~isfield(opts,'window'),   opts.window   = 'stimoff_to_arrival'; end
if ~isfield(opts,'thr'),      opts.thr      = 0.60; end
if ~isfield(opts,'kConsec'),  opts.kConsec  = 2; end
if ~isfield(opts,'minN'),     opts.minN     = 10; end
if ~isfield(opts,'plot'),     opts.plot     = true; end

S = RESob.Sess;
Tb = BRES.T;

if isempty(S) || isempty(Tb)
    error('Empty RESob.Sess or BRES.T');
end

S.session = string(S.session);
S.band    = string(S.band);
S.window  = string(S.window);

Tb.session = string(Tb.session);

% session ordering
sessOrder = unique(S.session,'stable');
if isfield(RESob,'SessInfo') && ~isempty(RESob.SessInfo) && ismember('session', RESob.SessInfo.Properties.VariableNames)
    sessOrder = string(RESob.SessInfo.session);
end

% pupil vector over sessions (aligned to sessOrder)
pAU = nan(numel(sessOrder),1);
pNT = nan(numel(sessOrder),1);
pNR = nan(numel(sessOrder),1);

for i = 1:numel(sessOrder)
    m = (Tb.session == sessOrder(i));
    if any(m)
        r = Tb(find(m,1,'first'),:);
        pAU(i) = r.AUROC;
        if ismember('nTar',Tb.Properties.VariableNames), pNT(i) = r.nTar; end
        if ismember('nRef',Tb.Properties.VariableNames), pNR(i) = r.nRef; end
    end
end

% enforce minN
goodP = isfinite(pAU);
if ismember('nTar',Tb.Properties.VariableNames) && ismember('nRef',Tb.Properties.VariableNames)
    goodP = goodP & (pNT >= opts.minN) & (pNR >= opts.minN);
end
pAU(~goodP) = nan;

pEmerg = ra_first_consec_crossing(pAU, opts.thr, opts.kConsec);

% OB bands emergence
rows = {};
for bi = 1:numel(opts.bandList)
    b = string(opts.bandList{bi});

    T = S(strcmpi(S.band,b) & strcmpi(S.window,string(opts.window)), :);

    % map AUROC to sessOrder
    x = nan(numel(sessOrder),1);
    nT = nan(numel(sessOrder),1);
    nR = nan(numel(sessOrder),1);

    for i = 1:numel(sessOrder)
        m = (T.session == sessOrder(i));
        if any(m)
            r = T(find(m,1,'first'),:);
            x(i) = r.AUROC;
            if ismember('nTar',T.Properties.VariableNames), nT(i) = r.nTar; end
            if ismember('nRef',T.Properties.VariableNames), nR(i) = r.nRef; end
        end
    end

    good = isfinite(x);
    if ismember('nTar',T.Properties.VariableNames) && ismember('nRef',T.Properties.VariableNames)
        good = good & (nT >= opts.minN) & (nR >= opts.minN);
    end
    x(~good) = nan;

    bEmerg = ra_first_consec_crossing(x, opts.thr, opts.kConsec);

    rows(end+1,:) = {b, bEmerg, pEmerg, bEmerg - pEmerg}; %#ok<AGROW>

    if opts.plot
        % per-band overlay plot
        figure('Color','w','Units','pixels','Position',[200 120 900 320], ...
            'Name', sprintf('Emergence_%s_%s', char(b), char(opts.window)));
        plot(pAU,'LineWidth',1.5); hold on
        plot(x,'LineWidth',1.5);
        yline(opts.thr,'k:');
        yline(0.5,'k:');
        xlabel('session #'); ylabel('AUROC');
        legend({'Pupil','OB'},'Location','best','Box','off');
        title(sprintf('Window=%s | Band=%s | thr=%.2f k=%d | pupil=%d ob=%d', ...
            char(opts.window), char(b), opts.thr, opts.kConsec, pEmerg, bEmerg), 'Interpreter','none');
        ylim([0 1]); box off
        hold off
    end
end

Tout = cell2table(rows, 'VariableNames', {'band','obEmergenceSession','pupilEmergenceSession','obMinusPupil'});

OUT = struct();
OUT.sessionOrder = sessOrder;
OUT.pupilAUROC = pAU;
OUT.pupilEmergenceSession = pEmerg;
OUT.table = Tout;
OUT.opts_snapshot = opts;

% summary figure
if opts.plot
    figure('Color','w','Units','pixels','Position',[120 120 900 420], ...
        'Name', sprintf('EmergenceSummary_%s', char(opts.window)));

    % bar: obMinusPupil
    vals = Tout.obMinusPupil;
    x = 1:numel(vals);
    bar(x, vals); hold on
    yline(0,'k:');
    set(gca,'XTick',x,'XTickLabel',cellstr(Tout.band));
    ylabel('Emergence session (OB - Pupil)');
    title(sprintf('Positive => pupil earlier | thr=%.2f k=%d window=%s', opts.thr, opts.kConsec, char(opts.window)), 'Interpreter','none');
    box off
    hold off
end

end

function idx = ra_first_consec_crossing(x, thr, kConsec)
% Returns first index i where x(i:i+kConsec-1) are all >= thr.
% NaNs break runs. Returns NaN if never crosses.
idx = nan;
if isempty(x) || all(~isfinite(x)), return, end
n = numel(x);
for i = 1:(n-kConsec+1)
    seg = x(i:(i+kConsec-1));
    if all(isfinite(seg)) && all(seg >= thr)
        idx = i;
        return
    end
end
end