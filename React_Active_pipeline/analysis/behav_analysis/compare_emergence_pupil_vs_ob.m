function OUT = compare_emergence_pupil_vs_ob(pupilSessTable, obSessTable, opts)
% Compare emergence session index for pupil vs OB marker.
% Emergence criterion: AUROC >= thr for kConsec consecutive sessions.
%
% Inputs:
% - pupilSessTable: columns {session, AUROC} at least, one row per session (ordered)
% - obSessTable: columns {session, band, window, AUROC}
%
% opts fields:
%   .thr = 0.6
%   .kConsec = 2
%   .ob_band = "delta"
%   .ob_window = "stimOff_to_arr"
%
% Output: OUT with emergence indices and plotted curves.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'thr'), opts.thr = 0.6; end
if ~isfield(opts,'kConsec'), opts.kConsec = 2; end
if ~isfield(opts,'ob_band'), opts.ob_band = "delta"; end
if ~isfield(opts,'ob_window'), opts.ob_window = "stimOff_to_arr"; end

% Order sessions by appearance in pupil table
sessOrder = string(pupilSessTable.session(:));

% Pupil AUROC vector
pA = nan(numel(sessOrder),1);
for i = 1:numel(sessOrder)
    r = pupilSessTable(string(pupilSessTable.session)==sessOrder(i),:);
    if ~isempty(r), pA(i) = r.AUROC(1); end
end

% OB AUROC vector (select band+window)
obSel = obSessTable(string(obSessTable.band)==string(opts.ob_band) & string(obSessTable.window)==string(opts.ob_window),:);

oA = nan(numel(sessOrder),1);
for i = 1:numel(sessOrder)
    r = obSel(string(obSel.session)==sessOrder(i),:);
    if ~isempty(r), oA(i) = r.AUROC(1); end
end

pEmerg = first_consec_crossing(pA, opts.thr, opts.kConsec);
oEmerg = first_consec_crossing(oA, opts.thr, opts.kConsec);

OUT = struct();
OUT.sessionOrder = sessOrder;
OUT.pupil_AUROC = pA;
OUT.ob_AUROC = oA;
OUT.pupil_emergence_idx = pEmerg;
OUT.ob_emergence_idx = oEmerg;

figure('Color','w'); hold on
plot(pA,'LineWidth',1);
plot(oA,'LineWidth',1);
yline(0.5,'k:'); yline(opts.thr,'k--');
xlabel('session #'); ylabel('AUROC');
legend({'pupil','OB'},'Box','off','Location','southoutside');
title(sprintf('Emergence: pupil=%s, OB(%s,%s)=%s', num2str(pEmerg), char(opts.ob_band), char(opts.ob_window), num2str(oEmerg)), 'Interpreter','none');
box off
hold off

end

function idx = first_consec_crossing(x, thr, k)
idx = NaN;
if numel(x) < k, return, end
for i = 1:(numel(x)-k+1)
    if all(x(i:i+k-1) >= thr)
        idx = i;
        return
    end
end
end
