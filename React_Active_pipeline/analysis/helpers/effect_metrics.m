function met = effect_metrics(xTar, xRef)

xTar = xTar(:);
xRef = xRef(:);

vTar = xTar(isfinite(xTar));
vRef = xRef(isfinite(xRef));

met = struct();
met.nA = numel(vTar);
met.nB = numel(vRef);

met.meanTar = nan; met.meanRef = nan;
met.stdTar  = nan; met.stdRef  = nan;

met.AshmanD  = nan;
met.AUROC    = nan;
met.CohenD   = nan;   % NEW canonical field
met.CohensD  = nan;   % backward-compatible alias
met.meanDiff = nan;

met.rocFPR = [];
met.rocTPR = [];

if met.nA < 3 || met.nB < 3
    return
end

met.meanTar = mean(vTar);
met.meanRef = mean(vRef);
met.stdTar  = std(vTar);
met.stdRef  = std(vRef);

met.meanDiff = met.meanTar - met.meanRef;

vr = var(vRef);
vt = var(vTar);

if isfinite(vr) && isfinite(vt) && (vr+vt) > 0
    met.AshmanD = abs(met.meanRef - met.meanTar) / sqrt(0.5*(vr+vt));
end

den = (met.nB + met.nA - 2);
if den > 0
    sPool = sqrt(((met.nB-1)*vr + (met.nA-1)*vt) / den);
    if isfinite(sPool) && sPool > 0
        d = (met.meanTar - met.meanRef) / sPool;
        met.CohenD  = d;
        met.CohensD = d;
    end
end

allv = [vRef; vTar];
if numel(unique(allv)) > 1
    lab = [zeros(numel(vRef),1); ones(numel(vTar),1)];
    try
        [fpr,tpr,~,auc] = perfcurve(lab, allv, 1);
        met.AUROC = auc;
        met.rocFPR = fpr;
        met.rocTPR = tpr;
    catch
        % perfcurve may be unavailable
    end
end

end
