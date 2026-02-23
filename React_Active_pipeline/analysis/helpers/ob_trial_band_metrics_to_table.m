function T = ob_trial_band_metrics_to_table(sessname, alignName, idx_used, trialType, bandRanges, bandNames, winNames, valMat, slpMat)
% Build long-format per-trial table from eventlocked_spec_grid_plot outputs.
% valMat/slpMat: [nTrialsUsed x nBands x nWins]

nT = size(valMat,1);
nB = size(valMat,2);
nW = size(valMat,3);

rows = nT*nB*nW;

session = repmat(string(sessname), rows, 1);
align   = repmat(string(alignName), rows, 1);

trialIdx   = zeros(rows,1);
trialTypeS = strings(rows,1);
band       = strings(rows,1);
bandLo     = nan(rows,1);
bandHi     = nan(rows,1);
window     = strings(rows,1);
meanPower  = nan(rows,1);
slope      = nan(rows,1);

r = 0;
for i = 1:nT
    for b = 1:nB
        for w = 1:nW
            r = r + 1;
            trialIdx(r)   = idx_used(i);
            trialTypeS(r) = string(trialType(i));
            band(r)       = string(bandNames{b});
            bandLo(r)     = bandRanges(b,1);
            bandHi(r)     = bandRanges(b,2);
            window(r)     = string(winNames{w});
            meanPower(r)  = valMat(i,b,w);
            slope(r)      = slpMat(i,b,w);
        end
    end
end

T = table(session, align, trialIdx, trialTypeS, band, bandLo, bandHi, window, meanPower, slope, ...
    'VariableNames', {'session','align','trialIdx','trialType','band','bandLo','bandHi','window','meanPower','slope'});
end
