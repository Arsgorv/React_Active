function T = RAA_make_OB_as_bodypart_metrics(obTrialCsvLow, obTrialCsvMid, outCsv)
% Builds a behaviour_analysis-like metrics table for OB bands using the same
% columns as your Pupil-EyeCam metrics.csv (bodypart/marker/subset/window/AUROC/etc).
%
% Markers: delta/theta/gamma/highgamma
% Window: stimOn_to_stimOff, stimOff_to_arr, arr_to_stop
% Subset: "TarRef" (only compares Target vs Reference)

Tlow = readtable(obTrialCsvLow);
Tmid = readtable(obTrialCsvMid);
Tall = [Tlow; Tmid];

Tall.trialType = string(Tall.trialType);
Tall.band      = string(Tall.band);
Tall.window    = string(Tall.window);

bands = unique(Tall.band,'stable');
wins  = unique(Tall.window,'stable');

rows = {};
for bi = 1:numel(bands)
    b = bands(bi);
    for wi = 1:numel(wins)
        w = wins(wi);

        Tk = Tall(Tall.band==b & Tall.window==w,:);
        xTar = Tk.meanPower(contains(lower(Tk.trialType),'target'));
        xRef = Tk.meanPower(contains(lower(Tk.trialType),'ref'));

        xTar = xTar(isfinite(xTar));
        xRef = xRef(isfinite(xRef));

        if numel(xTar) < 5 || numel(xRef) < 5
            continue
        end

        M = effect_metrics(xTar, xRef); % your existing helper used elsewhere

        % Flatten ROC to fixed-length columns like behaviour tables do
        rocFPR = M.rocFPR(:)'; rocTPR = M.rocTPR(:)';
        nR = numel(rocFPR);

        % Build row
        base = { "OB", char(b), "TarRef", char(w), ...
            numel(xTar), numel(xRef), M.AshmanD, M.AUROC, M.CohensD, M.meanDiff, ...
            mean(xRef), mean(xTar), std(xRef), std(xTar) };

        % Append roc columns
        rocCols = cell(1, 2*nR);
        for i = 1:nR
            rocCols{i} = rocFPR(i);
            rocCols{nR+i} = rocTPR(i);
        end

        rows(end+1,:) = [base, rocCols]; %#ok<AGROW>
    end
end

% Build variable names
% Match your existing behaviour metrics.csv convention
if isempty(rows)
    T = table(); return
end

nR = (size(rows,2) - 14)/2;
vn = {'bodypart','marker','subset','window','nTar','nRef','AshmanD','AUROC','CohensD','meanDiff','meanRef','meanTar','stdRef','stdTar'};
for i = 1:nR, vn{end+1} = sprintf('rocFPR_%d', i); end %#ok<AGROW>
for i = 1:nR, vn{end+1} = sprintf('rocTPR_%d', i); end %#ok<AGROW>

T = cell2table(rows, 'VariableNames', vn);

if nargin >= 3 && ~isempty(outCsv)
    writetable(T, outCsv);
end
end
