function BigT = RAA_collect_all_behaviour_tables(csv_files, outdir)
% RAA_merge_BigT_csvs
% Merge multiple BigT_*.csv tables with possibly different column sets
% (e.g., different numbers of rocFPR_/rocTPR_ columns).
%
% INPUT
%   csv_files : cell array of csv file paths
%   outdir    : output directory (default: pwd)
%
% OUTPUT
%   BigT : merged table (also saved to outdir)
%
% Saves:
%   BigT_allmarkers.mat  (full table, includes ROC columns)
%   BigT_allmarkers_scalar.csv (drops rocFPR_/rocTPR_ columns)

if nargin < 2 || isempty(outdir)
    outdir = pwd;
end
if ~exist(outdir,'dir')
    mkdir(outdir);
end
if ~iscell(csv_files)
    error('csv_files must be a cell array of file paths');
end

BigT = table();
allVars = {};

% 1) First pass: collect union of variable names
for i = 1:numel(csv_files)
    f = csv_files{i};
    if ~exist(f,'file')
        error('Missing file: %s', f);
    end
    T = readtable(f);
    allVars = union(allVars, T.Properties.VariableNames, 'stable');
end

% 2) Second pass: align and concatenate
for i = 1:numel(csv_files)
    f = csv_files{i};
    T = readtable(f);

    % Drop __MEAN__ always
    if ismember('marker', T.Properties.VariableNames)
        T = T(~strcmpi(T.marker,'__MEAN__'), :);
    end

    % Add missing variables
    miss = setdiff(allVars, T.Properties.VariableNames);
    for k = 1:numel(miss)
        vn = miss{k};
        % default: numeric NaN column
        T.(vn) = nan(height(T),1);
    end

    % Ensure same column order
    T = T(:, allVars);

    BigT = [BigT; T]; %#ok<AGROW>
    fprintf('Merged %d rows from %s\n', height(T), f);
end

% Save full table
save(fullfile(outdir,'BigT_allmarkers.mat'), 'BigT', '-v7.3');

% Save scalar-only CSV (drop ROC curves)
isRoc = startsWith(BigT.Properties.VariableNames,'rocFPR_') | ...
        startsWith(BigT.Properties.VariableNames,'rocTPR_');
BigT_scalar = BigT(:, ~isRoc);

writetable(BigT_scalar, fullfile(outdir,'BigT_allmarkers_scalar.csv'));

fprintf('Saved:\n  %s\n  %s\n', ...
    fullfile(outdir,'BigT_allmarkers.mat'), ...
    fullfile(outdir,'BigT_allmarkers_scalar.csv'));
end
