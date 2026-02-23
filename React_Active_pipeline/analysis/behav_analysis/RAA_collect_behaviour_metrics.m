function BigT = RAA_collect_behaviour_metrics(SessT, opts)
% RAA_collect_behaviour_metrics
% Build a long table of metrics across sessions for a given bodypart.
%
% REQUIRED INPUTS
%   SessT : from RAA_get_training_sessions
%   opts.bodypart : e.g. 'Pupil-EyeCam'
%
% OPTIONAL
%   opts.subsets : {'regular','nosound','nomotor'} (default)
%   opts.window_patterns : cellstr patterns to keep
%       default for your case: {'stimoff_to_arrival'}
%   opts.SW_list : numeric array of SW values to keep, [] = keep all
%   opts.markers : cell array of marker names to keep, {} = keep all found

if ~isfield(opts,'bodypart'), error('opts.bodypart is required'); end

if ~isfield(opts,'subsets') || isempty(opts.subsets)
    opts.subsets = {'regular','nosound','nomotor'};
end
if ~isfield(opts,'window_patterns') || isempty(opts.window_patterns)
    % your default: [stim off : arrival]
    opts.window_patterns = {'stimoff_to_arrival'};
end
if ~isfield(opts,'SW_list'),  opts.SW_list  = []; end
if ~isfield(opts,'markers'),  opts.markers  = {}; end

varNames = { ...
    'Animal','SoundPair','BlockIndex','BlockSessionIndex','SessionName','SessionPath','DateNum','GlobalSessionIndex', ...
    'MotorCondition','Bodypart','Subset','Window','SW','Marker', ...
    'AUROC','meanDiff','meanTar','meanRef','stdTar','stdRef','CohensD','AshmanD','nA','nB'};

nCols = numel(varNames);
rows  = cell(0, nCols);

badRowCount = 0;

for s = 1:height(SessT)
    spath = SessT.SessionPath{s};
    metricsDir = fullfile(spath, 'analysis', 'behaviour', opts.bodypart);
    if ~exist(metricsDir,'dir')
        continue
    end

    D = dir(fullfile(metricsDir, '*_metrics.mat'));
    if isempty(D), continue; end

    for f = 1:numel(D)
        fname = D(f).name;
        fpath = fullfile(D(f).folder, fname);

        % subset
        subset = parse_subset_from_name(fname, opts.subsets);
        if isempty(subset), continue; end

        % SW
        SW = parse_sw_from_name(fname);
        if ~isempty(opts.SW_list)
            if ~isfinite(SW) || ~any(abs(opts.SW_list - SW) < 1e-12)
                continue
            end
        end

        % window string (raw token chunk)
        win = parse_window_from_name(fname, subset);
        if strlength(win)==0, win = "unknown"; end

        % window filter
        if ~window_matches_patterns(win, opts.window_patterns)
            continue
        end

        % load metrics
        L = load(fpath);
        if ~isfield(L,'allMet') || isempty(L.allMet)
            continue
        end
        allMet = L.allMet;

        for m = 1:numel(allMet)
            if ~isfield(allMet(m),'marker') || ~isfield(allMet(m),'metrics')
                continue
            end
            marker = string(allMet(m).marker);

            if ~isempty(opts.markers) && ~any(strcmp(marker, string(opts.markers)))
                continue
            end

            M = allMet(m).metrics;
            M = normalize_metrics_struct(M);
            M = ensure_metric_fields(M);

            row = cell(1, nCols);
            row(:) = {[]};

            row{1}  = SessT.Animal{s};
            row{2}  = SessT.SoundPair{s};
            row{3}  = SessT.BlockIndex(s);
            row{4}  = SessT.BlockSessionIndex(s);
            row{5}  = SessT.SessionName{s};
            row{6}  = SessT.SessionPath{s};
            row{7}  = SessT.DateNum(s);
            row{8}  = SessT.GlobalSessionIndex(s);
            row{9}  = char(SessT.MotorCondition(s));   % R2018b-safe
            row{10} = opts.bodypart;
            row{11} = subset;
            row{12} = char(win);
            row{13} = SW;
            row{14} = char(marker);

            row{15} = M.AUROC;
            row{16} = M.meanDiff;
            row{17} = M.meanTar;
            row{18} = M.meanRef;
            row{19} = M.stdTar;
            row{20} = M.stdRef;
            row{21} = M.CohensD;
            row{22} = M.AshmanD;
            row{23} = M.nA;
            row{24} = M.nB;

            if numel(row) ~= nCols
                badRowCount = badRowCount + 1;
                if badRowCount <= 10
                    fprintf('BAD ROW (%d): %s | %s | %s\n', badRowCount, spath, fname, marker);
                end
                continue
            end

            rows(end+1,:) = row; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    BigT = cell2table(cell(0,nCols), 'VariableNames', varNames);
    return
end

BigT = cell2table(rows, 'VariableNames', varNames);

end

function subset = parse_subset_from_name(fname, subsetList)
subset = '';
low = lower(fname);
for i = 1:numel(subsetList)
    if contains(low, lower(subsetList{i}))
        subset = subsetList{i};
        return
    end
end
end

function SW = parse_sw_from_name(fname)
tok = regexp(fname, 'sw([0-9]+(\.[0-9]+)?)', 'tokens', 'once');
if isempty(tok), SW = NaN; return; end
SW = str2double(tok{1});
end

function win = parse_window_from_name(fname, subset)
% window is whatever sits between subset and "_sw"
low = fname;
p1 = strfind(lower(low), ['_' lower(subset) '_']);
p2 = strfind(lower(low), '_sw');
if isempty(p1) || isempty(p2) || p2(1) <= p1(1)
    win = "";
    return
end
startIdx = p1(1) + numel(subset) + 2;
endIdx   = p2(1) - 1;
win = string(low(startIdx:endIdx));
end

function ok = window_matches_patterns(win, patterns)
% If pattern contains '_' (like 'stimoff_to_arrival'), enforce exact match.
% Otherwise keep "contains" logic for flexible matching.
w = lower(char(win));
ok = false;
for i = 1:numel(patterns)
    p = lower(patterns{i});
    if contains(p, '_')
        if strcmp(w, p)
            ok = true;
            return
        end
    else
        if contains(w, p)
            ok = true;
            return
        end
    end
end
end

function M = normalize_metrics_struct(M)
% Make sure M is a 1x1 struct. Convert table/cell/struct array to 1x1.
if iscell(M)
    if isempty(M), M = struct(); return; end
    M = M{1};
end
if istable(M)
    M = table2struct(M);
end
if isstruct(M) && numel(M) > 1
    M = M(1);
end
if ~isstruct(M)
    M = struct();
end

% Force scalar numeric values; empty numeric -> NaN
fn = fieldnames(M);
for i = 1:numel(fn)
    v = M.(fn{i});
    if isnumeric(v)
        if isempty(v)
            M.(fn{i}) = NaN;
        elseif ~isscalar(v)
            M.(fn{i}) = v(1);
        end
    end
end
end

function M = ensure_metric_fields(M)
need = {'AUROC','meanDiff','meanTar','meanRef','stdTar','stdRef','CohensD','AshmanD','nA','nB'};
for i = 1:numel(need)
    if ~isfield(M, need{i}) || (isnumeric(M.(need{i})) && isempty(M.(need{i})))
        M.(need{i}) = NaN;
    end
end
end

