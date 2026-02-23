function RAA_run_behaviour_across_sessions(remove_sess, opts)
% RAA_run_behaviour_across_sessions
%
% Close-to-original style:
% - Build BigT for ONE selected marker (fast).
% - Optionally call your existing plotter with minimal assumptions.
%
% opts.marker : char (default 'Pupil-EyeCam')
% opts.use_mat: true/false (default true)
% opts.do_plot: true/false (default true)

if nargin < 1, remove_sess = {}; end
if nargin < 2, opts = struct(); end
if ~isfield(opts,'marker') || isempty(opts.marker), opts.marker = 'Pupil-EyeCam'; end
if ~isfield(opts,'use_mat'), opts.use_mat = true; end
if ~isfield(opts,'do_plot'), opts.do_plot = true; end

SessT = RAA_get_training_sessions(remove_sess);
if isempty(SessT) || height(SessT)==0
    warning('RAA_run_behaviour_across_sessions:NoSessions','No sessions returned by RAA_get_training_sessions.');
    return
end
if ~ismember('SessionPath', SessT.Properties.VariableNames)
    error('SessT must contain SessionPath');
end

% Collect metrics (ONE marker)
BigT = table();
for i = 1:height(SessT)
    dp = SessT.SessionPath{i};
    if ~ischar(dp) && ~isstring(dp), continue; end
    dp = char(string(dp));
    if ~exist(dp,'dir'), continue; end

    % load ONLY this marker's metrics from this session
    T = ra_load_session_metrics_marker(dp, opts.marker, opts.use_mat);

    if isempty(T) || height(T)==0
        continue
    end

    % Ensure session bookkeeping columns exist (stable across old/new tables)
    if ~ismember('session_path', T.Properties.VariableNames)
        T.session_path = repmat({dp}, height(T), 1);
    end
    if ~ismember('session', T.Properties.VariableNames)
        [~,sn] = fileparts(dp);
        T.session = repmat({sn}, height(T), 1);
    end

    BigT = ra_vertcat_align(BigT, T);
end

if isempty(BigT) || height(BigT)==0
    warning('RAA_run_behaviour_across_sessions:NoData','No metrics found for marker %s.', opts.marker);
    return
end

BigT = ra_join_session_metadata(BigT, SessT);

% Safety net: remove exact duplicates by key
BigT = ra_alias_table(BigT, 'Marker', {'Marker','marker'});
BigT = ra_alias_table(BigT, 'Subset', {'Subset','subset'});
BigT = ra_alias_table(BigT, 'Window', {'Window','window'});

k = strcat(string(BigT.session_path),"|",string(BigT.Marker),"|",string(BigT.Subset),"|",string(BigT.Window));
[~, ia] = unique(k, 'stable');
if numel(ia) < height(BigT)
    BigT = BigT(ia,:);
end

% Save aggregate
baseRoot = fileparts(SessT.SessionPath{1});
outRoot = fullfile(baseRoot, 'DS_figures', 'behaviour', ra_safename(opts.marker));
if ~exist(outRoot,'dir'), mkdir(outRoot); end
save(fullfile(outRoot, sprintf('BigT_%s.mat', ra_safename(opts.marker))), 'BigT','-v7.3');

% export 
writetable(BigT, fullfile(outRoot,[sprintf('BigT_%s', ra_safename(opts.marker)) '.csv']));

% Optional plotting
if opts.do_plot && exist('RAA_plot_behaviour_B1_to_B4','file') == 2
    try
        RAA_plot_behaviour_B1_to_B4(BigT, SessT, opts.marker, outRoot);
    catch ME
        warning('Plot failed for marker %s: %s', opts.marker, ME.message);
    end
end

end

% ======================================================================
function BigT = ra_join_session_metadata(BigT, SessT)
% Join SessT metadata into BigT by matching session path (or session name fallback)

if isempty(BigT) || height(BigT)==0, return; end
if ~ismember('SessionPath', SessT.Properties.VariableNames)
    error('SessT must contain SessionPath');
end

% Join key as string
BigT.key = string(BigT.session_path);
SessT.key = string(SessT.SessionPath);

% Select useful columns if present
keep = {'key'};
cand = {'Animal','SessionName','SoundPair','MotorCondition','DateNum'};
for i = 1:numel(cand)
    if ismember(cand{i}, SessT.Properties.VariableNames)
        keep{end+1} = cand{i}; %#ok<AGROW>
    end
end

S = SessT(:, keep);

% --- PATCH: enforce one row per session_path key (prevents doubling BigT) ---
[~, ia] = unique(string(S.key), 'stable');   % unique by key only
S = S(ia, :);

BigT = outerjoin(BigT, S, 'Keys','key', 'MergeKeys', true);

% SessionName fallback
if ~ismember('SessionName', BigT.Properties.VariableNames)
    BigT.SessionName = cellstr(string(BigT.session));
else
    % handle string/cell/categorical safely
    if all(cellfun(@isempty, ra_to_cellstr(BigT.SessionName)))
        BigT.SessionName = cellstr(string(BigT.session));
    end
end

% Create BlockSessionIndex within each (Animal, SoundPair), ordered by DateNum if present
if ismember('Animal', BigT.Properties.VariableNames) && ismember('SoundPair', BigT.Properties.VariableNames)
    BigT.BlockSessionIndex = ra_block_index(BigT);
end

% ---- plotter-friendly aliases (ONLY if missing) ----
BigT = ra_alias_table(BigT, 'Marker', {'Marker','marker'});
BigT = ra_alias_table(BigT, 'Subset', {'Subset','subset'});
BigT = ra_alias_table(BigT, 'Window', {'Window','window'});
BigT = ra_alias_table(BigT, 'BodyPart', {'BodyPart','bodypart'});

% Clean
BigT.key = [];

end

function idx = ra_block_index(T)
n = height(T);
idx = nan(n,1);

A = string(T.Animal);
P = string(T.SoundPair);

if ismember('DateNum', T.Properties.VariableNames)
    D = T.DateNum;
    if iscell(D) || isstring(D) || iscategorical(D)
        % try to coerce to numeric if possible
        Dn = nan(n,1);
        Ds = ra_to_cellstr(D);
        for i = 1:n
            di = str2double(Ds{i});
            if isfinite(di), Dn(i) = di; end
        end
        D = Dn;
    end
    [~,ord] = sortrows([A P D], [1 2 3]);
else
    ord = (1:n)';
end

u = unique([A P], 'rows');
for g = 1:size(u,1)
    m = (A==u(g,1) & P==u(g,2));
    ordg = ord(m(ord));
    idx(ordg) = (1:numel(ordg))';
end
end

% ======================================================================
% LOAD ONLY ONE MARKER (FAST)
% ======================================================================
function T = ra_load_session_metrics_marker(datapath, marker, use_mat)
% Load metrics for a given bodypart token (e.g. 'Pupil-EyeCam').
% IMPORTANT: do NOT load both MAT and CSV for the same metrics (duplicates).
% Preference order: MAT (full) -> CSV (fallback)

T = table();
anaDir = fullfile(datapath, 'analysis', 'behaviour');
if ~exist(anaDir,'dir'), return; end

pat_mat = ['*_' marker '_*_metrics*.mat']; % e.g. 20250722_n_Pupil-EyeCam_regular_..._metrics.mat
pat_csv = ['*_' marker '_*_metrics*.csv'];

% 1) Prefer MAT if allowed and present
if use_mat
    m = dir(fullfile(anaDir, '**', pat_mat));
    if ~isempty(m)
        for i = 1:numel(m)
            Ti = ra_load_one_metrics_mat(fullfile(m(i).folder, m(i).name));
            T  = ra_vertcat_align(T, Ti);
        end
        return
    end
end

% 2) Fallback to CSV (only if no MAT used/found)
c = dir(fullfile(anaDir, '**', pat_csv));
for i = 1:numel(c)
    Ti = ra_load_one_metrics_csv(fullfile(c(i).folder, c(i).name));
    T  = ra_vertcat_align(T, Ti);
end
end

function Ti = ra_load_one_metrics_mat(fpath)
Ti = table();
if ~exist(fpath,'file'), return; end
S = load(fpath);

if isfield(S,'T_full') && istable(S.T_full), Ti = S.T_full; return; end
if isfield(S,'T')      && istable(S.T),      Ti = S.T;      return; end
if isfield(S,'tbl')    && istable(S.tbl),    Ti = S.tbl;    return; end
if isfield(S,'S') && isstruct(S.S) && isfield(S.S,'table') && istable(S.S.table)
    Ti = S.S.table; return
end

fn = fieldnames(S);
for k = 1:numel(fn)
    x = S.(fn{k});
    if istable(x)
        Ti = x;
        return
    end
end
end

function Ti = ra_load_one_metrics_csv(fpath)
Ti = table();
if ~exist(fpath,'file'), return; end
Ti = readtable(fpath);
Ti = ra_reconstruct_roc_from_wide_columns(Ti);
end

function T = ra_reconstruct_roc_from_wide_columns(T)
v = T.Properties.VariableNames;

fpr_cols = v(startsWith(v,'rocFPR_'));
tpr_cols = v(startsWith(v,'rocTPR_'));

if ~isempty(fpr_cols) && ~any(strcmp(v,'rocFPR'))
    X = T{:, fpr_cols};
    roc = cell(height(T),1);
    for i = 1:height(T)
        xi = X(i,:);
        xi = xi(isfinite(xi));
        roc{i} = xi(:)';
    end
    T.rocFPR = roc;
    T = removevars(T, fpr_cols);
end

if ~isempty(tpr_cols) && ~any(strcmp(v,'rocTPR'))
    X = T{:, tpr_cols};
    roc = cell(height(T),1);
    for i = 1:height(T)
        xi = X(i,:);
        xi = xi(isfinite(xi));
        roc{i} = xi(:)'; 
    end
    T.rocTPR = roc;
    T = removevars(T, tpr_cols);
end

if any(strcmp(T.Properties.VariableNames,'rocFPR')) && ~iscell(T.rocFPR)
    T.rocFPR = ra_numeric_col_to_cellvec(T.rocFPR);
end
if any(strcmp(T.Properties.VariableNames,'rocTPR')) && ~iscell(T.rocTPR)
    T.rocTPR = ra_numeric_col_to_cellvec(T.rocTPR);
end
end

function c = ra_numeric_col_to_cellvec(x)
% Accept scalar/column numeric, returns {[]} for NaN, or {scalar} for finite
n = numel(x);
c = cell(n,1);
for i = 1:n
    xi = x(i);
    if ~isfinite(xi)
        c{i} = [];
    else
        c{i} = xi;
    end
end
end
% ======================================================================
% CONCAT WITH ALIGNMENT
% ======================================================================
function A = ra_vertcat_align(A, B)
if isempty(B) || height(B)==0, return; end
if isempty(A) || height(A)==0, A = B; return; end

va = A.Properties.VariableNames;
vb = B.Properties.VariableNames;
v_all = [va, setdiff(vb, va, 'stable')];

A = ra_add_missing_vars(A, v_all);
B = ra_add_missing_vars(B, v_all);

[A, B] = ra_make_types_compatible(A, B, v_all);

A = A(:, v_all);
B = B(:, v_all);
A = [A; B];
end

function T = ra_add_missing_vars(T, v_all)
vt = T.Properties.VariableNames;
missing = setdiff(v_all, vt, 'stable');
if isempty(missing), return; end

for k = 1:numel(missing)
    name = missing{k};
    fill = nan(height(T),1);

    if endsWith(name, {'_path','_name'}, 'IgnoreCase', true) || ...
       any(strcmpi(name, {'session','animal','pair','subset','window','bodypart','marker','SessionName','SoundPair','MotorCondition'}))
        fill = repmat({''}, height(T), 1);
    end

    if any(strcmpi(name, {'rocFPR','rocTPR'}))
        fill = repmat({[]}, height(T), 1);
    end

    T.(name) = fill;
end
end

function [A, B] = ra_make_types_compatible(A, B, v_all)
for k = 1:numel(v_all)
    name = v_all{k};
    if ~ismember(name, A.Properties.VariableNames) || ~ismember(name, B.Properties.VariableNames)
        continue
    end

    a = A.(name);
    b = B.(name);

    if strcmp(class(a), class(b))
        continue
    end

    if any(strcmpi(name, {'rocFPR','rocTPR'}))
        A.(name) = ra_to_cell_vector(a, height(A));
        B.(name) = ra_to_cell_vector(b, height(B));
        continue
    end

    A.(name) = ra_to_cell_any(a, height(A));
    B.(name) = ra_to_cell_any(b, height(B));
end
end

function c = ra_to_cell_vector(x, nrows)
if iscell(x), c = x; return; end

if isnumeric(x)
    if isvector(x)
        c = repmat({x(:)'}, nrows, 1);
    else
        if size(x,1) == nrows
            c = cell(nrows,1);
            for i = 1:nrows
                c{i} = x(i,:);
            end
        else
            c = repmat({x}, nrows, 1);
        end
    end
    return
end

c = repmat({[]}, nrows, 1);
end

function c = ra_to_cell_any(x, nrows)
if iscell(x), c = x; return; end
if iscategorical(x), c = cellstr(x(:)); return; end
if isstring(x), c = cellstr(x(:)); return; end
if ischar(x), c = repmat({x}, nrows, 1); return; end

if isnumeric(x) && isvector(x) && numel(x)==nrows
    c = num2cell(x(:));
    return
end

c = repmat({x}, nrows, 1);
end

function T = ra_alias_table(T, outName, candidates)
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

function s = ra_safename(s)
s = regexprep(s,'[^a-zA-Z0-9_]+','_');
end
