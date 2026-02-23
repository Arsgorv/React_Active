function SessT = RAA_get_training_sessions(sessions_to_exclude)
% RAA_get_training_sessions
% Build session metadata table for React_Active training, Mochi + Tvorozhok,
% grouped by sound-pair blocks using PathForExperimentsReactActive.
%
% INPUT
%   sessions_to_exclude : cell array of full session paths to exclude
%                         (optional)
%
% OUTPUT
%   SessT : table with per-session metadata

if nargin < 1 || isempty(sessions_to_exclude)
    sessions_to_exclude = {};
end
if ~iscell(sessions_to_exclude)
    error('sessions_to_exclude must be a cell array of full session paths');
end

animals = {'Tvorozhok','Mochi'};

pairs.Tvorozhok = {'11_10','13_12'}; % training only
pairs.Mochi     = {'11_10','13_12','14_15','16_17','18_19'};

R = {};
for a = 1:numel(animals)
    an = animals{a};
    sp = pairs.(an);
    for b = 1:numel(sp)
        pair = sp{b};
        Dir = PathForExperimentsReactActive(an, 'training', 'none', pair);
        if ~isfield(Dir,'path') || isempty(Dir.path), continue; end

        for k = 1:numel(Dir.path)
            sess_path = Dir.path{k};
            sess_name = Dir.name{k};

            % Parse date from session name: first 8 digits in the name
            dt = parse_session_date(sess_name);

            R(end+1,:) = {an, pair, b, k, sess_name, sess_path, dt}; %#ok<AGROW>
        end
    end
end

SessT = cell2table(R, 'VariableNames', ...
    {'Animal','SoundPair','BlockIndex','BlockSessionIndex','SessionName','SessionPath','DateNum'});

% Exclude requested sessions
if ~isempty(sessions_to_exclude)
    mask_ex = false(height(SessT),1);
    for i = 1:numel(sessions_to_exclude)
        mask_ex = mask_ex | strcmp(SessT.SessionPath, sessions_to_exclude{i});
    end
    SessT(mask_ex,:) = [];
end

% Global ordering (chronological within animal)
SessT = sortrows(SessT, {'Animal','DateNum','SessionName'});

SessT.GlobalSessionIndex = nan(height(SessT),1);
uA = unique(SessT.Animal,'stable');
for i = 1:numel(uA)
    idx = strcmp(SessT.Animal, uA{i});
    SessT.GlobalSessionIndex(idx) = (1:sum(idx))';
end

% Motor condition (Mochi: operand/equal from 20251128_m onward)
SessT.MotorCondition = repmat("classic", height(SessT), 1);
for i = 1:height(SessT)
    if strcmp(SessT.Animal{i}, 'Mochi')
        % date threshold: 20251128
        if SessT.DateNum(i) >= datenum('2025-11-28','yyyy-mm-dd')
            SessT.MotorCondition(i) = "operand";
        end
    end
end

end

function dn = parse_session_date(sess_name)
tok = regexp(sess_name, '(\d{8})', 'tokens', 'once');
if isempty(tok)
    dn = NaN;
    return
end
s = tok{1};
dn = datenum(s, 'yyyymmdd');
end
