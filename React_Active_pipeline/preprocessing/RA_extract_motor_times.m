function [t_move_ts, t_arr_ts] = RA_extract_motor_times(datapath, trial_start_ts, trial_stop_ts)

n = numel(trial_start_ts);
t_move_ts = nan(n,1);
t_arr_ts  = nan(n,1);

dlc_file = fullfile(datapath,'video','DLC_data.mat');
if ~exist(dlc_file,'file')
    return
end

S = load(dlc_file, 'spout_likelihood', 'spout_center_mvt', 'time_face');

if ~isfield(S,'time_face') || ~isfield(S,'spout_likelihood')
    return
end

t = S.time_face(:);
dt = nanmedian(diff(t));
if dt < 1
    t_ts = t * 1e4;
else
    t_ts = t;
end

lk = S.spout_likelihood;
if isa(lk,'tsd'), lk = Data(lk); end
lk = double(lk(:));

mv = [];
if isfield(S,'spout_center_mvt')
    mv = S.spout_center_mvt;
    if isa(mv,'tsd'), mv = Data(mv); end
    mv = double(mv(:));
end

% parameters
lk_thr = 0.6;
mv_thr = 0.02;   % adjust if needed
hold_n = 3;

for tr = 1:n
    ii = find(t_ts>=trial_start_ts(tr) & t_ts<=trial_stop_ts(tr));
    if numel(ii) < 10, continue; end

    % arrival from likelihood
    t_arr_ts(tr) = first_sustained_crossing(t_ts(ii), lk(ii), lk_thr, hold_n);

    % movement start from spout movement if available
    if ~isempty(mv)
        t_move_ts(tr) = first_sustained_crossing(t_ts(ii), mv(ii), mv_thr, hold_n);
    end
end

end
