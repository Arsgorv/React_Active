function E = RA_build_epochs_active(datapath, trigOE, B)

TsRate = 1e4;

% --- checks ---
if ~isfield(B,'n_trials') || isempty(B.n_trials)
    error('RA_build_epochs_active:MissingNTrials','B.n_trials missing/empty.');
end
if ~isfield(B,'trial')
    error('RA_build_epochs_active:MissingTrialStruct','B.trial missing.');
end

req = {'abs_trialstart_s','abs_trialend_s','abs_stim_start_s','abs_stim_stop_s','spout_arrival_abs_s','type','is_nosound','is_nomotor'};
for i = 1:numel(req)
    if ~isfield(B.trial, req{i})
        error('RA_build_epochs_active:MissingField','Missing B.trial.%s', req{i});
    end
end

n = B.n_trials;

% --- absolute seconds -> timestamps (ticks) ---
trial_start_ts = B.trial.abs_trialstart_s(1:n) * TsRate;
trial_stop_ts  = B.trial.abs_trialend_s(1:n)   * TsRate;

stim_start_ts  = B.trial.abs_stim_start_s(1:n) * TsRate;
stim_stop_ts   = B.trial.abs_stim_stop_s(1:n)  * TsRate;

spout_arrival_ts = B.trial.spout_arrival_abs_s(1:n) * TsRate;

% --- validity masks ---
goodTrial = isfinite(trial_start_ts) & isfinite(trial_stop_ts) & (trial_stop_ts > trial_start_ts);
goodStim  = isfinite(stim_start_ts)  & isfinite(stim_stop_ts)  & (stim_stop_ts > stim_start_ts);
goodArr   = isfinite(spout_arrival_ts);

isRef = strcmpi(B.trial.type(1:n), 'Reference');
isTar = strcmpi(B.trial.type(1:n), 'Target');

isNoSound = logical(B.trial.is_nosound(1:n));

% nomotor:
% - B.trial.is_nomotor in your parser is already "exclusive" (NoMotor minus NoSound)
isNoMotorExclusive = logical(B.trial.is_nomotor(1:n));

% If you also want raw NoMotor (including overlap with NoSound), try to use idx.NoMotor when present
isNoMotorRaw = isNoMotorExclusive;
if isfield(B,'idx') && isfield(B.idx,'NoMotor') && ~isempty(B.idx.NoMotor)
    isNoMotorRaw = false(n,1);
    nm = B.idx.NoMotor(:);
    nm = nm(nm>=1 & nm<=n);
    isNoMotorRaw(nm) = true;
end

% --- pack output ---
E = struct();

E.trial_start_ts = trial_start_ts;
E.trial_stop_ts  = trial_stop_ts;
E.stim_start_ts  = stim_start_ts;
E.stim_stop_ts   = stim_stop_ts;
E.spout_arrival_ts = spout_arrival_ts;

% 1) trial_all
E.trial_all = intervalSet(trial_start_ts(goodTrial), trial_stop_ts(goodTrial));

% 2) stim_all
m = goodStim;
E.stim_all = intervalSet(stim_start_ts(m), stim_stop_ts(m));

% 3) stim_target
m = isTar & goodStim;
E.stim_target = intervalSet(stim_start_ts(m), stim_stop_ts(m));

% 4) stim_reference
m = isRef & goodStim;
E.stim_reference = intervalSet(stim_start_ts(m), stim_stop_ts(m));

% 5) trial_nosound
m = isNoSound & goodTrial;
E.trial_nosound = intervalSet(trial_start_ts(m), trial_stop_ts(m));

% 6) trial_nomotor (raw, may include overlap with NoSound if idx.NoMotor exists)
m = isNoMotorRaw & goodTrial;
E.trial_nomotor = intervalSet(trial_start_ts(m), trial_stop_ts(m));

% 7) trial_nomotorexclusive (NoMotor AND NOT NoSound)
m = isNoMotorExclusive & goodTrial;
E.trial_nomotorexclusive = intervalSet(trial_start_ts(m), trial_stop_ts(m));

% 8) [stim offset : spout arrival]
m = goodStim & goodArr & (spout_arrival_ts > stim_stop_ts);
E.stimoff_to_arrival = intervalSet(stim_stop_ts(m), spout_arrival_ts(m));

% 9) [stim onset : spout arrival]
m = goodStim & goodArr & (spout_arrival_ts > stim_start_ts);
E.stimon_to_arrival = intervalSet(stim_start_ts(m), spout_arrival_ts(m));

% Motor timing from DLC (ticks)
[E.motor_move_start_ts, E.motor_arrival_ts] = RA_extract_motor_times(datapath, trial_start_ts, trial_stop_ts);

% --- optional sanity check vs trigOE if provided ---
if nargin >= 2 && isstruct(trigOE) && isfield(trigOE,'baphy') && isfield(trigOE.baphy,'trial_start_ts') ...
        && ~isempty(trigOE.baphy.trial_start_ts)
    t_oe = trigOE.baphy.trial_start_ts(:);
    m = min(numel(t_oe), n);
    d = double(trial_start_ts(1:m)) - double(t_oe(1:m));
    if median(abs(d(isfinite(d)))) > 5  % >0.5 ms at 1e4 Hz
        warning('RA_build_epochs_active:StartMismatch', ...
            'B.abs_trialstart_s*1e4 differs from trigOE.baphy.trial_start_ts (median |diff|=%.1f ticks).', ...
            median(abs(d(isfinite(d)))));
    end
end

end
