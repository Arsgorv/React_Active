function Baphy = RA_parse_baphy_active(datapath, trigBaphy, varargin)
%{
% RA_parse_baphy_active
%
% % Uses trigBaphy (from Master_data_sync_preproc; trigOE.baphy) as ground-truth
% absolute timing. exptevents StartTime/StopTime are TRIAL-RELATIVE and reset
% every trial. Absolute times are:
%   abs_*_s = abs_trialstart_s + rel_*.
%
% INPUT
%   datapath  : session folder
%   trigBaphy : struct (trigOE.baphy) with at least:
%               - n_trials
%               - t_raw_s OR trial_start_ts (seconds) for trial onsets
%               - (optional) trial_stop_ts (seconds) for trial ends
%
% NAME-VALUE (subset)
%   'KeepRaw'                (true)
%   'SaveMat'                (true)
%   'SaveFig'                (true)
%   'SmoothWin'              (5)
%   'MinRunLen'              (20)
%   'ForceEarlyToHit'        (true)
%   'ForceNoMotorTargetToNaN'(true)
%
%   'AddSpoutArrival'        (true)
%   'SpoutThr'               (0.6)
%   'SpoutMinDur'            (1*1e4)  % ticks
%   'SpoutMergeGap'          (1*1e4)  % ticks
%   'ArrLagNom'              (3.95)   % seconds
%   'ArrLagTol'              (0.5)    % seconds
%
% OUTPUT
%   Baphy.trial.rel_* : trial-relative seconds (from exptevents)
%   Baphy.trial.abs_* : absolute seconds (session clock, from trigBaphy)
%}

%% -------------------- defaults --------------------
KeepRaw = true;
SaveMat = true;
SaveFig = true;

MFile = '';
DLCSubdir = '';
OutTag = '';

SmoothWin = 5;
MinRunLen = 20;
ForceEarlyToHit = true;
ForceNoMotorTargetToNaN = true;

AddSpoutArrival = true;
SpoutThr = 0.6;
SpoutMinDur = 1*1e4;
SpoutMergeGap = 1*1e4;
ArrLagNom = 3.95;
ArrLagTol = 0.5;

TsRate = 1e4;

%% -------------------- parse varargin --------------------
if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        key = varargin{k};
        val = varargin{k+1};
        switch lower(key)
            case 'keepraw'
                KeepRaw = logical(val);
            case 'savemat'
                SaveMat = logical(val);
            case 'savefig'
                SaveFig = logical(val);

            case 'smoothwin'
                SmoothWin = val;
            case 'minrunlen'
                MinRunLen = val;
            case 'forceearlytohit'
                ForceEarlyToHit = logical(val);
            case 'forcenomotortargettonan'
                ForceNoMotorTargetToNaN = logical(val);

            case 'addspoutarrival'
                AddSpoutArrival = logical(val);
            case 'spoutthr'
                SpoutThr = val;
            case 'spoutmindur'
                SpoutMinDur = val;
            case 'spoutmergegap'
                SpoutMergeGap = val;
            case 'arrlagnom'
                ArrLagNom = val;
            case 'arrlagtol'
                ArrLagTol = val;
                
            case 'tsrate'
                TsRate = val;
            case 'mfile'
                MFile = val;
            case 'dlcsubdir'
                DLCSubdir = val;
            case 'outtag'
                OutTag = val;

            otherwise
                error('Unknown parameter: %s', key);
        end
    end
end

%% -------------------- basic checks --------------------
if nargin < 2 || isempty(trigBaphy) || ~isstruct(trigBaphy)
    error('trigBaphy must be provided (e.g., trigOE.baphy).');
end

if isfield(trigBaphy,'baphy') && isstruct(trigBaphy.baphy)
    TB = trigBaphy.baphy;
else
    TB = trigBaphy;
end

if ~isfield(TB,'n_trials') || isempty(TB.n_trials)
    error('n_trials missing in trigBaphy.');
end

stimdir = fullfile(datapath,'stim');
if ~exist(stimdir,'dir')
    error('stim folder not found: %s', stimdir);
end

%% ---------------- load baphy mfile ----------------
if isempty(MFile)
    d = dir(fullfile(stimdir,'*.m'));
    if isempty(d), error('No .m files found in: %s', stimdir); end
    if numel(d) > 1
        warning('RA_parse_baphy_active:MultipleMFiles', ...
            'Multiple .m files found in %s. Using latest only. Use Master_data_sync_preproc for stitched multi-mfile parsing.', stimdir);
    end
    [~,ix] = max([d.datenum]);
    mfile = fullfile(stimdir, d(ix).name);
else
    mfile = MFile;
end

clear globalparams exptparams exptevents
run(mfile);

if ~exist('exptevents','var') || ~exist('exptparams','var')
    error('exptevents/exptparams not found after running %s', mfile);
end

%% ---------------- canonical nTrials ----------------
% nTrials = double(TB.n_trials);
% if ~isfinite(nTrials) || nTrials<=0
%     error('Invalid trigBaphy.n_trials: %s', mat2str(TB.n_trials));
% end
% trial_id = (1:nTrials)';
% abs_trialstart_s = nan(nTrials,1);
% abs_trialend_s   = nan(nTrials,1);
% 
% % Prefer t_raw_s if it matches nTrials (common in your trigOE struct)
% if isfield(TB,'t_raw_s') && ~isempty(TB.t_raw_s)
%     t0 = TB.t_raw_s(:);
%     if numel(t0) >= nTrials
%         abs_trialstart_s = t0(1:nTrials);
%     end
% end
% 
% % Fallback: trial_start_ts (ticks) -> seconds
% if ~any(isfinite(abs_trialstart_s)) && isfield(TB,'trial_start_ts') && ~isempty(TB.trial_start_ts)
%     t0_ts  = TB.trial_start_ts(:);
%     if numel(t0_ts ) >= nTrials
%         abs_trialstart_s = t0_ts(1:nTrials) ./ TsRate;
%     end
% end
% 
% % Extra fallback: t_raw_ts (ticks) -> seconds
% if ~any(isfinite(abs_trialstart_s)) && isfield(TB,'t_raw_ts') && ~isempty(TB.t_raw_ts)
%     t0_ts = TB.t_raw_ts(:);
%     if numel(t0_ts) >= nTrials
%         abs_trialstart_s = t0_ts(1:nTrials) ./ TsRate;
%     end
% end
% 
% if ~any(isfinite(abs_trialstart_s))
%     error('Could not obtain abs trial starts from trigBaphy (expected t_raw_s or trial_start_ts).');
% end
% 
% % Trial ends: if provided, use them ONLY if they are not identical to starts
% if isfield(TB,'trial_stop_s') && ~isempty(TB.trial_stop_s)
%     t1 = TB.trial_stop_s(:);
%     if numel(t1) >= nTrials
%         abs_trialend_s = t1(1:nTrials);
%     end
% end
% 
% if isfield(TB,'trial_stop_ts') && ~isempty(TB.trial_stop_ts)
%     t1_ts = TB.trial_stop_ts(:);
%     if numel(t1_ts) >= nTrials
%         tmp_end_s = t1_ts(1:nTrials) ./ TsRate;
% 
%         % RA case: trial_stop_ts repeats trial_start_ts (one pulse per trial) -> ignore
%         if all(isfinite(abs_trialstart_s)) && all(isfinite(tmp_end_s))
%             if median(abs(tmp_end_s - abs_trialstart_s)) < 1e-3
%                 tmp_end_s(:) = NaN;
%             end
%         end
% 
%         if ~any(isfinite(abs_trialend_s))
%             abs_trialend_s = tmp_end_s;
%         end
%     end
% end

trAll = [exptevents.Trial];
trAll = trAll(trAll > 0 & isfinite(trAll));
if isempty(trAll)
    error('exptevents.Trial is empty/invalid; cannot define number of trials.');
end
nTrials = double(max(trAll));   % ground truth
trial_id = (1:nTrials)';

% optional: log mismatch vs TTL
nTrialsTTL = NaN;
if isfield(trigBaphy,'baphy') && isfield(trigBaphy.baphy,'n_trials') && ~isempty(trigBaphy.baphy.n_trials)
    nTrialsTTL = double(trigBaphy.baphy.n_trials);
end
if isfinite(nTrialsTTL) && nTrialsTTL ~= nTrials
    warning('Trial count mismatch: exptevents=%d, trigBaphyTTL=%d. Restricting to exptevents.', nTrials, nTrialsTTL);
end

abs_trialstart_s = nan(nTrials,1);
abs_trialend_s   = nan(nTrials,1);

% Prefer t_raw_s, else trial_start_ts (ticks->s)
t0 = [];
if isfield(trigBaphy.baphy,'t_raw_s') && ~isempty(trigBaphy.baphy.t_raw_s)
    t0 = trigBaphy.baphy.t_raw_s(:);
elseif isfield(trigBaphy.baphy,'trial_start_ts') && ~isempty(trigBaphy.baphy.trial_start_ts)
    t0 = trigBaphy.baphy.trial_start_ts(:) ./ 1e4;
end

if ~isempty(t0)
    useN = min(numel(t0), nTrials);
    abs_trialstart_s(1:useN) = t0(1:useN);
end

%% ---------------- parse exptevents (REL times, trial clock) ----------------
trial_type = repmat({'Unknown'}, nTrials, 1);
sound_name = repmat({''}, nTrials, 1);

rel_trialstart = nan(nTrials,1);     % TRIALSTART StartTime (usually 0)
rel_trialend   = nan(nTrials,1);     % TRIALSTOP/END time if available
rel_prestim_start = nan(nTrials,1);
rel_stim_start = nan(nTrials,1);
rel_stim_stop  = nan(nTrials,1);
rel_reward     = nan(nTrials,1);

prestim_silence  = nan(nTrials,1);
poststim_silence = nan(nTrials,1);
is_nosound = false(nTrials,1);

nEvents = numel(exptevents);
notesAll = cell(nEvents,1);
for i = 1:nEvents
    note = exptevents(i).Note;
    if iscell(note)
        if isempty(note), note = '';
        else, note = note{1};
        end
    end
    if isstring(note)
        if numel(note) > 1, note = note(1); end
        note = char(note);
    end
    if ~ischar(note), note = ''; end
    notesAll{i} = note;
end
trial_start_idx = find(contains(notesAll,'TRIALSTART'));
trial_seen = false(nTrials,1);

for ii = 1:numel(trial_start_idx)
    i0 = trial_start_idx(ii);
    tr0 = exptevents(i0).Trial;
    if ~isscalar(tr0) || ~isfinite(tr0) || tr0<1 || tr0>nTrials
        continue
    end
    tr = double(tr0);
    trial_seen(tr) = true;
    
    rel_trialstart(tr) = exptevents(i0).StartTime;
    
    k = i0 + 1;
    
    stimT = NaN; stimStopT = NaN; rewT = NaN;
    trType = ''; sName = '';
    preT = NaN; preVal = NaN; postVal = NaN;
    endT = NaN;
    
    while k <= nEvents && exptevents(k).Trial == tr
        trk = exptevents(k).Trial;
        if ~isscalar(trk) || ~isfinite(trk) || double(trk) ~= tr
            break
        end
        
        note = notesAll{k};
        
        if contains(note,'PreStimSilence')
            preT = exptevents(k).StartTime;
            
            if contains(note,'Reference')
                trType = 'Reference';
            elseif contains(note,'Target')
                trType = 'Target';
            end
            
            parts = strsplit(note, ',');
            parts = strtrim(parts);
            if numel(parts) >= 2
                cand = parts{2};
                if ~isempty(cand) && ~contains(cand,'Reference') && ~contains(cand,'Target')
                    sName = cand;
                end
            end
            
            tok = regexp(note, 'value\s*=?\s*([-+]?\d*\.?\d+)', 'tokens', 'once');
            if isempty(tok)
                tok = regexp(note, '([-+]?\d*\.?\d+)', 'tokens'); % fallback: any number
                if ~isempty(tok), tok = tok{end}; end             % take last number
            end
            if ~isempty(tok)
                preVal = str2double(tok{1});
            end
        end
        
        if isnan(stimT) && contains(note,'Stim') && contains(note,'snippet_sequence') && ~contains(note,'PreStim') && ~contains(note,'PostStim')
            stimT = exptevents(k).StartTime;
            stimStopT = exptevents(k).StopTime;
        elseif isnan(stimT) && contains(note,'Stim') && ~contains(note,'PreStim') && ~contains(note,'PostStim')
            stimT = exptevents(k).StartTime;
            stimStopT = exptevents(k).StopTime;
        end
        
        if isnan(rewT) && (contains(note,'AUTOMATIC REWARD') || contains(note,'LICK,EARLY') || contains(note,'LICK,HIT,RIGHT'))
            rewT = exptevents(k).StartTime;
        end
        
        if contains(note,'Stim') && isfield(exptevents,'Rove') && ~contains(note,'PreStim') && ~contains(note,'PostStim')
            rv = exptevents(k).Rove;
            if isscalar(rv) && isequal(rv,0)
                is_nosound(tr) = true;
            end
        end
        
        if contains(note,'PostStimSilence') && isnan(postVal)
            tok = regexp(note, 'value\s*=?\s*([-+]?\d*\.?\d+)', 'tokens', 'once');
            if isempty(tok)
                tok = regexp(note, '([-+]?\d*\.?\d+)', 'tokens'); % fallback: any number
                if ~isempty(tok), tok = tok{end}; end             % take last number
            end
            if ~isempty(tok)
                postVal = str2double(tok{1});
            end
        end
        
        if contains(note,'TRIALSTOP') || contains(note,'TRIALEND')
            endT = exptevents(k).StartTime;
        end
        
        stp = NaN;
        if isfield(exptevents,'StopTime')
            stp = exptevents(k).StopTime;
        end
        if isempty(stp)
            stp = NaN;
        elseif ~isscalar(stp)
            stp = max(stp(:));
        end
        if isfinite(stp)
            if isnan(endT) || stp > endT
                endT = stp;
            end
        end
        
        k = k + 1;
    end
    
    if ~isempty(trType), trial_type{tr} = trType; end
    if ~isempty(sName),  sound_name{tr} = sName; end
    
    if isfinite(preT),      rel_prestim_start(tr) = preT; end
    if isfinite(stimT),     rel_stim_start(tr) = stimT; end
    if isfinite(stimStopT), rel_stim_stop(tr)  = stimStopT; end
    if isfinite(rewT),      rel_reward(tr)     = rewT; end
    if isfinite(preVal),    prestim_silence(tr)  = preVal; end
    if isfinite(postVal),   poststim_silence(tr) = postVal; end
    if isfinite(endT),      rel_trialend(tr) = endT; end
end

idxRef = find(strcmpi(trial_type,'Reference'));
idxTar = find(strcmpi(trial_type,'Target'));

idxNoSound = find(is_nosound);
noSoundType = trial_type(idxNoSound);

%% ---------------- ProbeTrialLst (NoMotor) ----------------
idxNoMotor = [];
noMotorType = {};
idxNoMotorExclusive = [];
noMotorExclusiveType = {};

if isfield(exptparams,'ProbeTrialLst') && ~isempty(exptparams.ProbeTrialLst)
    s = exptparams.ProbeTrialLst;
    if isstring(s), s = char(s); end
    if ~ischar(s), s = ''; end
    nums = regexp(s, '[-+]?\d+(\.\d+)?', 'match');
    probeVec = str2double(nums);
    if ~isempty(probeVec)
        idxNoMotor = find(probeVec==1);
        idxNoMotor = idxNoMotor(:);
        idxNoMotor = idxNoMotor(idxNoMotor>=1 & idxNoMotor<=nTrials);

        noMotorType = trial_type(idxNoMotor);

        idxNoMotorExclusive = setdiff(idxNoMotor, idxNoSound);
        idxNoMotorExclusive = idxNoMotorExclusive(:);
        noMotorExclusiveType = trial_type(idxNoMotorExclusive);
    end
end

%% -------------------- ABS event times (seconds) --------------------
abs_stim_start_s = abs_trialstart_s + rel_stim_start;
abs_stim_stop_s  = abs_trialstart_s + rel_stim_stop;
abs_reward_s     = abs_trialstart_s + rel_reward;

% If trial ends are missing OR useless (e.g. abs_end == abs_start), use exptevents rel_trialend
dur0 = abs_trialend_s - abs_trialstart_s;
if ~any(isfinite(abs_trialend_s)) || mean(isfinite(abs_trialend_s)) < 0.5 || median(dur0(isfinite(dur0))) < 1e-3
    abs_trialend_s = abs_trialstart_s + rel_trialend;

    % Fill missing ends with next start (seconds domain)
    for tr = 1:nTrials-1
        if ~isfinite(abs_trialend_s(tr)) && isfinite(abs_trialstart_s(tr+1))
            abs_trialend_s(tr) = abs_trialstart_s(tr+1);
        end
    end
end

% last trial end fallback if still missing
if ~isfinite(abs_trialend_s(nTrials))
    relDur = rel_trialend - rel_trialstart;
    relDur = relDur(isfinite(relDur) & relDur>0);
    if ~isempty(relDur)
        abs_trialend_s(nTrials) = abs_trialstart_s(nTrials) + median(relDur);
    end
end

trial_dur_s = abs_trialend_s - abs_trialstart_s;
stim_dur_s  = abs_stim_stop_s - abs_stim_start_s;

%% -------------------- Performance --------------------
Perf = struct();
Perf.n_trials = nTrials;

Perf.lick_rate_all = nan(nTrials,1);
Perf.hit_rate_all  = nan(nTrials,1);
Perf.hit_all       = nan(nTrials,1);
Perf.this_trial    = repmat({''}, nTrials, 1);

if isfield(exptparams,'Performance') && ~isempty(exptparams.Performance)
    P = exptparams.Performance;
    nP = numel(P);

    % handle trailing summary
    useN = min(nTrials, nP);
    if nP == nTrials + 1
        useN = nTrials;
    else
        % common: last element has empty ThisTrial => summary
        if nP >= 2 && isfield(P(end),'ThisTrial') && isempty(P(end).ThisTrial)
            useN = min(nTrials, nP-1);
        end
    end

    for i = 1:useN
        if isfield(P(i),'LickRate') && ~isempty(P(i).LickRate)
            v = P(i).LickRate;
            if ~isscalar(v), v = v(1); end
            Perf.lick_rate_all(i) = double(v);
        end
        if isfield(P(i),'HitRate') && ~isempty(P(i).HitRate)
            v = P(i).HitRate;
            if ~isscalar(v), v = v(1); end
            Perf.hit_rate_all(i) = double(v);
        end
        if isfield(P(i),'Hit') && ~isempty(P(i).Hit)
            v = P(i).Hit;
            if ~isscalar(v), v = v(1); end
            Perf.hit_all(i) = double(v);
        end
        if isfield(P(i),'ThisTrial') && ~isempty(P(i).ThisTrial)
            tt = P(i).ThisTrial;
            if iscell(tt), tt = tt{1}; end
            if isstring(tt), tt = char(tt(1)); end
            if ~ischar(tt), tt = ''; end
            Perf.this_trial{i} = tt;
        end
    end
end

Perf.firstlick_ref = [];
Perf.firstlick_tar = [];
if isfield(exptparams,'FirstLick') && ~isempty(exptparams.FirstLick)
    if isfield(exptparams.FirstLick,'Ref'), Perf.firstlick_ref = exptparams.FirstLick.Ref(:); end
    if isfield(exptparams.FirstLick,'Tar'), Perf.firstlick_tar = exptparams.FirstLick.Tar(:); end
end

Perf.hit_target_raw = Perf.hit_all(idxTar);
Perf.hit_target = Perf.hit_target_raw;

if ForceNoMotorTargetToNaN && ~isempty(idxNoMotorExclusive)
    isNMTarget = strcmpi(noMotorExclusiveType,'Target');
    nm_abs_idx = idxNoMotorExclusive(isNMTarget);
    [inTar, nm_pos_in_tar] = ismember(nm_abs_idx, idxTar);
    nm_pos_in_tar = nm_pos_in_tar(inTar);
    if ~isempty(nm_pos_in_tar)
        Perf.hit_target(nm_pos_in_tar) = NaN;
    end
end

if ForceEarlyToHit && ~isempty(Perf.this_trial)
    earlyGlobIdx = find(strcmpi(Perf.this_trial,'Early'));
    [isEarlyInTar, locInTar] = ismember(earlyGlobIdx, idxTar);
    tarPos = locInTar(isEarlyInTar);
    if ~isempty(tarPos)
        Perf.hit_target(tarPos) = 1;
    end
end

%% -------------------- fatigue cutoff -> good trials --------------------
goodTrials = (1:nTrials)';
fatigue_cut_trial = NaN;

hr = Perf.hit_rate_all(:);
if any(isfinite(hr))
    hr2 = hr;
    if SmoothWin > 1
        hr2 = movmean(hr2, SmoothWin, 'omitnan');
    end

    dHr = diff(hr2);
    isDown = dHr < 0;

    switchRuns = [true; diff(isDown)~=0];
    runStart = find(switchRuns);
    runLen   = diff([runStart; numel(isDown)+1]);
    runVal   = isDown(runStart);

    ix2 = find(runVal & runLen >= MinRunLen, 1, 'first');
    if ~isempty(ix2)
        firstBad = runStart(ix2);
        goodTrials = (1:firstBad)';
        fatigue_cut_trial = firstBad;
    end
end

%% -------------------- spout arrival (optional) --------------------
spout_arrival_abs_s = nan(nTrials,1);
spout_arrival_rel_s = nan(nTrials,1);
spout_arrival_abs_s_raw = nan(nTrials,1);
spout_arrival_rel_s_raw = nan(nTrials,1);
spout_arrival_detected = false(nTrials,1);
spout_arrival_imputed = false(nTrials,1);
spout_arrival_usable = false(nTrials,1);

spoutQC = struct();
spoutQC.source = '';
spoutQC.used_affine_remap = false;
spoutQC.n_raw = 0;
spoutQC.n_imputed = 0;
spoutQC.n_missing = nTrials;
spoutQC.temp = [];
spoutQC.t_s = [];
spoutQC.arrLagMin = ArrLagNom - ArrLagTol;
spoutQC.arrLagMax = ArrLagNom + ArrLagTol;
spoutQC.thr_like = nan(nTrials,1);
spoutQC.thr_mvt = nan(nTrials,1);

if AddSpoutArrival
    if isempty(DLCSubdir)
        dlcfile = fullfile(datapath,'video','DLC_data.mat');
    else
        dlcfile = fullfile(datapath,'video',DLCSubdir,'DLC_data.mat');
    end

    if exist(dlcfile,'file')
        S = load(dlcfile, 'spout_likelihood', 'spout_center_mvt', 'time_face_source');

        if isfield(S,'spout_likelihood') && ~isempty(S.spout_likelihood)
            spout_likelihood = S.spout_likelihood;
            if isfield(S,'time_face_source') && ~isempty(S.time_face_source)
                spoutQC.source = S.time_face_source;
            else
                spoutQC.source = 'unknown';
            end

            like_raw = Data(spout_likelihood);
            like_raw = like_raw(:);
            t_sp = Range(spout_likelihood);
            t_sp = ensure_seconds_vector(t_sp);

            if strcmpi(spoutQC.source, 'csv_column') && ...
                    isfield(trigBaphy,'face_cam') && isfield(trigBaphy.face_cam,'t_raw_ts') && ...
                    ~isempty(trigBaphy.face_cam.t_raw_ts)
                [t_sp, didMap] = maybe_affine_map_to_face_ttl(t_sp, trigBaphy.face_cam.t_raw_ts(:) ./ 1e4);
                spoutQC.used_affine_remap = didMap;
            end

            like_smooth = runmean(like_raw, 20);
            like_smooth = fillmissing_linear_local(like_smooth);

            hasMvt = false;
            if isfield(S,'spout_center_mvt') && ~isempty(S.spout_center_mvt)
                spout_mvt = S.spout_center_mvt;
                mvt_raw = Data(spout_mvt);
                mvt_raw = mvt_raw(:);
                if numel(mvt_raw) == numel(t_sp)
                    mvt_smooth = runmean(mvt_raw, 10);
                    mvt_smooth = fillmissing_linear_local(mvt_smooth);
                    hasMvt = true;
                end
            end

            arrLagMin = ArrLagNom - ArrLagTol;
            arrLagMax = ArrLagNom + ArrLagTol;

            for tr = 1:nTrials
                if ~isfinite(abs_trialstart_s(tr))
                    continue
                end

                winBeg = abs_trialstart_s(tr) + arrLagMin;
                winEnd = abs_trialstart_s(tr) + arrLagMax;
                baseBeg = abs_trialstart_s(tr) + max(0, arrLagMin - 1.2);
                baseEnd = abs_trialstart_s(tr) + max(0, arrLagMin - 0.15);

                kWin = find(t_sp >= winBeg & t_sp <= winEnd);
                if isempty(kWin)
                    continue
                end

                kBase = find(t_sp >= baseBeg & t_sp <= baseEnd);
                if isempty(kBase)
                    kBase = find(t_sp < winBeg);
                    if numel(kBase) > 100
                        kBase = kBase(max(1, end-99):end);
                    end
                end

                xw = like_smooth(kWin);
                tw = t_sp(kWin);
                xb = like_smooth(kBase);

                thLike = choose_trial_threshold(xb, SpoutThr);
                spoutQC.thr_like(tr) = thLike;

                dt = nanmedian(diff(tw));
                if ~isfinite(dt) || dt <= 0
                    hold_n = 3;
                else
                    hold_n = max(3, round(0.04 / dt));
                end

                t0 = first_sustained_crossing(tw, xw, thLike, hold_n);

                if ~isfinite(t0) && hasMvt
                    mw = mvt_smooth(kWin);
                    mb = mvt_smooth(kBase);
                    thMvt = choose_trial_threshold(mb, nan);
                    spoutQC.thr_mvt(tr) = thMvt;
                    t0 = first_sustained_crossing(tw, mw, thMvt, hold_n);
                end

                if isfinite(t0)
                    spout_arrival_abs_s_raw(tr) = t0;
                    spout_arrival_rel_s_raw(tr) = t0 - abs_trialstart_s(tr);
                    spout_arrival_abs_s(tr) = t0;
                    spout_arrival_rel_s(tr) = t0 - abs_trialstart_s(tr);
                    spout_arrival_detected(tr) = true;
                    spout_arrival_usable(tr) = true;
                end
            end

            for tr = 1:nTrials
                if spout_arrival_detected(tr)
                    continue
                end
                relImp = local_median_arrival(spout_arrival_rel_s_raw, tr, 10);
                if isfinite(relImp) && relImp >= arrLagMin && relImp <= arrLagMax && isfinite(abs_trialstart_s(tr))
                    spout_arrival_rel_s(tr) = relImp;
                    spout_arrival_abs_s(tr) = abs_trialstart_s(tr) + relImp;
                    spout_arrival_imputed(tr) = true;
                    spout_arrival_usable(tr) = true;
                end
            end

            spoutQC.temp = like_smooth;
            spoutQC.t_s = t_sp(:);
            spoutQC.n_raw = nnz(spout_arrival_detected);
            spoutQC.n_imputed = nnz(spout_arrival_imputed);
            spoutQC.n_missing = nnz(~spout_arrival_usable);

            if spoutQC.n_raw == 0
                warning('No raw spout arrivals detected in %s.', dlcfile);
            elseif spoutQC.n_raw < max(10, round(0.2 * nTrials))
                warning('Sparse raw spout arrivals: %d/%d trials.', spoutQC.n_raw, nTrials);
            end

            ok = isfinite(spout_arrival_rel_s_raw) & isfinite(rel_stim_start);
            if nnz(ok) > 10
                medArr = median(spout_arrival_rel_s_raw(ok));
                if medArr < 2 || medArr > 8
                    warning('Spout arrival rel looks off (median=%.3fs). Check face/Baphy alignment.', medArr);
                end
            end

            refRel = spout_arrival_rel_s_raw(idxRef);
            tarRel = spout_arrival_rel_s_raw(idxTar);
            if any(isfinite(refRel)) && any(isfinite(tarRel))
                mRef = median(refRel(isfinite(refRel)));
                mTar = median(tarRel(isfinite(tarRel)));
                if isfinite(mRef) && isfinite(mTar) && (mRef >= mTar)
                    warning('Spout arrival sanity: median Ref (%.3fs) >= median Tar (%.3fs).', mRef, mTar);
                end
            end
        else
            warning('DLC_data.mat loaded but spout_likelihood missing/empty.');
        end
    else
        warning('DLC_data.mat not found, skipping spout arrival.');
    end
end

%% -------------------- indices + masks --------------------
idx = struct();
idx.All = (1:nTrials)';
idx.Reference = idxRef(:);
idx.Target = idxTar(:);

idx.NoSound = idxNoSound(:);
idx.NoSoundType = noSoundType(:);

idx.NoMotor = idxNoMotor(:);
idx.NoMotorType = noMotorType(:);

idx.NoMotorExclusive = idxNoMotorExclusive(:);
idx.NoMotorExclusiveType = noMotorExclusiveType(:);

idx.goodTrials = goodTrials(:);
idx.fatigue_cut_trial = fatigue_cut_trial;

M = struct();
M.mskTarget = false(nTrials,1); M.mskTarget(idx.Target) = true;
M.mskRef    = false(nTrials,1); M.mskRef(idx.Reference) = true;

M.mskNS  = false(nTrials,1); M.mskNS(idx.NoSound) = true;
M.mskNS_T = false(nTrials,1);
M.mskNS_R = false(nTrials,1);
if ~isempty(idx.NoSound)
    M.mskNS_T(idx.NoSound(strcmpi(idx.NoSoundType,'Target'))) = true;
    M.mskNS_R(idx.NoSound(strcmpi(idx.NoSoundType,'Reference'))) = true;
end

M.mskNM  = false(nTrials,1); M.mskNM(idx.NoMotorExclusive) = true;
M.mskNM_T = false(nTrials,1);
M.mskNM_R = false(nTrials,1);
if ~isempty(idx.NoMotorExclusive)
    M.mskNM_T(idx.NoMotorExclusive(strcmpi(idx.NoMotorExclusiveType,'Target'))) = true;
    M.mskNM_R(idx.NoMotorExclusive(strcmpi(idx.NoMotorExclusiveType,'Reference'))) = true;
end

M.mskReg = ~(M.mskNS | M.mskNM);
M.mskGood = false(nTrials,1); M.mskGood(idx.goodTrials) = true;

%% -------------------- pack output --------------------
Baphy = struct();

Baphy.meta = struct();
Baphy.meta.datapath = datapath;
Baphy.meta.mfile = mfile;
Baphy.meta.parse_time = datestr(now);
Baphy.meta.clock_rel = 'baphy_trial_seconds';
Baphy.meta.clock_abs = 'OE_seconds_from_trigBaphy';
Baphy.meta.TsRate = TsRate;
Baphy.meta.trigBaphy = trigBaphy;

Baphy.meta.assumptions = struct();
Baphy.meta.assumptions.ForceNoMotorTargetToNaN = ForceNoMotorTargetToNaN;
Baphy.meta.assumptions.ForceEarlyToHit = ForceEarlyToHit;
Baphy.meta.fatigue = struct('SmoothWin',SmoothWin,'MinRunLen',MinRunLen);
Baphy.meta.spout = spoutQC;

Baphy.n_trials = nTrials;
Baphy.trial_id = trial_id(:);

Baphy.trial = struct();
Baphy.trial.type = trial_type(:);
Baphy.trial.sound_name = sound_name(:);

Baphy.trial.rel_trialstart = rel_trialstart(:);
Baphy.trial.rel_trialend   = rel_trialend(:);
Baphy.trial.rel_prestim_start = rel_prestim_start(:);
Baphy.trial.rel_stim_start = rel_stim_start(:);
Baphy.trial.rel_stim_stop  = rel_stim_stop(:);
Baphy.trial.rel_reward     = rel_reward(:);

Baphy.trial.abs_trialstart_s = abs_trialstart_s(:);
Baphy.trial.abs_trialend_s   = abs_trialend_s(:);
Baphy.trial.abs_stim_start_s = abs_stim_start_s(:);
Baphy.trial.abs_stim_stop_s  = abs_stim_stop_s(:);
Baphy.trial.abs_reward_s     = abs_reward_s(:);

Baphy.trial.trial_dur_s = trial_dur_s(:);
Baphy.trial.stim_dur_s  = stim_dur_s(:);

Baphy.trial.prestim_silence  = prestim_silence(:);
Baphy.trial.poststim_silence = poststim_silence(:);

Baphy.trial.is_nosound = is_nosound(:);
Baphy.trial.is_nomotor = false(nTrials,1);
if ~isempty(idx.NoMotorExclusive)
    Baphy.trial.is_nomotor(idx.NoMotorExclusive) = true;
end

Baphy.trial.spout_arrival_abs_s = spout_arrival_abs_s(:);
Baphy.trial.spout_arrival_rel_s = spout_arrival_rel_s(:);
Baphy.trial.spout_arrival_abs_s_raw = spout_arrival_abs_s_raw(:);
Baphy.trial.spout_arrival_rel_s_raw = spout_arrival_rel_s_raw(:);
Baphy.trial.spout_arrival_detected = spout_arrival_detected(:);
Baphy.trial.spout_arrival_imputed = spout_arrival_imputed(:);
Baphy.trial.spout_arrival_usable = spout_arrival_usable(:);

Baphy.idx = idx;
Baphy.M = M;
Baphy.Perf = Perf;

Baphy.task = struct();
if isfield(exptparams,'TrialObject') && ~isempty(exptparams.TrialObject), Baphy.task.TrialObject = exptparams.TrialObject; end
if isfield(exptparams,'ProbeTrialLst'), Baphy.task.ProbeTrialLst = exptparams.ProbeTrialLst; end
if isfield(exptparams,'StimTag') && ~isempty(exptparams.StimTag), Baphy.task.StimTag = exptparams.StimTag; end
if isfield(exptparams,'Reference') && ~isempty(exptparams.Reference), Baphy.task.Reference = exptparams.Reference; end
if isfield(exptparams,'Target') && ~isempty(exptparams.Target), Baphy.task.Target = exptparams.Target; end

if KeepRaw
    Baphy.raw = struct();
    if exist('globalparams','var'), Baphy.raw.globalparams = globalparams; end
    Baphy.raw.exptparams = exptparams;
    Baphy.raw.exptevents = exptevents;
end

%% -------------------- sanity checks --------------------
if any(isfinite(abs_trialstart_s)) && any(isfinite(abs_trialend_s))
    if median(abs_trialstart_s(isfinite(abs_trialstart_s))) > 1e6
        warning('abs_trialstart_s median is very large; looks like ticks not seconds. Check TsRate or inputs.');
    end
    if median(abs_trialend_s(isfinite(abs_trialend_s))) > 1e6
        warning('abs_trialend_s median is very large; looks like ticks not seconds. Check TsRate or inputs.');
    end
end

bad = find(isfinite(stim_dur_s) & stim_dur_s <= 0);
if ~isempty(bad)
    warning('Non-positive stim duration for trials: %s', mat2str(bad(1:min(10,end))'));
end

bad = find(isfinite(trial_dur_s) & trial_dur_s <= 0);
if ~isempty(bad)
    warning('Non-positive trial duration for trials: %s', mat2str(bad(1:min(10,end))'));
end

fracStimMissing = mean(~isfinite(rel_stim_start));
if fracStimMissing > 0.05
    warning('Stim start missing (rel) for %.1f%% of trials', 100*fracStimMissing);
end

missing = find(~trial_seen);
if ~isempty(missing)
    warning('exptevents defines nTrials=%d but TRIALSTART is missing for some trial ids: %s', ...
        nTrials, mat2str(missing(1:min(20,end))'));
end

trAll = [exptevents.Trial];
trAll = trAll(trAll>0);
mx = max(trAll);
if mx > nTrials
    warning('exptevents has trials up to %d, but n_trials=%d. Ignoring >nTrials.', mx, nTrials);
end

f1 = [];
fSpout = [];

%% -------------------- figure (legacy-style hit rate vs ABS time) --------------------
if SaveFig
    tAll = Baphy.trial.abs_trialstart_s(:);
    hr = Baphy.Perf.hit_rate_all(:);

    mskTarget = false(nTrials,1); mskTarget(idx.Target) = true;
    mskRef    = false(nTrials,1); mskRef(idx.Reference) = true;
    mskGood   = false(nTrials,1); mskGood(idx.goodTrials) = true;

    mskNoSound = false(nTrials,1); mskNoSound(idx.NoSound) = true;
    isNSTarget = strcmpi(idx.NoSoundType,'Target');
    isNSRef    = strcmpi(idx.NoSoundType,'Reference');
    mskNoSound_Target = false(nTrials,1);
    mskNoSound_Reference = false(nTrials,1);
    if ~isempty(idx.NoSound)
        mskNoSound_Target(idx.NoSound(isNSTarget)) = true;
        mskNoSound_Reference(idx.NoSound(isNSRef)) = true;
    end

    if isfield(exptparams,'ProbeTrialLst') && ~isempty(exptparams.ProbeTrialLst)
        idxNoMotorX = idx.NoMotorExclusive(:);
        nmType = idx.NoMotorExclusiveType;
        isNMTarget = strcmpi(nmType,'Target');
        isNMRef    = strcmpi(nmType,'Reference');

        mskNoMotor = false(nTrials,1); mskNoMotor(idxNoMotorX) = true;
        mskRegular = ~(mskNoSound | mskNoMotor);

        mskNoMotor_Target = false(nTrials,1);
        mskNoMotor_Reference = false(nTrials,1);
        if ~isempty(idxNoMotorX)
            mskNoMotor_Target(idxNoMotorX(isNMTarget)) = true;
            mskNoMotor_Reference(idxNoMotorX(isNMRef)) = true;
        end
    else
        mskRegular = ~mskNoSound;
        mskNoMotor_Target = false(nTrials,1);
        mskNoMotor_Reference = false(nTrials,1);
    end

    clr.Target = [0.98 0.50 0.45];
    clr.Reference = [0.00 0.39 0.00];
    shp.Regular = '.';
    shp.NoSound = 'x';
    shp.NoMotor = '*';
    sz.Good  = 10;
    sz.Other = 5;

    f1 = figure('Color','w','visible','on'); hold on; ylim([0 1]);
    set(f1, 'Units', 'Pixels', 'Position', [1 1 1920 1080/3]);

    plot(tAll(mskRegular & mskTarget &  mskGood), hr(mskRegular & mskTarget &  mskGood), shp.Regular, 'Color',clr.Target,    'MarkerSize',sz.Good,  'LineStyle','none', 'DisplayName','Target – good – regular');
    plot(tAll(mskRegular & mskTarget & ~mskGood), hr(mskRegular & mskTarget & ~mskGood), shp.Regular, 'Color',clr.Target,    'MarkerSize',sz.Other, 'LineStyle','none', 'DisplayName','Target – other – regular');
    plot(tAll(mskRegular & mskRef    &  mskGood), hr(mskRegular & mskRef    &  mskGood), shp.Regular, 'Color',clr.Reference, 'MarkerSize',sz.Good,  'LineStyle','none', 'DisplayName','Reference – good – regular');
    plot(tAll(mskRegular & mskRef    & ~mskGood), hr(mskRegular & mskRef    & ~mskGood), shp.Regular, 'Color',clr.Reference, 'MarkerSize',sz.Other, 'LineStyle','none', 'DisplayName','Reference – other – regular');

    plot(tAll(mskNoSound_Target    &  mskGood), hr(mskNoSound_Target    &  mskGood), shp.NoSound, 'Color',clr.Target,    'MarkerSize',sz.Good,  'LineStyle','none', 'DisplayName','Target – good – no-sound');
    plot(tAll(mskNoSound_Target    & ~mskGood), hr(mskNoSound_Target    & ~mskGood), shp.NoSound, 'Color',clr.Target,    'MarkerSize',sz.Other, 'LineStyle','none', 'DisplayName','Target – other – no-sound');
    plot(tAll(mskNoSound_Reference &  mskGood), hr(mskNoSound_Reference &  mskGood), shp.NoSound, 'Color',clr.Reference, 'MarkerSize',sz.Good,  'LineStyle','none', 'DisplayName','Reference – good – no-sound');
    plot(tAll(mskNoSound_Reference & ~mskGood), hr(mskNoSound_Reference & ~mskGood), shp.NoSound, 'Color',clr.Reference, 'MarkerSize',sz.Other, 'LineStyle','none', 'DisplayName','Reference – other – no-sound');

    if any(mskNoMotor_Target) || any(mskNoMotor_Reference)
        plot(tAll(mskNoMotor_Target    &  mskGood), hr(mskNoMotor_Target    &  mskGood), shp.NoMotor, 'Color',clr.Target,    'MarkerSize',sz.Good,  'LineStyle','none', 'DisplayName','Target – good – no-motor');
        plot(tAll(mskNoMotor_Target    & ~mskGood), hr(mskNoMotor_Target    & ~mskGood), shp.NoMotor, 'Color',clr.Target,    'MarkerSize',sz.Other, 'LineStyle','none', 'DisplayName','Target – other – no-motor');
        plot(tAll(mskNoMotor_Reference &  mskGood), hr(mskNoMotor_Reference &  mskGood), shp.NoMotor, 'Color',clr.Reference, 'MarkerSize',sz.Good,  'LineStyle','none', 'DisplayName','Reference – good – no-motor');
        plot(tAll(mskNoMotor_Reference & ~mskGood), hr(mskNoMotor_Reference & ~mskGood), shp.NoMotor, 'Color',clr.Reference, 'MarkerSize',sz.Other, 'LineStyle','none', 'DisplayName','Reference – other – no-motor');
    end

    ylabel('Hit rate');
    xlabel('Time from session start (s)');
    [~,session_name] = fileparts(datapath);
    sgtitle(['Session ' strrep(session_name,'_','\_')]);
    legend('Location','eastoutside', 'Interpreter','none');
    grid on; box on;
end

if SaveFig && AddSpoutArrival && ~isempty(spoutQC.t_s)
    try
        fSpout = figure('Color','w','visible','off');
        set(fSpout, 'Units', 'Pixels', 'Position', [50 50 1400 850]);

        subplot(3,1,1);
        plot(spoutQC.t_s, spoutQC.temp, 'k-'); hold on;
        yline(SpoutThr, 'r--');
        if any(spout_arrival_detected)
            plot(spout_arrival_abs_s_raw(spout_arrival_detected), ...
                interp1_unique_safe(spoutQC.t_s, spoutQC.temp, spout_arrival_abs_s_raw(spout_arrival_detected)), ...
                'go', 'MarkerSize', 4, 'LineWidth', 1);
        end
        if any(spout_arrival_imputed)
            plot(spout_arrival_abs_s(spout_arrival_imputed), ...
                interp1_unique_safe(spoutQC.t_s, spoutQC.temp, spout_arrival_abs_s(spout_arrival_imputed)), ...
                'rx', 'MarkerSize', 4, 'LineWidth', 1);
        end
        xlabel('Time (s)');
        ylabel('Spout likelihood');
        title('Spout likelihood and arrivals');
        grid on;

        subplot(3,1,2);
        plot(1:nTrials, spout_arrival_rel_s_raw, 'go'); hold on;
        plot(1:nTrials, spout_arrival_rel_s, 'k.');
        plot(find(spout_arrival_imputed), spout_arrival_rel_s(spout_arrival_imputed), 'rx');
        xlabel('Trial');
        ylabel('Arrival rel. time (s)');
        legend({'raw','final','imputed'}, 'Location', 'best');
        title(sprintf('raw=%d, imputed=%d, missing=%d', ...
            spoutQC.n_raw, spoutQC.n_imputed, spoutQC.n_missing));
        grid on;

        subplot(3,1,3);
        histogram(spout_arrival_rel_s_raw(isfinite(spout_arrival_rel_s_raw)), 30); hold on;
        histogram(spout_arrival_rel_s(spout_arrival_imputed), 30);
        xlabel('Arrival rel. time (s)');
        ylabel('Count');
        legend({'raw','imputed'}, 'Location', 'best');
        title('Arrival distribution');
        grid on;
    catch MEsp
        warning('Spout QC plot failed: %s', MEsp.message);
        if exist('fSpout','var') && ~isempty(fSpout)
            try, close(fSpout); end %#ok<TRYNC>
        end
        fSpout = [];
    end
end

%% -------------------- save --------------------
tag = OutTag;
if isempty(tag), tag = ''; else, tag = ['_' tag]; end

if SaveMat
    out = fullfile(stimdir, ['Baphy_RA' tag '.mat']);
    save(out,'Baphy','-v7.3');
end
if SaveFig && ~isempty(f1)
    try
        saveas(f1, fullfile(stimdir, ['session_structure' tag '.png']));
    catch
    end
end
if SaveFig && ~isempty(fSpout)
    try
        saveas(fSpout, fullfile(stimdir, ['spout_QC' tag '.png']));
    catch
    end
end

end

function yq = interp1_unique_safe(x, y, xq)
x = double(x(:));
y = double(y(:));
xq = double(xq(:));
yq = nan(size(xq));
if isempty(x) || isempty(y) || isempty(xq)
    return
end
n = min(numel(x), numel(y));
x = x(1:n);
y = y(1:n);
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
if numel(x) < 2
    return
end
[xu, ia] = unique(x, 'stable');
yu = y(ia);
if numel(xu) < 2
    return
end
yq = interp1(xu, yu, xq, 'linear', 'extrap');
end

function t_s = ensure_seconds_vector(t)
t_s = double(t(:));
if isempty(t_s), return; end
if nanmedian(diff(t_s)) > 1
    t_s = t_s ./ 1e4;
end
end

function [t_out_s, didMap] = maybe_affine_map_to_face_ttl(t_in_s, face_ttl_s)
t_out_s = t_in_s(:);
didMap = false;

m = min(numel(t_out_s), numel(face_ttl_s));
if m < 2
    return
end

x = double(t_out_s(1:m));
y = double(face_ttl_s(1:m));
[x, iu] = unique(x, 'stable');
y = y(iu);

if numel(x) >= 2
    p = polyfit(x, y, 1);
    t_out_s = polyval(p, double(t_out_s));
    didMap = true;
end
end

function x = fillmissing_linear_local(x)
x = double(x(:));
if isempty(x), return; end
good = find(isfinite(x));
if isempty(good)
    return
end
if numel(good) == 1
    x(~isfinite(x)) = x(good);
    return
end
bad = find(~isfinite(x));
if isempty(bad)
    return
end
x(bad) = interp1(good, x(good), bad, 'linear', 'extrap');
end

function th = choose_trial_threshold(xb, defaultThr)
xb = double(xb(:));
xb = xb(isfinite(xb));

if isempty(xb)
    if nargin < 2 || ~isfinite(defaultThr)
        th = 0;
    else
        th = defaultThr;
    end
    return
end

mu = median(xb);
sig = median(abs(xb - mu));
sig = 1.4826 * sig;

if ~isfinite(sig) || sig == 0
    sig = std(xb);
end
if ~isfinite(sig)
    sig = 0;
end

th = mu + 4 * sig;
if nargin >= 2 && isfinite(defaultThr)
    th = max(th, defaultThr);
end
end

function relImp = local_median_arrival(relRaw, tr, halfWidth)
relImp = nan;
ix = max(1, tr-halfWidth):min(numel(relRaw), tr+halfWidth);
ix(ix == tr) = [];
vals = relRaw(ix);
vals = vals(isfinite(vals));
if numel(vals) >= 3
    relImp = median(vals);
end
end
