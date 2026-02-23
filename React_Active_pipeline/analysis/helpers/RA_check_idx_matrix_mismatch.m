function report = RA_check_idx_matrix_mismatch(sessions, opts)
% RA_check_idx_matrix_mismatch
% Quickly scan sessions for the "second bug": trial-index vectors (idxA/idxB)
% not matching the number of extracted rows in matA/matB for any marker/subset.
%
% Usage:
%   sessions = {... full paths ...};
%   report = RA_check_idx_matrix_mismatch(sessions);
%
% Optional:
%   opts.fs_video = 50;
%   opts.baseline_pre_s = 0.2;
%   opts.smoothing_win_s = 0.1;
%   opts.epoch_post_reward_s = 2.0;   % just for tvec sizing
%   opts.max_markers_per_session = Inf; % to speed up for debugging, e.g. 5
%   opts.only_subsets = {};            % {}, or {'nomotor'} etc
%
% Output:
%   report: struct array with fields:
%     .session, .subset, .group, .marker, .nIdx, .nOn, .nMat, .note

if nargin < 2, opts = struct(); end
if ~isfield(opts,'fs_video'), opts.fs_video = 50; end
if ~isfield(opts,'baseline_pre_s'), opts.baseline_pre_s = 0.20; end
if ~isfield(opts,'smoothing_win_s'), opts.smoothing_win_s = 0.10; end
if ~isfield(opts,'epoch_post_reward_s'), opts.epoch_post_reward_s = 2.0; end
if ~isfield(opts,'max_markers_per_session'), opts.max_markers_per_session = inf; end
if ~isfield(opts,'only_subsets'), opts.only_subsets = {}; end

fs = opts.fs_video;
baseline_samples = max(1, round(opts.baseline_pre_s * fs));
smoothSamples    = max(1, round(opts.smoothing_win_s * fs));

bodyparts = {
    'Pupil', {'pupil_area_004','pupil_center_004','pupil_center_004_mvt'};
    'Eye',   {'eye_area_004'};
    'Nostril',{'nostril_area','nostril_center','nostril_center_mvt'};
    'Nose',  {'nose_area','nose_center','nose_center_mvt'};
    'Cheek', {'cheek_center','cheek_center_mvt'};
    'Ear',   {'ear_center','ear_center_mvt'};
    'Jaw',   {'jaw_center','jaw_center_mvt'};
    'Tongue',{'tongue_center','tongue_center_mvt'};
    'Pupil-EyeCam', {'pupil_area_007','pupil_center_007','pupil_center_007_mvt'};
    'Eye-EyeCam',   {'eye_area_007'};
    'Spout', {'spout_likelihood'};

    % Physio (from physio_for_behaviour.mat)
    'OB-delta-power', {'OB_delta_rs','OB_delta_log10_rs','OB_delta_z_rs'};
    'OB-gamma-power', {'OB_gamma_rs','OB_gamma_log10_rs','OB_gamma_z_rs'};
    'OB-gamma-fast',  {'OB_gamma_fast_rs','OB_gamma_fast_log10_rs'};
    'OB-delta-fast',  {'OB_delta_fast_rs','OB_delta_fast_log10_rs'};

    'Respiration', {'RespPower_rs','RespPower_log10_rs','RespPhase_rs'};
    'EMG', {'EMGPower_rs','EMGPower_log10_rs','EMGPhase_rs'};
    'Accelerometer', {'AccPower_rs','AccPower_log10_rs'};
    'Heart', {'HeartPower_rs','HeartPower_log10_rs','HeartPhase_rs'};
};

report = struct('session',{},'subset',{},'group',{},'marker',{}, ...
    'nIdx',{},'nOn',{},'nMat',{},'note',{});

for si = 1:numel(sessions)
    datapath = sessions{si};
    [~,sessname] = fileparts(datapath);
    fprintf('[%d/%d] %s\n', si, numel(sessions), sessname);

    try
        Baphy = load_baphy(datapath);
    catch
        addrep(sessname,'__load__','__load__','__load__',NaN,NaN,NaN,'FAILED: load_baphy');
        continue
    end

    dlc_file = fullfile(datapath,'video','DLC_data.mat');
    if ~exist(dlc_file,'file')
        addrep(sessname,'__load__','__load__','__load__',NaN,NaN,NaN,'FAILED: missing DLC_data.mat');
        continue
    end
    D = load(dlc_file);

    P = struct();
    physio_file = fullfile(datapath,'ephys','physio_for_behaviour.mat');
    if exist(physio_file,'file')
        P = load(physio_file);
    end

    n = Baphy.n_trials;
    trial_start_abs_s = Baphy.trial.abs_trialstart_s(:);
    stim_on_abs_s  = Baphy.trial.abs_stim_start_s(:);
    stim_off_abs_s = Baphy.trial.abs_stim_stop_s(:);

    arrival_abs_s = nan(n,1);
    if isfield(Baphy.trial,'spout_arrival_abs_s')
        arrival_abs_s = Baphy.trial.spout_arrival_abs_s(:);
    end

    stim_on_rel_s  = stim_on_abs_s  - trial_start_abs_s;
    stim_off_rel_s = stim_off_abs_s - trial_start_abs_s;
    arrival_rel_s  = arrival_abs_s  - trial_start_abs_s;

    arr_med_all = nanmedian(arrival_rel_s(isfinite(arrival_rel_s)));
    if isfinite(arr_med_all)
        arrival_rel_s(~isfinite(arrival_rel_s)) = arr_med_all;
    end

    reward_rel_s_pertrial = nan(n,1);
    if isfield(Baphy.trial,'abs_reward_s')
        reward_rel_s_pertrial = Baphy.trial.abs_reward_s(:) - trial_start_abs_s;
    elseif isfield(Baphy.trial,'abs_target_s')
        reward_rel_s_pertrial = Baphy.trial.abs_target_s(:) - trial_start_abs_s;
    end
    reward_rel_s = nanmedian(reward_rel_s_pertrial(isfinite(reward_rel_s_pertrial)));
    if ~isfinite(reward_rel_s)
        reward_rel_s = nanmedian(stim_off_rel_s(isfinite(stim_off_rel_s))) + 2.0;
    end
    trial_dur_s = reward_rel_s + opts.epoch_post_reward_s;
    n_samples = round(trial_dur_s*fs)+1;
    tvec = (0:n_samples-1)/fs; %#ok<NASGU>

    % trial masks (same as your script)
    mskGood = false(n,1);
    if isfield(Baphy,'idx') && isfield(Baphy.idx,'goodTrials')
        mskGood(Baphy.idx.goodTrials(:)) = true;
    else
        mskGood = isfinite(trial_start_abs_s);
    end

    mskTar = false(n,1); mskRef = false(n,1);
    if isfield(Baphy,'idx') && isfield(Baphy.idx,'Target')
        mskTar(Baphy.idx.Target(:)) = true;
        mskRef(Baphy.idx.Reference(:)) = true;
    elseif isfield(Baphy.trial,'type')
        mskTar(strcmpi(Baphy.trial.type,'Target')) = true;
        mskRef(strcmpi(Baphy.trial.type,'Reference')) = true;
    end

    mskNS = false(n,1);
    if isfield(Baphy,'idx') && isfield(Baphy.idx,'NoSound')
        mskNS(Baphy.idx.NoSound(:)) = true;
    elseif isfield(Baphy.trial,'is_nosound')
        mskNS = logical(Baphy.trial.is_nosound(:));
    end

    mskNM = false(n,1);
    if isfield(Baphy,'idx') && isfield(Baphy.idx,'NomotorExclusive')
        mskNM(Baphy.idx.NomotorExclusive(:)) = true;
    elseif isfield(Baphy.trial,'is_nomotor')
        mskNM = logical(Baphy.trial.is_nomotor(:));
    end

    mskReg = ~(mskNS | mskNM);

    subsets = struct();
    subsets(1).name = 'regular';
    subsets(1).idxA = find(mskTar & mskGood & mskReg);
    subsets(1).idxB = find(mskRef & mskGood & mskReg);

    subsets(2).name = 'nosound';
    subsets(2).idxA = find(mskTar & mskGood & mskNS);
    subsets(2).idxB = find(mskRef & mskGood & mskNS);

    subsets(3).name = 'nomotor';
    subsets(3).idxA = find(mskTar & mskGood & mskNM);
    subsets(3).idxB = find(mskRef & mskGood & mskNM);

    if ~isempty(opts.only_subsets)
        keep = false(1,numel(subsets));
        for k = 1:numel(subsets)
            keep(k) = any(strcmpi(subsets(k).name, opts.only_subsets));
        end
        subsets = subsets(keep);
    end

    % fixed arrival anchor as in your script
    allA = find(mskTar & mskGood);
    allB = find(mskRef & mskGood);
    medA = nanmedian(arrival_rel_s(allA));
    medB = nanmedian(arrival_rel_s(allB));
    arr_anchor = nan;
    if isfinite(medA) && isfinite(medB), arr_anchor = min(medA, medB);
    elseif isfinite(medA), arr_anchor = medA;
    elseif isfinite(medB), arr_anchor = medB;
    end %#ok<NASGU>

    % scan markers quickly
    mk_count = 0;
    for b = 1:size(bodyparts,1)
        markers = bodyparts{b,2};
        for m = 1:numel(markers)
            mk = markers{m};
            [t_s, y] = get_marker_xy(mk, D, P);
            if isempty(t_s)
                continue
            end

            mk_count = mk_count + 1;
            if mk_count > opts.max_markers_per_session
                break
            end

            for ss = 1:numel(subsets)
                subsetName = subsets(ss).name;

                idxA = subsets(ss).idxA;
                idxB = subsets(ss).idxB;

                idxA = idxA(isfinite(trial_start_abs_s(idxA)));
                idxB = idxB(isfinite(trial_start_abs_s(idxB)));

                onA = trial_start_abs_s(idxA);
                onB = trial_start_abs_s(idxB);

                % this second filter is the one that can desync in your current code
                onA2 = onA(isfinite(onA));
                onB2 = onB(isfinite(onB));

                % try extraction; if extract_matrix drops rows, catch it here
                try
                    matA = extract_matrix(t_s, y, onA2, (0:n_samples-1)'/fs, baseline_samples, smoothSamples);
                    matB = extract_matrix(t_s, y, onB2, (0:n_samples-1)'/fs, baseline_samples, smoothSamples);
                catch ME
                    addrep(sessname, subsetName, 'Tar/Ref', mk, numel(idxA)+numel(idxB), numel(onA2)+numel(onB2), NaN, ...
                        ['FAILED: extract_matrix: ' ME.identifier]);
                    continue
                end

                % report mismatches relevant to the bug
                if numel(idxA) ~= numel(onA2)
                    addrep(sessname, subsetName, 'Tar', mk, numel(idxA), numel(onA2), size(matA,1), 'idxA vs onA mismatch (post-filter)');
                end
                if size(matA,1) ~= numel(onA2)
                    addrep(sessname, subsetName, 'Tar', mk, numel(idxA), numel(onA2), size(matA,1), 'onA vs matA mismatch (extract drops rows?)');
                end
                if numel(idxB) ~= numel(onB2)
                    addrep(sessname, subsetName, 'Ref', mk, numel(idxB), numel(onB2), size(matB,1), 'idxB vs onB mismatch (post-filter)');
                end
                if size(matB,1) ~= numel(onB2)
                    addrep(sessname, subsetName, 'Ref', mk, numel(idxB), numel(onB2), size(matB,1), 'onB vs matB mismatch (extract drops rows?)');
                end
            end
        end
        if mk_count > opts.max_markers_per_session
            break
        end
    end
end

% print summary
if isempty(report)
    fprintf('No mismatches detected.\n');
else
    fprintf('\nMISMATCHES DETECTED: %d\n', numel(report));
    [uSess,~,ic] = unique({report.session});
    for k = 1:numel(uSess)
        fprintf('  %s : %d\n', uSess{k}, sum(ic==k));
    end
end


    function addrep(sess, subset, group, marker, nIdx, nOn, nMat, note)
        r.session = sess;
        r.subset  = subset;
        r.group   = group;
        r.marker  = marker;
        r.nIdx    = nIdx;
        r.nOn     = nOn;
        r.nMat    = nMat;
        r.note    = note;
        report(end+1) = r; %#ok<AGROW>
    end
end
