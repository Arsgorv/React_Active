function React_Active_Behaviour_AG(session_dlc, opts)
% React_Active_Behaviour_AG
% Master wrapper for per-session behaviour analysis (React_Active).
%
% Required session contents (per datapath):
%   video/DLC_data.mat
%   Master_sync.mat  (preferred) OR Baphy_RA.mat
%   ephys/SleepScoring_OBGamma.mat (optional, used if present)
%
% opts fields (all optional):
%   opts.smoothing_win_s   (default 0)
%   opts.trialSet          'all' | 'regular' | 'nosound' | 'nomotor' (default 'all')
%   opts.windows           cellstr subset of {'stimon_to_arrival','stimoff_to_arrival'} (default both)
%   opts.comparisons       cellstr subset of:
%       {'Tar_vs_Ref_regular','Tar_vs_Ref_NoSound','Tar_vs_Ref_NoMotor','Sound_vs_NoSound','Motor_vs_NoMotor'}
%   opts.save_fig          true/false (default true)
%   opts.visible_fig       true/false (default false)
%
% Example:
%   opts = struct(); opts.smoothing_win_s = 0.10;
%   opts.trialSet = 'all';
%   React_Active_Behaviour_AG(session_dlc, opts)
% Steps:
%   - Generation of the basic figures (behaviour_analysis)
% Under construction  - Study the correlation between OB/Cortical/Hippocampal gamma and pupil area (gamma_pupil_corr)
% Under construction  - Producing the composition video with all variables synced (composition_video_OB_DLC_ferret)

if nargin < 2, opts = struct(); end
opts = ra_default_opts(opts);

%% SESS: per-trial behavioural analysis
for sess = 1:numel(session_dlc)
    datapath = session_dlc{sess};
    disp(['[RA_behaviour] ' datapath])

    % 1) physio preprocessing (OB gamma resampled to video)
    if opts.do_physio_preproc
        try
            preproc_physio_for_behaviour(datapath, opts);
        catch ME
            warning('Physio preproc failed: %s', ME.message);
        end
    end

    % 2) behaviour analysis (figures + metrics)
    try
        behaviour_analysis(datapath, opts);
    catch ME
        warning('behaviour_analysis failed: %s', ME.message);
    end

    % 3) OB event interactions (spectrogram + phase locking)
    if opts.do_ob_event_analysis
        try
            ob_event_analysis(datapath, opts);
        catch ME
            warning('ob_event_analysis failed: %s', ME.message);
        end
    end
    
    close all
    
end

%% Check if all sessions are well processed
rootDir = 'Z:\Arsenii\React_Active\training\Tvorozhok';

opts = struct();
opts.ignore_tilde_sessions = true;
opts.save_csv = fullfile(rootDir, 'required_outputs_check_Tvorozhok.csv');

[Report, remove_sess] = RA_check_required_outputs(rootDir, opts);

bad = find(~[Report.ok]);
for i = 1:numel(bad)
    k = bad(i);
    fprintf('\nMISSING in %s\n', Report(k).session_path);
    for j = 1:numel(Report(k).missing_files)
        fprintf('  - %s\n', Report(k).missing_files{j});
    end
end
% tv_r = remove_sess;
% mo_r = remove_sess;
% remove_sess = [tv_r; mo_r];

%% Check idx trial mismatch
opts = struct();
opts.max_markers_per_session = 8;     
opts.only_subsets = {'nomotor','nosound'};  
rep = RA_check_idx_matrix_mismatch(sessions, opts);

if ~isempty(rep)
    T = struct2table(rep);
    disp(T);
end

%% DATASET: behavioural analysis
markers = {'Accelerometer', 'Cheeck', 'EMG', 'Eye-EyeCam', 'Jaw', 'Nose', 'Nostril', 'OB-delta-fast', 'OB-gamma-fast', 'Pupil-EyeCam', 'Respiration'};
for i = 1:numel(markers)
    opts.marker = markers{i};
    RAA_run_behaviour_across_sessions(remove_sess,opts)
end

path_csv = 'Z:\Arsenii\React_Active\training\Mochi\DS_figures\behaviour';
csv_files = {
    fullfile('Pupil_EyeCam', 'BigT_Pupil_EyeCam.csv'),...
    fullfile('Eye_EyeCam', 'BigT_Eye_EyeCam.csv'),...
    fullfile('Nostril', 'BigT_Nostril.csv'),...
    fullfile('Nose', 'BigT_Nose.csv'),...
    fullfile('Jaw', 'BigT_Jaw.csv'),...
    fullfile('EMG', 'BigT_EMG.csv'),...
    fullfile('Accelerometer', 'BigT_Accelerometer.csv'),...
    fullfile('OB_delta_fast', 'BigT_OB_delta_fast.csv'),...
    fullfile('OB_gamma_fast', 'BigT_OB_gamma_fast.csv'),...
};

RAA_collect_all_behaviour_tables(csv_files, fullfile(pwd,'across_sessions'));

%% DATASET: OB Events analysis
window = {'stimon_to_stimoff', 'stimoff_to_arrival', 'arr_to_stop'};
window_selection = 2;

% OB session summaries
RESob = RAA_collect_ob_trialmetrics_across_sessions(sessions, struct('alignName','stimOn'));
% RAA_plot_ob_across_sessions(RESob);

% Behaviour AUROC (pupil)
BRES = RAA_collect_behaviour_marker_AUROC_across_sessions(sessions, struct( ...
    'bodypart','Pupil-EyeCam', 'marker','pupil_area_007', ...
    'subset','regular', 'window',window{window_selection}));
% Scatter (OB theta vs pupil)
RAA_scatter_pupil_vs_ob_theta_AUROC(RESob, BRES, struct( ...
    'band','gamma', 'window',window{window_selection}, 'minN',10, 'useSpearman',false));
RAA_scatter_pupil_vs_ob_specificity(RESob, BRES, opts)
OUT = RAA_compare_emergence_speed_pupil_vs_ob(RESob, BRES, struct( ...
    'bandList',{{'delta','theta','gamma','highgamma'}}, ...
    'window',window{window_selection}, ...
    'thr',0.60, 'kConsec',2, 'minN',10, 'plot',true));

% Band ans resp
RESband = RAA_collect_ob_bandtraces_across_sessions(sessions, struct('alignName','stimOn'));
RAA_plot_ob_bandtrace_heatmaps(RESband);


RESresp = RAA_collect_resp_coupling_across_sessions(sessions, opts);
RAA_plot_resp_coupling_across_sessions(RESresp, opts)

%% %%%%%%%%%%%%%%%%%%%%%%%
%     makePupilSummaryVideo(session_dlc{sess} , smoothing_win)
%     decoding_trial_identity(session_dlc{sess})

for b = 1:nbp
%     bodypart_name = bodyparts{b,1};
    marker_names = bodyparts{b,2};
%     baseCol   = baseColors(b,:);
    n_markers = numel(marker_names);
    for m = 1:n_markers
        disp('Analysing behavioural data across sessions...')
%         behaviour_analysis_dataset(session_dlc, smoothing_win, AW.(marker_names{m}), bodyparts{b}, marker_names{m})
        behaviour_analysis_dataset(session_dlc, smoothing_win, bodyparts{b}, marker_names{m})
    end
end

%% Performance
% careful


% 
% % good 
% 'Z:\Arsenii\React_Active\training\Tvorozhok\20250723_m'
% 'Z:\Arsenii\React_Active\training\Tvorozhok\20250724_m'
% 'Z:\Arsenii\React_Active\training\Tvorozhok\20250724_n'

% remove
remove = {'Z:\Arsenii\React_Active\training\Tvorozhok\20250712' ,...
          'Z:\Arsenii\React_Active\training\Tvorozhok\20250716_n',...
          'Z:\Arsenii\React_Active\training\Tvorozhok\20250718_n'};
remove = remove';

keepIdx = ~ismember(session_dlc, remove);
session_dlc = session_dlc(keepIdx);

% SESS: Performance/pupil/OB gamma corr
sm_win_behav = 3; % in trials.
sm_win_ephys = 250; % in samples. fps = 50 -> 200/50 = 4s
disp('Correlating pergormance, pupil and OB gamma power ...')
for sess = 8:(length(session_dlc)-1)
    disp([num2str(length(session_dlc)-sess + 1) '/' num2str(length(session_dlc)) ' left'])
    disp(['Running session: ' session_dlc{sess}])
    pupil_gamma_performance_corr(session_dlc{sess}, sm_win_behav, sm_win_ephys)
    close all
end
% DATASET: 
signalType = 'area';%'pos';%
fig_pupil_gamma_performance_dataset(session_dlc,signalType)

%% Relevant for pupil data. 

% % SESS: Pupil physiology
% for sess = 1:length(session_dlc)
%     disp('Running physiology analysis...')    
%     pupil_physio(session_dlc{sess})
%     %     Ferret_Paper_CrossCorr_GammaPupil
% end
% % DATASET: Pupil physiology
% fig_pupil_physio_dataset(session_dlc);
% 
% % SESS: Study correlation between pupil and brain signal as a function of frequency band passed
% for sess = 10:length(session_dlc)
%     disp('Running pupil-brain correlation vs f.band analysis...')
%     stats{sess} = brain_pupil_fband(session_dlc{sess});
% end
% % DATASET: Correlate brain signals with pupil size
% fig_brain_pupil_fband_dataset(session_dlc)

% SESS: Correlate brain signals with pupil size
for sess = 1:length(session_dlc)
    disp('Running pupil-brain correlation analysis...')
    brain_pupil_corr(session_dlc{sess})
end
% DATASET: Correlate brain signals with pupil size
fig_brain_pupil_corr_dataset(session_dlc)

% Under construction: Generate composition video
composition_video_OB_DLC_ferret

end