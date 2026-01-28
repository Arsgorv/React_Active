function React_Active_Behaviour_AG(session_dlc, opts)
% This is a master script to process the DLC behavioural ferret data for
% ReactActive project
% Steps:
%   - Generation of the basic figures (behaviour_analysis)
% Under construction  - Study the correlation between OB/Cortical/Hippocampal gamma and pupil area (gamma_pupil_corr)
% Under construction  - Producing the composition video with all variables synced (composition_video_OB_DLC_ferret)

if nargin < 2, opts = struct(); end
if ~isfield(opts,'smoothing_win'), opts.smoothing_win = 0; end

%% optional exclusion upfront
remove = {
    'Z:\Arsenii\React_Active\training\Tvorozhok\20250712'
    'Z:\Arsenii\React_Active\training\Tvorozhok\20250716_n'
    'Z:\Arsenii\React_Active\training\Tvorozhok\20250718_n'
};
session_dlc = session_dlc(~ismember(session_dlc, remove));

%% SESS: per-trial behavioural analysis
for sess = 1:numel(session_dlc)
    datapath = session_dlc{sess};
    disp(['[behav_analysis] ' datapath])
    behaviour_analysis(datapath, opts.smoothing_win, 'TrialSet', 'all');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%     makePupilSummaryVideo(session_dlc{sess} , smoothing_win)
%     decoding_trial_identity(session_dlc{sess})


% DATASET: behavioural analysis
bodyparts = {
    'Pupil', {'pupil_area_004','pupil_center_004','pupil_center_004_mvt'};
%     'Pupil', {'pupil_center_004','pupil_center_004_mvt'};
    'Eye',   {'eye_area_004'};
    'Nostril',{'nostril_area','nostril_center','nostril_center_mvt'};
    'Nose',  {'nose_area','nose_center','nose_center_mvt'};
    'Cheek', {'cheek_center','cheek_center_mvt'};
    'Ear',   {'ear_center','ear_center_mvt'};
    'Jaw',   {'jaw_center','jaw_center_mvt'};
    'Tongue',{'tongue_center','tongue_center_mvt'};
    'Pupil_EyeCam', {'pupil_area_007','pupil_center_007','pupil_center_007_mvt'};
    'Eye_EyeCam',   {'eye_area_007'};
    'Spout', {'spout_likelihood'};
    'OB_gamma_power', {'SmoothGamma','SmoothGamma_rs', 'Phase_OB', 'Phase_OB_rs'};
    'Respiration', {'BreathingPower_rs', 'BreathingPower_log10_rs', 'Phase','Phase_rs'};
%     'Respiration_other', {'BreathingPower_var_rs', 'Breathing', 'Breathing_rs', 'Breathing_var_rs'};
    'EMG', {'EMGData_rs', 'EMGData_log10_rs', 'Phase_EMG', 'Phase_EMG_rs'};
    'Accelerometer', {'MovAcctsd_rs', 'MovAcctsd_log10_rs'};
    };

nbp = size(bodyparts,1);
% baseColors = lines(nbp);

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