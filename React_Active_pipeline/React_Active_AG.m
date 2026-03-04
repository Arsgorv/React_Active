function React_Active_AG()
%{
This is a pipeline for React_Active project
Work in progres...

Jul 2025
Arsenii Goriachenkov
Paris, France
%}

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Settings: Select sessions
project = 'RA'; % 'RP' 'Tonotopy'
mode = 'experiment';     % 'training' or 'experiment'

% Form the list of sessions
selection = 14;

% Tvorozhok training
Dir{1} = PathForExperimentsReactActive({'Tvorozhok'}, 'training', 'none', '11_10');
Dir{2} = PathForExperimentsReactActive({'Tvorozhok'}, 'training', 'none', '13_12');
Dir{3} = PathForExperimentsReactActive({'Tvorozhok'}, 'training', 'none', 'all');

% Mochi training
Dir{4} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '11_10');
Dir{5} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '13_12');
Dir{6} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '14_15');
Dir{7} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '16_17');
Dir{8} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '18_19');
Dir{9} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '20_21');
Dir{10} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', 'all');

% Mix training
Dir{11} = PathForExperimentsReactActive({'Mochi', 'Tvorozhok'}, 'training', 'none', 'all');

% Tvorozhok experiment
Dir{12} = PathForExperimentsReactActive({'Tvorozhok'}, 'experiment', 'none', '14_15');
Dir{13} = PathForExperimentsReactActive({'Tvorozhok'}, 'experiment', 'none', '16_17');
Dir{14} = PathForExperimentsReactActive({'Tvorozhok'}, 'experiment', 'none', 'all');

% Mochi experiment
Dir{15} = PathForExperimentsReactActive({'Mochi'}, 'experiment', 'none', 'all');

% Mix experiment
Dir{16} = PathForExperimentsReactActive({'Mochi', 'Tvorozhok'}, 'experiment', 'none', 'all');


sessions = Dir{selection}.path';

% sessions with DLC
session_dlc = filter_sessions_with_dlc(sessions);

%% Settings:  Fixing sessions: Mochi
% remove sessions
remove_sess = {...
        'Z:\Arsenii\React_Active\training\Mochi\20251106_m',... % remove completely: no camera , no OB
        'Z:\Arsenii\React_Active\training\Mochi\20251112_m',... % remove completely: very small number of trialss
        'Z:\Arsenii\React_Active\training\Mochi\20251118_m',... % remove completely: short/no behav
        'Z:\Arsenii\React_Active\training\Mochi\20251118_n',... % remove completely: short/no behav
        'Z:\Arsenii\React_Active\training\Mochi\20251222_m',... % remove completely: no OB, no behaviour
        'Z:\Arsenii\React_Active\training\Mochi\20251104_m',... % fix? spout
        'Z:\Arsenii\React_Active\training\Mochi\20251111_n',... % fix? spout
        'Z:\Arsenii\React_Active\training\Mochi\20251211_n',... % consider separately? 
        'Z:\Arsenii\React_Active\training\Mochi\20251031_n',... % TEMPORARY REMOVE ; Weird OB for ob_events 
        'Z:\Arsenii\React_Active\training\Mochi\20251104_m',... % TEMPORARY REMOVE ; Weird OB for ob_events 
        'Z:\Arsenii\React_Active\training\Mochi\20251117_n',... % TEMPORARY REMOVE ; Weird OB for ob_events 
        'Z:\Arsenii\React_Active\training\Mochi\20251128_n',... % TEMPORARY REMOVE ; Weird OB for ob_events 
        'Z:\Arsenii\React_Active\training\Mochi\20251209_n',... % TEMPORARY REMOVE ; Weird OB for ob_events 
        'Z:\Arsenii\React_Active\training\Mochi\20251219_m',... % TEMPORARY REMOVE ; Weird OB for ob_events 
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250723_m',... % fix? spout
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250731_m',...% fix? spout + face/eye mismatch
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250801_m',...% fix? spout + eye mismatch
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250818_n',...% fix? spout + face mismatch. Although it looks fine, spout is indeed misdetected
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250808_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250813_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250814_m',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250814_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250815_m',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250815_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250816_m',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250818_m',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250818_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250820_m',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250820_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250821_m',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250821_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250822_m',... % TEMPORARY REMOVE ; Weird OB for ob_events
        'Z:\Arsenii\React_Active\training\Tvorozhok\20250822_n',... % TEMPORARY REMOVE ; Weird OB for ob_events
};

%% Settings:  Fixing sessions: Tvorozhok
%     'Z:\Arsenii\React_Active\training\Tvorozhok\20250712',... % it looks alright, but let's see later
%     'Z:\Arsenii\React_Active\training\Tvorozhok\20250716_n',... % it looks alright, but let's see later
%     'Z:\Arsenii\React_Active\training\Tvorozhok\20250718_n',... % face mismatch. I see no problems with processing, but keep an eye

%% Settings:  Final session selection
session_dlc = session_dlc((~ismember(session_dlc, remove_sess)));

keepIdx = ~ismember(sessions, remove_sess);
sessions = sessions(keepIdx);

%% Settings:  opts
opts = struct();
opts.mode = mode;
opts.project = project;

opts.do_sync = true;
opts.do_basic_preproc = true;
opts.do_motion_energy = false;
opts.do_physio_preproc = 1;
opts.do_ob_event_analysis = 1;

opts.smoothing_win = 1;
opts.block_gap_s = 600;      % split Conditioning vs PostTest (seconds)
opts.video_gap_s = 1.0;      % split camera trains (seconds)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% PreProcessing:  Sleep Scoring
% Ferret_ProcessData_BM % Original version
Master_SleepScoring_preproc(sessions)

%% PreProcessing:  Behaviour 
% Move DLC videos 
% to_dir = 'Z:\Arsenii\React_Active\experiment\Processed_data\Tvorozhok\'; 
to_dir = 'Z:\Arsenii\React_Active\training\Mochi\'; 

from_dir = 'E:\DLC';
copy_dlc_related_files(to_dir, from_dir)

% check if you have '_filtered.csv' everywhere
for sess = 1:numel(sessions)
    datapath = sessions{sess};
    disp('------------------------------------------')
    disp(['Working on ' datapath])
    missing = check_missing_DLC_filtered_csv(datapath);
end

Master_DLC_preproc(sessions, opts)

%% PreProcessing:  align all datastreams -> extract trial information -> reconstruct missed triggers -> make Epochs
Master_data_sync_preproc(sessions, false, project)

% Check data sync
for sess = 1:numel(sessions)
    datapath = sessions{sess};
    disp('------------------------------------------')
    disp(['Working on ' datapath])
    if contains(project, 'RA')
        QC = RA_sanity_check_datasync_active(datapath);
    elseif contains(project, 'RP')
        QC = RP_sanity_check_datasync_passive(datapath);
    end
end

%% PreProcessing:  FMA spikesorting
sessions = {'Z:\Arsenii\React_Active\experiment\Processed_data\Tvorozhok\20260213'};
opts = struct();
if contains(sessions{1}, 'React_Active')
    opts.phases = {'PreSleep', 'Conditioning', 'PostSleep', 'PostTest'};
elseif contains(sessions{1}, 'React_Passive_ephys')
    opts.phases = {'PreExp', 'Exp', 'PostExp'};
end
opts.force_sort = false;      % true if you want to re-run wave_clus
opts.do_analysis = true;
opts.save_unit_figs = true;  % set true only if you really want per-unit PNGs
opts.do_sanity_figs = true;    % will create the epoch-column sanity figure

opts.groups = {1:32};
opts.commonRef = {1:32};

Master_FMA_preproc(sessions, opts);

%% PreProcessing:  NP spikesorting

%% PreProcessing:  sync OneBox and OpenEphys
% LFP = load_channel_as_tsd(fullfile(datapath,'ephys','LFPDataNP',['LFPmat_' segName '.mat']), 10);

MAP = RAE_fit_onebox_to_master_map(datapath, segName, struct('onebox_adc_idx',3));
OUT = RAE_apply_map_to_kilosort(datapath, segName, fullfile(datapath,'ephys','NP_spikes',segName), ...
                                fullfile(datapath,'analysis','sync_onebox',['map_' segName '.mat']), struct('fs_ap',30000));
                            sessions = {'Z:\Arsenii\React_Active\experiment\Processed_data\Tvorozhok\20260213'};

opts.onebox_adc_idx = 3;    % ADC2
opts.fs_ap = 30000;
opts.force = false;
opts.phase_filter = {'PreSleep', 'Conditioning', 'PostSleep', 'PostTest'}; % optional
Master_NP_spikes_sync_preproc(sessions, opts);

%% PreProcessing:  ECG and Respi
% MakeHeartRateForSession_BM(varargin)
% [Times,Template,HeartRate,GoodEpoch]=DetectHeartBeats_EmbReact_BM(LFP,NoiseEpoch,Options,plo)
% MakeInstFreqForSession_BM

%% PreProcessing:  ripples
% [Spindles, meanVal, stdVal] = FindRipplesKJ(LFP, Epoch, varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Analysis:  Behaviour analysis
React_Active_Behaviour_AG(session_dlc);

%% Analysis:  Sound-evoked response analysis

%% Analysis:  Sound decoding

%% Analysis:  Reactivation analysis

%% Analysis:  Brain-Behaviour analysis

%% Analysis:  Control analysis

%% Analysis:  Dataset figures

%% Analysis:  Paper figures


end