function React_Active_AG()
%{
This is a pipeline for React_Active project
Work in progres...

Jul 2025
Arsenii Goriachenkov
Paris, France
%}

%% Select sessions
% Form the list of sessions
selection = 10;

Dir{2} = PathForExperimentsReactActive({'Tvorozhok'}, 'training', 'none', '11_10');
Dir{3} = PathForExperimentsReactActive({'Tvorozhok'}, 'training', 'none', '13_12');
% Dir{4} = PathForExperimentsReactActive({'Tvorozhok'}, 'training', 'none', 'all', 'perf_drop');
Dir{5} = PathForExperimentsReactActive({'Tvorozhok'}, 'experiment', 'none', '14_15');

Dir{6} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '11_10');
Dir{7} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '13_12');
Dir{8} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '14_15');
Dir{9} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', '16_17');
Dir{10} = PathForExperimentsReactActive({'Mochi'}, 'training', 'none', 'all');
% Dir{1} = MergePathForExperiment(Dir{6},Dir{7},Dir{8},Dir{9});

sessions = Dir{selection}.path';

%% Select sessions with DLC
k = 1;
session_dlc = {};

for c = 1:length(sessions)
    dlc_path = fullfile(sessions{c}, 'video'); 
    files = dir(fullfile(dlc_path, '*_filtered.csv')); % Search for files ending with "_filtered.csv"
    
    if ~isempty(files) % Check if there are any matching files
        session_dlc{k} = sessions{c}; % Store the session path
        k = k + 1;
    else
        disp([Dir{selection}.path{c} ' - No DLC found']);
    end
end

session_dlc = session_dlc';


%% PreProcessing: Sleep Scoring
% Ferret_ProcessData_BM % Original version
Master_SleepScoring_preproc

%% PreProcessing: Behaviour 
React_Active_Behaviour_AG

%% PreProcessing: FMA spikesorting
run_wave_clus

%% PreProcessing: NP spikesorting


%% Sound-evoked response analysis


%% Reactivation analysis


%% Dataset figures


%% Paper figures


end