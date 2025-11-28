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
selection = 2;

Dir{1} = PathForExperimentsReactActive({'Tvorozhok'}, 'head-fixed', 'none', '11_10');
Dir{2} = PathForExperimentsReactActive({'Tvorozhok'}, 'head-fixed', 'none', '13_12');
Dir{3} = PathForExperimentsReactActive({'Tvorozhok'}, 'head-fixed', 'none', 'all', 'perf_drop');

% Dir{3} = PathForExperiments_React_Active_AG({'Kosichka'}, 'head-fixed', 'none');
% Dir{4} = PathForExperiments_React_Active_AG({'Shropshire', 'Labneh', 'Brynza'}, 'head-fixed', 'none');
Dir{5} = MergePathForExperiment(Dir{1},Dir{2});

% Dirs_names = {'Tvorozhok', 'Kosichka', 'All animals', 'Arbitrary_Mix'};

sessions = Dir{selection}.path';

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
Ferret_ProcessData_BM

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