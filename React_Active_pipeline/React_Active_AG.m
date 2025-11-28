function React_Active_AG()
%{
This is a pipeline for React_Active project
Work in progres...

Jul 2025
Arsenii Goriachenkov
Paris, France
%}

%% 
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