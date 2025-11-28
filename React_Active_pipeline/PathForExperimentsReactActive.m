function Dir = PathForExperiments_React_Active_AG(ferret_name, setup_type, pharma_protocol, stimuli_protocol, additional_filetring)

% PathForExperiments_React_Active_AG - Retrieves session paths and metadata for OB experiments.
% You can select several ferrets and several stimuli_protocols at the same time

% Inputs:
%   ferret_name     - {'all', 'Shropshire', 'Brynza', 'Labneh', 'Edel', ...} or cell array of ferret names
%   setup_type      - {'all', 'freely-moving', 'head-fixed'} (default: 'all')
%   pharma_protocol - {'all', 'atropine', 'domitor', 'saline', 'fluoxetine', 'none'} (default: 'all')
%   stimuli_protocol- {'all', 'TCI', 'puretones', 'yves', 'contstream', 'TORCs', 'LSP', 'resting', 'none', ...} or cell array of protocols (default: 'all')
%
%%%%%%% Examples of use %%%%%%%
% Dir = PathForExperimentsOB('Shropshire', 'freely-moving', 'all', {'LSP', 'TORCs'}); -- Will pull out both LSP and TORCs freely-moving sessions for Shropshire
% Dir = PathForExperimentsOB({'Shropshire', 'Brynza'}, 'freely-moving', 'atropine'); -- Will pull out atropine freely-moving sessions for Shropshire and Brynza

%%
% Default parameter handling
if nargin < 1 || isempty(ferret_name)
    ferret_name = 'all';
end
if nargin < 2 || isempty(setup_type)
    setup_type = 'all';
end
if nargin < 3 || isempty(pharma_protocol)
    pharma_protocol = 'all';
end
if nargin < 4 || isempty(stimuli_protocol)
    stimuli_protocol = 'all';
end

% Ensure ferret_name and stimuli_protocol are cell arrays
if ~iscell(ferret_name)
    ferret_name = {ferret_name};
end
if ~iscell(stimuli_protocol)
    stimuli_protocol = {stimuli_protocol};
end

% Initialize Dir structure
Dir = struct();
Dir.path = {};
Dir.ExpeInfo = {};
Dir.name = {};

% Session data (add paths and session names)
data = {
    
    % Tvorozhok head-fixed
    'Tvorozhok', 'head-fixed', 'none', '11_10', 'Z:\Arsenii\React_Active\training\Tvorozhok\' {'20250712', '20250716_m', '20250716_n', '20250717_m', '20250718_n', '20250722_m', '20250722_n'};
    'Tvorozhok', 'head-fixed', 'none', '13_12', 'Z:\Arsenii\React_Active\training\Tvorozhok\' {'20250723_m', '20250724_m', '20250724_n', '20250725_m',    '20250725_n'    ,    '20250728_m'    , '20250729_m', '20250729_n',    '20250730_m'    , '20250730_n',...
                                                                                                    '20250731_m'    , '20250731_n', '20250801_m', '20250805_n', '20250806_m', '20250806_n', '20250808_n', '20250813_n', '20250814_m',...
                                                                                               '20250814_n', '20250815_m', '20250815_n', '20250816_m', '20250818_m', '20250818_n', '20250820_m', '20250820_n', '20250821_m', '20250821_n',...
                                                                                               '20250822_m', '20250822_n'};
% removed: '20250813_m'
                                                                                               %     'Tvorozhok', 'head-fixed', 'none', 'all', 'perf_drop' 'Z:\Arsenii\React_Active\training\Tvorozhok\' {'20250723_m', '20250724_m', '20250724_n', '20250725_m', '20250725_n', '20250728_m', '20250729_m', '20250729_n',  '20250730_m', '20250730_n',  '20250731_m',  '20250731_n',  '20250801_m', '20250805_n', '20250806_m', '20250806_n'};

    
    
    
    
    
    
    
    
    
    %     'Shropshire', 'head-fixed', 'none', 'puretones', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', {'20241220_puretones'};
    'Shropshire', 'head-fixed', 'none', 'TCI', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', {'20241123_TCI', '20241125_TCI', '20241126_TCI', '20241128_TCI', '20241129_TCI', '20241130_TCI'};
    'Shropshire', 'head-fixed', 'none', 'resting', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', { '20241210_resting', '20241211_resting', '20241212_resting', '20241213_resting', '20241214_resting'};
    'Shropshire', 'head-fixed', 'none', 'yves', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', {'20241120_yves_train', '20241123_yves_train', '20241125_yves_train', '20241126_yves_train', '20241128_yves_train', '20241129_yves_test', '20241130_yves_test', '20241203_yves_test'}; % 20241220_resting
    
    'Shropshire', 'head-fixed', 'none', 'contstream', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', {'20241209_contstream', '20241210_contstream', '20241211_contstream', '20241212_contstream', '20241213_contstream', '20241214_contstream'};
    'Shropshire', 'head-fixed', 'none', 'TORCs', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', {'20241204_TORCs', '20241205_TORCs', '20241206_TORCs', '20241209_TORCs', '20241210_TORCs', '20241211_TORCs', '20241212_TORCs', '20241213_TORCs', '20241214_TORCs'};
    'Shropshire', 'head-fixed', 'atropine', 'TORCs', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', {'20241220_TORCs_atropine', '20241222_TORCs_atropine', '20241225_TORCs_atropine',  '20241227_TORCs_atropine', '20241231_TORCs_atropine', '20250102_TORCs_atropine', '20250105_TORCs_atropine'};
    'Shropshire', 'head-fixed', 'saline', 'TORCs', '/media/nas8/OB_ferret_AG_BM/Shropshire/head-fixed', {'20241224_TORCs_saline', '20241226_TORCs_saline', '20241228_TORCs_saline', '20250101_TORCs_saline', '20250103_TORCs_saline',  '20250104_TORCs_saline', '20250108_TORCs_saline'};
    
    % Shropshire freely-moving
    'Shropshire', 'freely-moving', 'atropine', 'LSP', '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving', {'20241122_LSP', '20241123_LSP', '20241125_LSP', '20241126_LSP', '20241128_LSP', '20241130_LSP', '20241203_LSP'}; % 20241129_LSP
    'Shropshire', 'freely-moving', 'saline', 'LSP', '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving', { '20241218_LSP_saline',  '20241224_LSP_saline' '20241228_LSP_saline', '20250103_LSP_saline', '20250107_LSP_saline'}; % 20241230_LSP_saline 20241217_LSP_saline
    'Shropshire', 'freely-moving', 'none', 'TORCs', '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving', {'20241204_TORCs', '20241205_TORCs', '20241206_TORCs', '20241209_TORCs', '20241210_TORCs', '20241211_TORCs', '20241212_TORCs', '20241213_TORCs', '20241214_TORCs', '20241221_TORCs', '20241223_TORCs', '20241224_TORCs_short'};
    'Shropshire', 'freely-moving', 'none', 'contstream', '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving', {'20241120_contstream', '20241121_contstream', '20241209_contstream', '20241210_contstream', '20241211_contstream', '20241213_contstream', '20241214_contstream'};
    'Shropshire', 'freely-moving', 'none', 'resting', '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving', {'20241210_resting', '20241211_resting', '20241214_resting'};
    'Shropshire', 'freely-moving', 'none', 'puretones', '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving', {'20241121_puretones'};
    'Shropshire', 'freely-moving', 'none', 'none', '/media/nas8/OB_ferret_AG_BM/Shropshire/freely-moving', {'20241120_no_sound'};
    
    % Brynza head-fixed
    'Brynza', 'head-fixed', 'none', 'none', '/media/nas7/React_Passive_AG/OBG/Brynza/head-fixed', {'20240124', '20240125', '20240126', '20240129', '20240204', '20240205', '20240305', '20240307', '20240308', '20240313'};
    'Brynza', 'head-fixed', 'atropine', 'none', '/media/nas7/React_Passive_AG/OBG/Brynza/head-fixed', {'20240131_atropine'};
    'Brynza', 'head-fixed', 'domitor', 'none', '/media/nas7/React_Passive_AG/OBG/Brynza/head-fixed', {'20240312_domitor', '20240410_domitor'};
    
    % Brynza freely-moving
    'Brynza', 'freely-moving', 'none', 'none', '/media/nas7/React_Passive_AG/OBG/Brynza/freely-moving', {'20240123_long', '20240124', '20240125', '20240205'}; % 20240123_short % 20240122 not long enough
    'Brynza', 'freely-moving', 'atropine', 'none', '/media/nas7/React_Passive_AG/OBG/Brynza/freely-moving', {'20240130_atropine', '20240201_atropine', '20240306_atropine'};
    'Brynza', 'freely-moving', 'saline', 'none', '/media/nas7/React_Passive_AG/OBG/Brynza/freely-moving', {'20240202_saline'};
    
    % Labneh head-fixed
    'Labneh', 'head-fixed', 'none', 'none', '/media/nas7/React_Passive_AG/OBG/Labneh/head-fixed', { '20230208', '20230225', '20230227', '20230303', '20230307', '20230308', '20230315', '20230321', '20230323', '20230407', '20230418', '20230419', '20230427', '20230504_1', '20230504_2', '20230505_1', '20230505_2', '20230505_3', '20230508_1', '20230508_2', '20230508_3'};
    % removed by BM '20230113', '20230113_1', '20230114', '20230118', '20230121',
    
    % Labneh freely-moving
    'Labneh', 'freely-moving', 'none', 'none', '/media/nas7/React_Passive_AG/OBG/Labneh/freely-moving', {'20221130', '20221201', '20221202_1', '20221202_2', '20221210', '20221212', '20221220', '20221221_long', '20221221_short', '20221222', '20221223', '20221227', '20230309_1', '20230309_2', '20230309_3', '20230621', '20230910_no_EMG_test'};
    'Labneh', 'freely-moving', 'atropine', 'none', '/media/nas7/React_Passive_AG/OBG/Labneh/freely-moving', {'20230206_atropine', '20230210_atropine', '20230218_mustbe_atropine', '20230317_atropine', '20230320_atropine', '20230705_atropine-30kss'};
    'Labneh', 'freely-moving', 'saline', 'none', '/media/nas7/React_Passive_AG/OBG/Labneh/freely-moving', {'20230302_saline', '20230704_saline'};
    'Labneh', 'freely-moving', 'fluoxetine', 'none', '/media/nas7/React_Passive_AG/OBG/Labneh/freely-moving', {'20230628_fluoxetine'};
    
    % Edel head-fixed
    'Edel', 'head-fixed', 'none', 'none', '/media/nas7/React_Passive_AG/OBG/Edel/head-fixed', {'20220415_m', '20220415_n', '20220418_m', '20220418_n', '20220419_m', '20220419_n', '20220420_m', '20220420_n', '20220421_m', '20220421_n', '20220422_m', '20220422_n', '20220426_n', '20220427_m', '20220427_n', '20220428_m', '20220428_n', '20220429_m', '20220429_n', '20220430_m', '20220502_n', '20220503_m', '20220503_n', '20220504_m', '20220504_n', '20220505_m', '20220505_n', '20220506_m', '20220509_n', '20220510_m', '20220510_n', '20220511_m', '20220511_n', '20220512_m', '20220512_n', '20220513_m', '20220517_m', '20220517_n', '20220518_m', '20220518_n', '20220519_m', '20220519_n', '20220520_m', '20220520_n', '20220523_m', '20220523_n', '20220524_m', '20220524_n'};
    };

% Filter sessions based on inputs
for i = 1:size(data, 1)
    if (any(strcmp('all', ferret_name)) || any(strcmp(data{i, 1}, ferret_name))) && ...
            (strcmp(setup_type, 'all') || strcmp(data{i, 2}, setup_type)) && ...
            (strcmp(pharma_protocol, 'all') || strcmp(data{i, 3}, pharma_protocol)) && ...
            (any(strcmp('all', stimuli_protocol)) || any(strcmp(data{i, 4}, stimuli_protocol)))
        
        % Add paths and names to Dir structure
        for j = 1:length(data{i, 6})
            session_name = data{i, 6}{j};
            session_path = fullfile(data{i, 5}, session_name);
            Dir.path{end+1} = session_path; 
            Dir.name{end+1} = session_name;
            
            % Load ExpeInfo if available
            expe_info_path = fullfile(session_path, 'ExpeInfo.mat');
            if exist(expe_info_path, 'file')
                load(expe_info_path, 'ExpeInfo');
                Dir.ExpeInfo{end+1} = ExpeInfo; 
            else
                Dir.ExpeInfo{end+1} = []; 
                warning(['ExpeInfo.mat not found for session: ' session_name]);
            end
        end
    end
end

end
