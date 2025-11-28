function Dir = PathForExperimentsReactActive(ferret_name, setup_type, pharma_protocol, stimuli_protocol)

% PathForExperimentsReactActive - Retrieves session paths and metadata for OB experiments.
% You can select several ferrets and several stimuli_protocols at the same time

% Inputs:
%   ferret_name     - {'all', 'Tvorozhok', 'Kosichka', 'Mochi', 'Brayon'...} or cell array of ferret names
%   setup_type      - {'all', 'training', 'experiment'} (default: 'all')
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
    % Tvorozhok ;     
        % removed: '20250813_m'
        % 'Tvorozhok', 'head-fixed', 'none', 'all', 'perf_drop' 'Z:\Arsenii\React_Active\training\Tvorozhok\' {'20250723_m', '20250724_m', '20250724_n', '20250725_m', '20250725_n', '20250728_m', '20250729_m', '20250729_n',  '20250730_m', '20250730_n',  '20250731_m',  '20250731_n',  '20250801_m', '20250805_n', '20250806_m', '20250806_n'};
    'Tvorozhok', 'training', 'none', '11_10', 'Z:\Arsenii\React_Active\training\Tvorozhok\' {'20250712', '20250716_m', '20250716_n', '20250717_m', '20250718_n', '20250722_m', '20250722_n'};
    'Tvorozhok', 'training', 'none', '13_12', 'Z:\Arsenii\React_Active\training\Tvorozhok\' {'20250723_m', '20250724_m', '20250724_n', '20250725_m',    '20250725_n'    ,    '20250728_m'    , '20250729_m', '20250729_n',    '20250730_m'    , '20250730_n',...
                                                                                                    '20250731_m'    , '20250731_n', '20250801_m', '20250805_n', '20250806_m', '20250806_n', '20250808_n', '20250813_n', '20250814_m',...
                                                                                               '20250814_n', '20250815_m', '20250815_n', '20250816_m', '20250818_m', '20250818_n', '20250820_m', '20250820_n', '20250821_m', '20250821_n',...
                                                                                               '20250822_m', '20250822_n'};
    'Tvorozhok', 'experiment', 'none', '14_15', 'Z:\Arsenii\React_Active\React_Active_AG\Tvorozhok\' {'20250919', '20250925', '20250929', '20250930', '20251001', '20251002', '20251003'};
                                                                                               
    
    % Mochi
    'Mochi', 'training', 'none', '11_10', 'Z:\Arsenii\React_Active\training\Mochi\' {'20251031_m', '20251031_n', '20251103_m', '20251103_n', '20251104_m', '20251106_m','20251106_n', '20251107_m', '20251107_n_1', '20251107_n_2'};
    'Mochi', 'training', 'none', '13_12', 'Z:\Arsenii\React_Active\training\Mochi\' {'20251111_n', '20251112_m','20251112_n', '20251113_m', '20251113_n', '20251114_m', '20251114_n', '20251115_m', '20251117_m', '20251117_n', '20251118_m', '20251118_n'};
    'Mochi', 'training', 'none', '14_15', 'Z:\Arsenii\React_Active\training\Mochi\' {'20251119_m', '20251119_n', '20251120_m','20251120_n', '20251124_m', '20251125_m', '20251125_n', '20251126_m', '20251126_n'};
    'Mochi', 'training', 'none', '16_17', 'Z:\Arsenii\React_Active\training\Mochi\' {'20251127_m', '20251128_n'};
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
