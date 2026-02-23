function [Report, remove_sess] = RA_check_required_outputs(rootDir, opts)
% RA_check_required_outputs_animal
% rootDir = animal folder containing session folders as direct children.
% Checks required outputs in:
%   session\analysis\ob_events\
%   session\analysis\behaviour\Pupil-EyeCam\ (variants)
%   session\ephys\
%
% INPUT
%   rootDir : e.g. 'Z:\Arsenii\React_Active\training\Tvorozhok'
%   opts (optional):
%     .ignore_tilde_sessions : true (default) -> ignore session folders starting with '~'
%     .save_csv              : '' (default) or full path to save CSV
%
% OUTPUT
%   Report struct array:
%     .session_path
%     .session_name
%     .ok
%     .missing_files   (cellstr of relative paths / messages)

if nargin < 2, opts = struct(); end
if ~isfield(opts,'ignore_tilde_sessions'), opts.ignore_tilde_sessions = true; end
if ~isfield(opts,'save_csv'),              opts.save_csv = ''; end

sessionPaths = RA_list_session_children(rootDir, opts.ignore_tilde_sessions);

nSess = numel(sessionPaths);
Report = repmat(struct('session_path','', 'session_name','', 'ok',false, 'missing_files',{{}}), nSess, 1);

for i = 1:nSess
    dp = sessionPaths{i};
    [~, sessname] = fileparts(dp);

    missing = {};

    % ------------------------------
    % A) analysis\ob_events
    % ------------------------------
    obDir = fullfile(dp, 'analysis', 'ob_events');
    if exist(obDir, 'dir') ~= 7
        missing{end+1,1} = 'Missing folder: analysis\ob_events'; %#ok<AGROW>
    else
        req_ob = { ...
            sprintf('%s_B_LowEvent_Spectrum_stimOn_obSpecMetrics.mat', sessname), ...
            sprintf('%s_B_Middle_Spectrum_stimOn_obSpecMetrics.mat', sessname), ...
            sprintf('%s_ITPC_stimOn_LOW_metrics.mat', sessname), ...
            sprintf('%s_ITPC_stimOn_MID_metrics.mat', sessname), ...
            sprintf('%s_PLV_stimOn_metrics.mat', sessname) ...
            };
        for k = 1:numel(req_ob)
            rel = fullfile('analysis','ob_events', req_ob{k});
            if exist(fullfile(dp, rel), 'file') ~= 2
                missing{end+1,1} = rel; %#ok<AGROW>
            end
        end
    end

    % ------------------------------
    % B) analysis\behaviour\Pupil-EyeCam (variants)
    % ------------------------------

    req_beh_file = sprintf('%s_Pupil-EyeCam_regular_stimoff_to_arrival_sw0.100_metrics.mat', sessname);
    
    behRoot = fullfile(dp, 'analysis', 'behaviour');
    if exist(behRoot, 'dir') ~= 7
        missing{end+1,1} = 'Missing folder: analysis\behaviour';
    else
        pupilDirs = { ...
            fullfile('analysis','behaviour','Pupil-EyeCam'), ...
            fullfile('analysis','behaviour','Pupil_EyeCam'), ...
            fullfile('analysis','behaviour','Pupil','EyeCam'), ...
            fullfile('analysis','behaviour','Pupil','Eye_EyeCam'), ...
            fullfile('analysis','behaviour','Pupil','Eye-EyeCam') ...
            };
        
        found = false;
        for k = 1:numel(pupilDirs)
            rel = fullfile(pupilDirs{k}, req_beh_file);
            if exist(fullfile(dp, rel), 'file') == 2
                found = true;
                break
            end
        end
        
        if ~found
            rel = fullfile('analysis','behaviour','Pupil-EyeCam', req_beh_file);
            missing{end+1,1} = rel;
        end
    end

% ------------------------------
    % C) ephys\
    % ------------------------------
    req_ephys = { ...
        fullfile('ephys','physio_for_behaviour.mat'), ...
        fullfile('ephys','B_LowEvent_Spectrum.mat') ...
        };
    for k = 1:numel(req_ephys)
        rel = req_ephys{k};
        if exist(fullfile(dp, rel), 'file') ~= 2
            missing{end+1,1} = rel; %#ok<AGROW>
        end
    end

    Report(i).session_path = dp;
    Report(i).session_name = sessname;
    Report(i).missing_files = missing;
    Report(i).ok = isempty(missing);
end

if ~isempty(opts.save_csv)
    RA_write_csv_report(Report, opts.save_csv);
end

% -------------------------------------------------
% Build remove_sess cell array (incomplete sessions)
% -------------------------------------------------

badIdx = find(~[Report.ok]);

remove_sess = cell(numel(badIdx),1);
for i = 1:numel(badIdx)
    remove_sess{i} = Report(badIdx(i)).session_path;
end

% Display in requested format
if ~isempty(remove_sess)
    fprintf('\nremove_sess = {\n');
    for i = 1:numel(remove_sess)
        fprintf('''%s'',...\n', remove_sess{i});
    end
    fprintf('};\n\n');
else
    fprintf('\nAll sessions complete. remove_sess = {};\n\n');
end

end

% ======================================================================

function sessionPaths = RA_list_session_children(rootDir, ignoreTilde)

D = dir(rootDir);
D = D([D.isdir]);
names = {D.name};
names = names(~ismember(names,{'.','..'}));

sessionPaths = cell(0,1);

for k = 1:numel(names)
    nm = names{k};

    if ignoreTilde && ~isempty(nm) && nm(1) == '~'
        continue
    end

    dp = fullfile(rootDir, nm);

    % session root signature in your screenshot:
    % must have both "analysis" and "ephys" dirs
    if exist(fullfile(dp,'analysis'),'dir') == 7 && exist(fullfile(dp,'ephys'),'dir') == 7
        sessionPaths{end+1,1} = dp; %#ok<AGROW>
    end
end

end

function RA_write_csv_report(Report, outCsv)

fid = fopen(outCsv, 'w');
if fid < 0
    warning('Could not open CSV for writing: %s', outCsv);
    return
end

fprintf(fid, 'session_path,session_name,ok,missing_files\n');
for i = 1:numel(Report)
    miss = '';
    if ~isempty(Report(i).missing_files)
        miss = strjoin(Report(i).missing_files, ' || ');
    end
    fprintf(fid, '%s,%s,%s,%s\n', ...
        RA_csv_escape(Report(i).session_path), ...
        RA_csv_escape(Report(i).session_name), ...
        RA_csv_escape(logical_to_str(Report(i).ok)), ...
        RA_csv_escape(miss));
end

fclose(fid);

end

function s = logical_to_str(tf)
if tf, s = 'true'; else, s = 'false'; end
end

function s = RA_csv_escape(s)
if isempty(s), s = ''; return; end
s = strrep(s, '"', '""');
if any(s==',') || any(s=='"') || any(s==char(10)) || any(s==char(13))
    s = ['"' s '"'];
end
end
