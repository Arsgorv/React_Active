function RES = RAA_collect_ob_bandtraces_across_sessions(sessions, opts)
% Collect session x time band-trace matrices from saved *_obSpecMetrics.mat
%
% Outputs:
% RES.band(b).name
% RES.band(b).group(g).name
% RES.band(b).group(g).t_rel
% RES.band(b).group(g).M   [nSess x nT] mean trace
% RES.band(b).group(g).N   [nSess x 1]  n events for that group

if nargin < 2, opts = struct(); end
if ~isfield(opts,'alignName'), opts.alignName = 'stimOn'; end
if ~isfield(opts,'groups'), opts.groups = {'all','tar','ref','nosound','nomotor'}; end

RES = struct();
RES.sessions = strings(numel(sessions),1);
RES.datapaths = strings(numel(sessions),1);
RES.alignName = string(opts.alignName);
RES.opts_snapshot = opts;

% We will collect from LOW and MID separately, then merge by band name.
allBands = struct(); % dynamic

for s = 1:numel(sessions)
    dp = sessions{s};
    RES.sessions(s) = string(get_sessname(dp));
    RES.datapaths(s) = string(dp);

    outDir = fullfile(dp,'analysis','ob_events');
    if ~exist(outDir,'dir'), continue, end

    % Find the saved Met files
    fLow = dir(fullfile(outDir, sprintf('*_B_LowEvent_Spectrum_%s_obSpecMetrics.mat', opts.alignName)));
    fMid = dir(fullfile(outDir, sprintf('*_B_Middle_Spectrum_%s_obSpecMetrics.mat', opts.alignName)));

    files = {};
    if ~isempty(fLow), files{end+1} = fullfile(fLow(1).folder, fLow(1).name); end %#ok<AGROW>
    if ~isempty(fMid), files{end+1} = fullfile(fMid(1).folder, fMid(1).name); end %#ok<AGROW>
    if isempty(files), continue, end

    for fi = 1:numel(files)
        A = load(files{fi}, 'Met');
        Met = A.Met;

        % groups present in Met.groupNames
        gnames = string(Met.groupNames(:));
        for gi = 1:numel(opts.groups)
            gWant = string(opts.groups{gi});
            gIdx = find(strcmpi(gnames, gWant), 1, 'first');
            if isempty(gIdx), continue, end

            nEv = Met.group(gIdx).nEvents;

            for b = 1:numel(Met.metricBandNames)
                bname = string(Met.metricBandNames{b});
                t_rel = Met.group(gIdx).band(b).trace.t_rel(:);
                y = Met.group(gIdx).band(b).trace.mean(:);

                key = char(bname);
                if ~isfield(allBands, key)
                    allBands.(key) = struct();
                    allBands.(key).name = bname;
                    allBands.(key).groups = struct();
                end

                gkey = char(gWant);
                if ~isfield(allBands.(key).groups, gkey)
                    allBands.(key).groups.(gkey) = struct();
                    allBands.(key).groups.(gkey).name = gWant;
                    allBands.(key).groups.(gkey).t_rel = t_rel;
                    allBands.(key).groups.(gkey).M = nan(numel(sessions), numel(t_rel));
                    allBands.(key).groups.(gkey).N = nan(numel(sessions),1);
                else
                    % sanity: ensure same time base; if not, interpolate to stored t_rel
                    t0 = allBands.(key).groups.(gkey).t_rel(:);
                    if numel(t_rel) ~= numel(t0) || any(abs(t_rel - t0) > 1e-9)
                        y = interp1(t_rel, y, t0, 'linear', NaN);
                        t_rel = t0;
                    end
                end

                allBands.(key).groups.(gkey).M(s,:) = y(:)';
                allBands.(key).groups.(gkey).N(s) = nEv;
            end
        end
    end
end

% pack into RES.band array
bn = fieldnames(allBands);
RES.band = struct([]);
for i = 1:numel(bn)
    B = allBands.(bn{i});
    RES.band(i).name = B.name;
    gfn = fieldnames(B.groups);
    RES.band(i).group = struct([]);
    for j = 1:numel(gfn)
        G = B.groups.(gfn{j});
        RES.band(i).group(j).name = G.name;
        RES.band(i).group(j).t_rel = G.t_rel(:);
        RES.band(i).group(j).M = G.M;
        RES.band(i).group(j).N = G.N;
    end
end

end

function sn = get_sessname(dp)
[~,sn] = fileparts(dp);
sn = char(sn);
end
