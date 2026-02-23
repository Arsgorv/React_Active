function ob_event_analysis(datapath, opts)
% ob_event_analysis
% (1) Event-locked mean spectrogram + band-power trace with baseline correction
% (2) Time-resolved ITPC (phase reset) scan around events (plots for stimOn only)
%     + compute-only ITPC for stimOff and arrival
% (3) PLV scan (stimOn)
% (4) Respiration overlay on low-OB spectrogram (stimOn, whole-trial view)

if nargin < 2, opts = struct(); end
opts = ra_default_opts(opts);

% -------- defaults / guards
if ~isfield(opts,'spec_twin_s'),          opts.spec_twin_s = [-0.5 7.0]; end

if ~isfield(opts,'spec_twin_whole_s'),    opts.spec_twin_whole_s = [-0.5 7.0]; end

if ~isfield(opts,'itpc_twin_s'),          opts.itpc_twin_s = [-0.5 7.0]; end
if ~isfield(opts,'itpc_f_centers'),       opts.itpc_f_centers = 2:2:120; end
if ~isfield(opts,'itpc_bw_hz'),           opts.itpc_bw_hz = 2; end

if ~isfield(opts,'relWinEdges_s'), opts.relWinEdges_s = []; end
if ~isfield(opts,'relWinNames'),   opts.relWinNames = {}; end



% -------- Load Baphy once
f1 = fullfile(datapath,'Master_sync.mat');
f2 = fullfile(datapath,'Baphy_RA.mat');
if exist(f1,'file')
    S = load(f1,'Baphy'); Baphy = S.Baphy;
elseif exist(f2,'file')
    S = load(f2,'Baphy'); Baphy = S.Baphy;
else
    error('Missing Master_sync.mat or Baphy_RA.mat');
end

n = Baphy.n_trials;

stim_on_abs_s  = nan(n,1);
stim_off_abs_s = nan(n,1);
arrival_abs_s  = nan(n,1);

if isfield(Baphy,'trial')
    if isfield(Baphy.trial,'abs_stim_start_s'),    stim_on_abs_s  = Baphy.trial.abs_stim_start_s(:); end
    if isfield(Baphy.trial,'abs_stim_stop_s'),     stim_off_abs_s = Baphy.trial.abs_stim_stop_s(:);  end
    if isfield(Baphy.trial,'spout_arrival_abs_s'), arrival_abs_s  = Baphy.trial.spout_arrival_abs_s(:); end
end

outDir = fullfile(datapath,'analysis','ob_events');
if ~exist(outDir,'dir'), mkdir(outDir); end
[~,sessname] = fileparts(datapath);

% -------- Trial groups
G = ra_build_trial_groups(Baphy);

% -------- Events struct (fields expected by ITPC/PLV routines)
isTar = false(n,1);
isRef = false(n,1);
if isfield(Baphy,'trial') && isfield(Baphy.trial,'type')
    tt = string(Baphy.trial.type(:));
    isTar = strcmpi(tt,'Target');
    isRef = strcmpi(tt,'Reference');
end

m_stim_on  = isfinite(stim_on_abs_s);
m_stim_off = isfinite(stim_off_abs_s);
m_arrival  = isfinite(arrival_abs_s);

ev = struct();
ev.stim_on.all   = stim_on_abs_s(m_stim_on);
ev.stim_off.all  = stim_off_abs_s(m_stim_off);
ev.arrival.all   = arrival_abs_s(m_arrival);

ev.stim_on.tar   = stim_on_abs_s(isTar & m_stim_on);
ev.stim_on.ref   = stim_on_abs_s(isRef & m_stim_on);
ev.stim_off.tar  = stim_off_abs_s(isTar & m_stim_off);
ev.stim_off.ref  = stim_off_abs_s(isRef & m_stim_off);
ev.arrival.tar   = arrival_abs_s(isTar & m_arrival);
ev.arrival.ref   = arrival_abs_s(isRef & m_arrival);

% default empty if not available
ev.stim_on.nosound = [];
ev.stim_on.nomotor = [];

if isfield(Baphy,'trial') && isfield(Baphy.trial,'is_nosound')
    ev.stim_on.nosound = stim_on_abs_s(logical(Baphy.trial.is_nosound(:)) & m_stim_on);
end
if isfield(Baphy,'trial') && isfield(Baphy.trial,'is_nomotor')
    ev.stim_on.nomotor = stim_on_abs_s(logical(Baphy.trial.is_nomotor(:)) & m_stim_on);
end

%% (1) Spectrogram grids (stimOn alignment)
alignName   = 'stimOn';
event_abs_s = stim_on_abs_s;

specFiles = { ...
    fullfile(datapath,'ephys','B_LowEvent_Spectrum.mat'), ...
    fullfile(datapath,'ephys','B_Middle_Spectrum.mat') ...
};

for i = 1:numel(specFiles)
    if ~exist(specFiles{i},'file'), continue; end
    [~,sn,~] = fileparts(specFiles{i});

    Ssp = load(specFiles{i},'Spectro');
    Spectro = Ssp.Spectro;

    tmp = opts;

    if contains(sn,'Low','IgnoreCase',true)
        tmp.spec_freq_xlim = tmp.spec_freq_xlim_low;
        tmp.metricBands = [tmp.ob_band_delta; tmp.ob_band_theta];
        tmp.metricBandNames = {'delta','theta'};
        tmp.tag = 'LOW';
    elseif contains(sn,'Middle','IgnoreCase',true)
        tmp.spec_freq_xlim = tmp.spec_freq_xlim_mid;
        tmp.metricBands = [tmp.ob_band_gamma; tmp.ob_band_highgamma];
        tmp.metricBandNames = {'gamma','highgamma'};
        tmp.tag = 'MID';
    else
        continue
    end
    
    tmp.relWinEdges_s = opts.relWinEdges_s;
    tmp.relWinNames   = opts.relWinNames;
    
    ttl = sprintf('%s_%s_%s', sessname, sn, alignName);
    Met = eventlocked_spec_grid_plot(Spectro, Baphy, event_abs_s, G, tmp, ttl, alignName);
    % ---- Export per-trial table for ALL-trials group only (group 1 = 'all')
    gi_all = find(strcmpi(string(Met.groupNames),'all'),1,'first');
    if isempty(gi_all), gi_all = 1; end
    idx_used = Met.group(gi_all).trialIdx_used(:);
    valMat  = Met.group(gi_all).valMat;
    slpMat  = Met.group(gi_all).slpMat;
    
    T = ob_trial_band_metrics_to_table(sessname, alignName, idx_used, Baphy.trial.type(idx_used), ...
    tmp.metricBands, tmp.metricBandNames, Met.winDefs, valMat, slpMat);
    
    save(fullfile(outDir, sprintf('%s_OB_trialMetrics_%s_%s.mat', sessname, alignName, tmp.tag)), 'T');
    writetable(T, fullfile(outDir, sprintf('%s_OB_trialMetrics_%s_%s.csv', sessname, alignName, tmp.tag)));

    saveas(gcf, fullfile(outDir, sprintf('%s_%s_%s_GRID.png', sessname, sn, alignName)));
    saveas(gcf, fullfile(outDir, sprintf('%s_%s_%s_GRID.svg', sessname, sn, alignName)));
    save(fullfile(outDir, sprintf('%s_%s_%s_obSpecMetrics.mat', sessname, sn, alignName)), ...
        'Met', 'tmp', 'alignName');

    close(gcf)
    % ---- (b) Tar vs Ref classification from ALL-trials suppression window
    try
        cls = ob_classify_tar_ref_fromMet(Met, Baphy, tmp);
        save(fullfile(outDir, sprintf('%s_%s_%s_TarRefClass.mat', sessname, sn, alignName)), 'cls', 'tmp');
    catch ME
        warning('Tar/Ref classification failed for %s: %s', sn, ME.message);
    end
end

%% (2) Respiration overlay (stimOn, whole-trial view; low OB spectrogram)
cfg = get_trigger_config(datapath);
opts.resp_lfp_chan = cfg.respi;
try
    lowFile = fullfile(datapath,'ephys','B_LowEvent_Spectrum.mat');
    if exist(lowFile,'file')
        Ssp = load(lowFile,'Spectro');
        SpectroLow = Ssp.Spectro;
        
        ttl = sprintf('%s_Low_respOverlay_stimOn', sessname);
        OutResp = ob_resp_overlay_plot(datapath, SpectroLow, Baphy, stim_on_abs_s, G, opts, ttl);
        for mi = 1:numel(OutResp.mode)
            modeName = OutResp.mode(mi).modeName;
            h1 = OutResp.mode(mi).figure1_handle;
            h2 = OutResp.mode(mi).figure2_handle;
            h3 = OutResp.mode(mi).figure3_handle;
            h4 = OutResp.mode(mi).figure4_handle;
            h5 = OutResp.mode(mi).figure5_handle;
            
            if ishghandle(h1), saveas(h1, fullfile(outDir, sprintf('%s_RESP_%s_fig1_overlay.png', sessname, modeName))); end
            if ishghandle(h2), saveas(h2, fullfile(outDir, sprintf('%s_RESP_%s_fig2_coupling.png', sessname, modeName))); end
            if ishghandle(h3), saveas(h3, fullfile(outDir, sprintf('%s_RESP_%s_fig3_trialxcorr.png', sessname, modeName))); end
            if ishghandle(h4), saveas(h4, fullfile(outDir, sprintf('%s_RESP_%s_fig4_windowxcorr.png', sessname, modeName))); end
            if ishghandle(h5), saveas(h5, fullfile(outDir, sprintf('%s_RESP_%s_fig5_phasepow.png', sessname, modeName))); end
            
            if ishghandle(h1), saveas(h1, fullfile(outDir, sprintf('%s_RESP_%s_fig1_overlay.svg', sessname, modeName))); end
            if ishghandle(h2), saveas(h2, fullfile(outDir, sprintf('%s_RESP_%s_fig2_coupling.svg', sessname, modeName))); end
            if ishghandle(h3), saveas(h3, fullfile(outDir, sprintf('%s_RESP_%s_fig3_trialxcorr.svg', sessname, modeName))); end
            if ishghandle(h4), saveas(h4, fullfile(outDir, sprintf('%s_RESP_%s_fig4_windowxcorr.svg', sessname, modeName))); end
            if ishghandle(h5), saveas(h5, fullfile(outDir, sprintf('%s_RESP_%s_fig5_phasepow.svg', sessname, modeName))); end
            
        end
        
        save(fullfile(outDir, sprintf('%s_RESP_overlay_stimOn_metrics.mat', sessname)), 'OutResp', 'opts', 'G', 'sessname');
        
        % compact CSV summary per group x band x window
        Tall = [];
        TtrAll = [];
        
        for mi = 1:numel(OutResp.mode)
            Tall = [Tall; OutResp.mode(mi).export.summaryTable]; %#ok<AGROW>
            TtrAll = [TtrAll; OutResp.mode(mi).export.trialTable]; %#ok<AGROW>
            TtrAll.trialType = string(Baphy.trial.type(TtrAll.trialIdx));
        end
        
        writetable(Tall,   fullfile(outDir, sprintf('%s_RESP_overlay_stimOn_summary.csv', sessname)));
        writetable(TtrAll, fullfile(outDir, sprintf('%s_RESP_coupling_trials.csv', sessname)));
        for mi = 1:numel(OutResp.mode)
            hs = [OutResp.mode(mi).figure1_handle, OutResp.mode(mi).figure2_handle, ...
                OutResp.mode(mi).figure3_handle, OutResp.mode(mi).figure4_handle, ...
                OutResp.mode(mi).figure5_handle];
            for h = hs
                if ishghandle(h), close(h); end
            end
        end
        
    end
catch ME
    warning('Respiration overlay failed: %s', ME.message);
end

% (e)
% optsB = struct();
% optsB.windowName = 'stimoff_to_arr';
% optsB.bands      = {'delta','theta','gamma','highgamma'};
% optsB.sw_s       = 0.10;
% optsB.ttl        = sprintf('%s | %s', sessname, optsB.windowName);
% fig = RAA_plot_OB_behaviour_style_session(datapath, SpectroLow, Baphy, stim_on_abs_s, G, optsB);
% saveas(fig, fullfile(saveRoot, sprintf('%s_OB_behaviour_%s.png', sessname, optsB.windowName)));

OUT = RAA_ob_bands_behaviour_style(datapath, SpectroLow, Baphy, stim_on_abs_s, G, opts);

close all
 
%% (3) ITPC + (3) PLV
% [ob_ch, ok] = load_ob_channel(datapath);
% if ~ok
%     warning('Could not determine OB channel -> skipping ITPC/PLV.');
%     return
% end
% 
% lfpFile = fullfile(datapath,'ephys','LFPData', sprintf('LFP%d.mat', ob_ch));
% if ~exist(lfpFile,'file')
%     warning('Missing %s -> skipping ITPC/PLV.', lfpFile);
%     return
% end
% 
% L = load(lfpFile,'LFP');
% OB_LFP = L.LFP;
% 
% % ---------------- ITPC MID (plot stimOn only)
% if numel(ev.stim_on.all) >= opts.min_events
%     tmp = opts;
%     tmp.itpc_f_centers   = 2:2:120;
%     tmp.itpc_bw_hz       = 2;
%     tmp.itpc_plot_f_xlim = [20 100];
%     tmp.itpc_band_hz     = tmp.ob_band_gamma;
%     tmp.itpc_band2_hz    = tmp.ob_band_highgamma;
% 
%     phaseBankMid = prepare_phase_bank(OB_LFP, tmp);
% 
%     MetITPC = itpc_scan_plot(OB_LFP, ev.stim_on, G, tmp, 'stim_on_mid', phaseBankMid);
%     saveas(gcf, fullfile(outDir, sprintf('%s_ITPC_stimOn_MID_GRID.png', sessname)));
%     save(fullfile(outDir, sprintf('%s_ITPC_stimOn_MID_metrics.mat', sessname)), ...
%         'MetITPC','tmp','G','sessname');
%     close(gcf)
% 
%     % compute-only stimOff + arrival (MID)
%     try
%         t_rel = (tmp.itpc_twin_s(1):1/phaseBankMid.fs:tmp.itpc_twin_s(2))';
%         MetITPC_mid = struct();
%         MetITPC_mid.t_rel = t_rel;
%         MetITPC_mid.fCenters = phaseBankMid.fCenters;
%         MetITPC_mid.opts_snapshot = tmp;
% 
%         [MetITPC_mid.stim_off.itpc, MetITPC_mid.stim_off.nUsed] = itpc_compute_from_phasebank(phaseBankMid, ev.stim_off.all, t_rel, tmp.min_events);
%         [MetITPC_mid.arrival.itpc,  MetITPC_mid.arrival.nUsed]  = itpc_compute_from_phasebank(phaseBankMid, ev.arrival.all,  t_rel, tmp.min_events);
% 
%         save(fullfile(outDir, sprintf('%s_ITPC_computeOnly_MID.mat', sessname)), 'MetITPC_mid','tmp','G','sessname');
%     catch ME
%         warning('ITPC compute-only MID failed: %s', ME.message);
%     end
% end
% 
% % ---------------- ITPC LOW (plot stimOn only)
% if numel(ev.stim_on.all) >= opts.min_events
%     tmp = opts;
%     tmp.itpc_f_centers   = 0.5:0.5:20;
%     tmp.itpc_bw_hz       = 1;
%     tmp.itpc_plot_f_xlim = [0.5 10];
%     tmp.itpc_band_hz     = tmp.ob_band_delta;
%     tmp.itpc_band2_hz    = tmp.ob_band_theta;
% 
%     phaseBankLow = prepare_phase_bank(OB_LFP, tmp);
% 
%     MetITPC = itpc_scan_plot(OB_LFP, ev.stim_on, G, tmp, 'stim_on_low', phaseBankLow);
%     saveas(gcf, fullfile(outDir, sprintf('%s_ITPC_stimOn_LOW_GRID.png', sessname)));
%     save(fullfile(outDir, sprintf('%s_ITPC_stimOn_LOW_metrics.mat', sessname)), ...
%         'MetITPC','tmp','G','sessname');
%     close(gcf)
% 
%     % compute-only stimOff + arrival (LOW)
%     try
%         t_rel = (tmp.itpc_twin_s(1):1/phaseBankLow.fs:tmp.itpc_twin_s(2))';
%         MetITPC_low = struct();
%         MetITPC_low.t_rel = t_rel;
%         MetITPC_low.fCenters = phaseBankLow.fCenters;
%         MetITPC_low.opts_snapshot = tmp;
% 
%         [MetITPC_low.stim_off.itpc, MetITPC_low.stim_off.nUsed] = itpc_compute_from_phasebank(phaseBankLow, ev.stim_off.all, t_rel, tmp.min_events);
%         [MetITPC_low.arrival.itpc,  MetITPC_low.arrival.nUsed]  = itpc_compute_from_phasebank(phaseBankLow, ev.arrival.all,  t_rel, tmp.min_events);
% 
%         save(fullfile(outDir, sprintf('%s_ITPC_computeOnly_LOW.mat', sessname)), 'MetITPC_low','tmp','G','sessname');
%     catch ME
%         warning('ITPC compute-only LOW failed: %s', ME.message);
%     end
% end
% 
% % ---------------- PLV scan (stimOn)
% if numel(ev.stim_on.all) >= opts.min_events
%     tmp = opts;
%     phaseBankBase = prepare_phase_bank(OB_LFP, tmp);
% 
%     MetPLV = plv_scan_grid_plot(OB_LFP, ev.stim_on, G, tmp, 'stim_on', phaseBankBase);
%     saveas(gcf, fullfile(outDir, sprintf('%s_PLV_GRID_stimOn.png', sessname)));
%     save(fullfile(outDir, sprintf('%s_PLV_stimOn_metrics.mat', sessname)), ...
%         'MetPLV','tmp','G','sessname');
%     close(gcf)
% end

end

function cls = ob_classify_tar_ref_fromMet(Met, Baphy, opts)
% Uses ALL group, window = stimOff_to_arr (default), and returns effect_metrics per band.

bandRanges = Met.metricBands;
bandNames  = Met.metricBandNames;

% find ALL group
gi = find(strcmpi(string(Met.groupNames),'all'),1,'first');
if isempty(gi), gi = 1; end

idx_used = Met.group(gi).trialIdxUsed(:);
if isempty(idx_used)
    cls = struct('winName','', 'band', struct([]));
    return
end

tt = lower(string(Baphy.trial.type(idx_used)));
isTar = contains(tt,'target');
isRef = contains(tt,'reference') | contains(tt,'ref');

% window
winName = 'stimOff_to_arr';
if isfield(opts,'supp_win_mode') && ~isempty(opts.supp_win_mode)
    winName = char(opts.supp_win_mode);
end

wn = string(Met.winDefs(:));
widx = find(strcmpi(wn, winName), 1, 'first');
if isempty(widx)
    error('Requested winName "%s" not found in Met.winDefs.', winName);
end

cls = struct();
cls.winName = winName;
cls.band = struct([]);

for b = 1:size(bandRanges,1)
    v = Met.group(gi).valMat(:, b, widx);
    met = effect_metrics(v(isTar), v(isRef));
    cls.band(b).name  = char(bandNames{b});
    cls.band(b).range = bandRanges(b,:);
    cls.band(b).met   = met;
    cls.band(b).nTar  = sum(isfinite(v(isTar)));
    cls.band(b).nRef  = sum(isfinite(v(isRef)));
end
end
