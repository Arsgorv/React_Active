function opts = ra_default_opts(opts)

if nargin < 1 || isempty(opts), opts = struct(); end

% General
if ~isfield(opts,'do_ob_event_analysis'),   opts.do_ob_event_analysis = true; end
if ~isfield(opts,'min_events'),             opts.min_events = 5; end

% ------------------------------
% Event-locked spectrogram
% ------------------------------
if ~isfield(opts,'spec_twin_s'),      opts.spec_twin_s = [-0.5 7]; end
if ~isfield(opts,'spec_twin_whole_s'),   opts.spec_twin_whole_s   = [-0.5 7.0]; end
if ~isfield(opts,'spec_baseline_s'),  opts.spec_baseline_s = [-0.5 -0.1]; end
if ~isfield(opts,'spec_smooth_t'),    opts.spec_smooth_t = 0; end
if ~isfield(opts,'spec_smooth_f'),    opts.spec_smooth_f = 0; end

if ~isfield(opts,'piezo_mode'),          opts.piezo_mode          = 'both'; end % 'resp'|'sniff'|'both'
if ~isfield(opts,'resp_band_hz'),        opts.resp_band_hz        = [0.1 1]; end
if ~isfield(opts,'sniff_band_hz'),       opts.sniff_band_hz       = [1 10]; end
if ~isfield(opts,'resp_smooth_s'),       opts.resp_smooth_s       = 0; end
if ~isfield(opts,'sniff_smooth_s'),      opts.sniff_smooth_s      = 0.2; end

if ~isfield(opts,'xcorr_maxLag_s'),     opts.xcorr_maxLag_s = 2; end
if ~isfield(opts,'xcorr_use_win'),      opts.xcorr_use_win = [0 6.5]; end % for legacy/global xcorr
if ~isfield(opts,'overlay_use_median'), opts.overlay_use_median = true; end

if ~isfield(opts,'phase_win_s'),   opts.phase_win_s = [0 6.5]; end
if ~isfield(opts,'phase_nBins'),   opts.phase_nBins = 18; end


% Frequency windows per spectrum type
if ~isfield(opts,'spec_freq_xlim_low'),  opts.spec_freq_xlim_low = [0.5 10]; end
if ~isfield(opts,'spec_freq_xlim_mid'),  opts.spec_freq_xlim_mid = [20 100]; end

% OB bands (explicit)
if ~isfield(opts,'ob_band_delta'),      opts.ob_band_delta = [0.5 4]; end
if ~isfield(opts,'ob_band_theta'),      opts.ob_band_theta = [4 8]; end
if ~isfield(opts,'ob_band_gamma'),      opts.ob_band_gamma = [40 60]; end
if ~isfield(opts,'ob_band_highgamma'),  opts.ob_band_highgamma = [60 80]; end

% ------------------------------
% ITPC / PLV
% ------------------------------
if ~isfield(opts,'itpc_twin_s'),    opts.itpc_twin_s = [-0.5 7]; end
if ~isfield(opts,'itpc_f_centers'), opts.itpc_f_centers = 2:2:120; end
if ~isfield(opts,'itpc_bw_hz'),     opts.itpc_bw_hz = 2; end
if ~isfield(opts,'itpc_dt_s'),      opts.itpc_dt_s = 1/1024; end

if ~isfield(opts,'plv_f_centers'),  opts.plv_f_centers = 2:2:120; end
if ~isfield(opts,'plv_bw_hz'),      opts.plv_bw_hz = 2; end

if ~isfield(opts,'relWinEdges_s'), opts.relWinEdges_s = []; end
if ~isfield(opts,'relWinNames'),   opts.relWinNames = {}; end

% behaviour-style figure settings
if ~isfield(opts,'make_behaviour_style'), opts.make_behaviour_style = false; end
if ~isfield(opts,'behav_window'),         opts.behav_window = 'stimoff_to_arr'; end % one of the winDefs
if ~isfield(opts,'behav_groups'),         opts.behav_groups = {'all'}; end         % which groups to plot
if ~isfield(opts,'behav_bands'),          opts.behav_bands  = {'delta','theta','gamma'}; end % rows
if ~isfield(opts,'behav_sw_s'),           opts.behav_sw_s   = 0.10; end

end
