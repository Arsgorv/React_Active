function OUT = ob_phase_uniformity_within_windows(datapath, opts)
% Computes within-window phase non-uniformity (Rayleigh) per trial, per band, per window.
% Output: per-trial table with MRL (mean resultant length) and Rayleigh p.
%
% Requires:
% - OB LFP channel detectable by load_ob_channel
% - Baphy timing fields (rel_stim_start, rel_stim_stop, spout_arrival_rel_s, rel_trialend)
%
% Windows: stimOn->stimOff, stimOff->arr, arr->stop
% Bands: delta/theta/gamma/highgamma (from opts)

if nargin < 2, opts = struct(); end
opts = ra_default_opts(opts);

[ob_ch, ok] = load_ob_channel(datapath);
if ~ok, error('Cannot determine OB channel'); end

lfpFile = fullfile(datapath,'ephys','LFPData', sprintf('LFP%d.mat', ob_ch));
L = load(lfpFile,'LFP'); LFP = L.LFP;

% Load Baphy
f1 = fullfile(datapath,'Master_sync.mat');
f2 = fullfile(datapath,'Baphy_RA.mat');
if exist(f1,'file'), S = load(f1,'Baphy'); else, S = load(f2,'Baphy'); end
Baphy = S.Baphy;
Ttr = Baphy.trial;
nTr = Baphy.n_trials;

% Trial-relative event times (seconds from stimOn=0)
stimOff = nan(nTr,1);
arr     = nan(nTr,1);
stopT   = nan(nTr,1);

if isfield(Ttr,'rel_stim_start') && isfield(Ttr,'rel_stim_stop')
    stimOff = Ttr.rel_stim_stop(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_stim_stop')
    stimOff = Ttr.rel_stim_stop(:);
end
if isfield(Ttr,'spout_arrival_rel_s') && isfield(Ttr,'rel_stim_start')
    arr = Ttr.spout_arrival_rel_s(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'spout_arrival_rel_s')
    arr = Ttr.spout_arrival_rel_s(:);
end
if isfield(Ttr,'rel_trialend') && isfield(Ttr,'rel_stim_start')
    stopT = Ttr.rel_trialend(:) - Ttr.rel_stim_start(:);
elseif isfield(Ttr,'rel_trialend')
    stopT = Ttr.rel_trialend(:);
end

winNames = {'stimOn_to_stimOff','stimOff_to_arr','arr_to_stop'};
w1 = [zeros(nTr,1), stimOff, arr];
w2 = [stimOff,      arr,     stopT];

% Absolute stimOn times
stimOn_abs = nan(nTr,1);
if isfield(Ttr,'abs_stim_start_s'), stimOn_abs = Ttr.abs_stim_start_s(:); end
if ~any(isfinite(stimOn_abs)), error('Missing abs_stim_start_s'); end

% Build phase bank only at the 4 canonical bands (fast)
bands = [opts.ob_band_delta; opts.ob_band_theta; opts.ob_band_gamma; opts.ob_band_highgamma];
bandNames = {'delta','theta','gamma','highgamma'};

fs = round(1/median(diff(double(Range(LFP,'s')))));
if ~isfinite(fs) || fs<=0, fs = 1024; end

t_abs = double(Range(LFP,'s')); t_abs = t_abs(:);
x = double(Data(LFP)); x = x(:);

% Filter+hilbert per band
phBand = cell(4,1);
for bi = 1:4
    Fil = FilterLFP(LFP, bands(bi,:), fs);
    xb = double(Data(Fil)); xb = xb(:);
    phBand{bi} = angle(hilbert(xb));
end

% Per-trial metrics
rows = {};
for tr = 1:nTr
    if ~isfinite(stimOn_abs(tr)), continue, end
    for wi = 1:3
        t1 = w1(tr,wi); t2 = w2(tr,wi);
        if ~isfinite(t1) || ~isfinite(t2) || t2<=t1, continue, end

        abs1 = stimOn_abs(tr) + t1;
        abs2 = stimOn_abs(tr) + t2;

        i1 = find(t_abs>=abs1, 1, 'first');
        i2 = find(t_abs<=abs2, 1, 'last');
        if isempty(i1) || isempty(i2) || i2<=i1, continue, end

        for bi = 1:4
            ph = phBand{bi}(i1:i2);
            ph = ph(isfinite(ph));
            if numel(ph) < 50, continue, end

            z = mean(exp(1i*ph));
            mrl = abs(z);

            n = numel(ph);
            R = n*mrl;
            zRay = (R^2)/n;
            pRay = exp(-zRay) * (1 + (2*zRay - zRay^2)/(4*n));

            rows(end+1,:) = {tr, string(Ttr.type(tr)), string(winNames{wi}), string(bandNames{bi}), mrl, pRay}; %#ok<AGROW>
        end
    end
end

OUT = struct();
OUT.T = cell2table(rows, 'VariableNames', {'trialIdx','trialType','window','band','MRL','pRayleigh'});
outDir = fullfile(datapath,'analysis','ob_events');
if ~exist(outDir,'dir'), mkdir(outDir); end
[~,sessname] = fileparts(datapath);
writetable(OUT.T, fullfile(outDir, sprintf('%s_OB_phaseUniformity_withinWindows.csv', sessname)));

end
