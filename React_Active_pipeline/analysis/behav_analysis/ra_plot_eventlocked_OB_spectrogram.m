function ra_plot_eventlocked_OB_spectrogram(spectrum_file, event_abs_s, twin_s, band_hz, freq_ylim, clim)
% ra_plot_eventlocked_OB_spectrogram
%
% spectrum_file : e.g. fullfile(datapath,'ephys','B_Middle_Spectrum.mat')
% event_abs_s   : vector of ABSOLUTE event times (seconds, same clock as Spectro{2})
% twin_s        : [t0 t1] relative to event, e.g. [-0.5 5]
% band_hz       : [f0 f1] e.g. [40 60]
% freq_ylim     : [fmin fmax] e.g. [20 100]
% clim          : [] or [cmin cmax] in log10 power

if nargin < 3 || isempty(twin_s), twin_s = [-0.5 5]; end
if nargin < 4 || isempty(band_hz), band_hz = [40 60]; end
if nargin < 5 || isempty(freq_ylim), freq_ylim = [20 100]; end
if nargin < 6, clim = []; end

S = load(spectrum_file,'Spectro');
Spectro = S.Spectro;

Pow = Spectro{1};      % time x freq
t   = Spectro{2}(:);   % sec
f   = Spectro{3}(:);   % Hz

Pow = log10(max(Pow, eps));

dt = median(diff(t));
t_rel = (twin_s(1):dt:twin_s(2))';
event_abs_s = event_abs_s(:);
event_abs_s = event_abs_s(isfinite(event_abs_s));

if isempty(event_abs_s)
    error('No finite event times provided.');
end

fmask = (f >= freq_ylim(1)) & (f <= freq_ylim(2));
f_plot = f(fmask);

X = nan(numel(event_abs_s), numel(t_rel), numel(f_plot));

for k = 1:numel(event_abs_s)
    tq = event_abs_s(k) + t_rel;
    X(k,:,:) = interp1(t, Pow(:,fmask), tq, 'linear', NaN);
end

M = squeeze(nanmean(X,1)); % time x freq

bandMask = (f_plot >= band_hz(1)) & (f_plot <= band_hz(2));
bandTrace = nanmean(M(:,bandMask),2);

figure('Color','w');

ax1 = subplot(2,1,1);
imagesc(t_rel, f_plot, M'); axis xy
ylim(freq_ylim); xlim([t_rel(1) t_rel(end)])
ylabel('Frequency (Hz)')
if ~isempty(clim), caxis(clim); end
hold on
plot([0 0], freq_ylim, 'k', 'LineWidth', 1)
hold off
box off

ax2 = subplot(2,1,2);
plot(t_rel, bandTrace, 'k', 'LineWidth', 2); hold on
plot([0 0], [min(bandTrace) max(bandTrace)], 'k', 'LineWidth', 1)
hold off
xlabel('time (s)')
ylabel(sprintf('OB %g-%g Hz power (log10, a.u.)', band_hz(1), band_hz(2)))
xlim([t_rel(1) t_rel(end)])
box off

linkaxes([ax1 ax2],'x');

end
