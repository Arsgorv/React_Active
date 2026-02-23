function eventlocked_spec_plot_fromSpectro(Spectro, event_abs_s, opts, ttl)

Praw = Spectro{1};
t0 = Spectro{2}(:);
f  = Spectro{3}(:);

% --- orient + clean + convert time to seconds ---
[t, P] = standardize_spectro_time_and_orient(t0, f, Praw); % t seconds, P time x freq

% log power
P = log10(max(double(P), eps));

% restrict frequency for plotting
fmask = (f >= opts.spec_freq_xlim(1)) & (f <= opts.spec_freq_xlim(2));
if ~any(fmask), error('No freqs in freq_xlim'); end
f_plot = f(fmask);
V = P(:, fmask); % time x freq

% validate interpolation grid
if numel(t) ~= size(V,1)
    error('t and V mismatch: numel(t)=%d size(V,1)=%d', numel(t), size(V,1));
end
if any(~isfinite(t)) || any(diff(t) <= 0)
    error('t is not strictly increasing finite after cleanup.');
end

dt = median(diff(t));
t_rel = (opts.spec_twin_s(1):dt:opts.spec_twin_s(2))';

event_abs_s = event_abs_s(:);
event_abs_s = event_abs_s(isfinite(event_abs_s));
if isempty(event_abs_s)
    warning('No events.'); return
end

X = nan(numel(event_abs_s), numel(t_rel), numel(f_plot));
nanFrac = nan(numel(event_abs_s),1);
nTotal = numel(event_abs_s);
nUsed  = sum(nanFrac < 1);

for k = 1:numel(event_abs_s)
    tq = event_abs_s(k) + t_rel;
    Xi = interp1(t, V, tq, 'linear', NaN); % Xi: numel(tq) x nF
    nanFrac(k) = mean(isnan(Xi(:)));
    X(k,:,:) = Xi;
end

M = squeeze(nanmean(X,1)); % t_rel x f_plot

% baseline correction per frequency
bmask = (t_rel >= opts.spec_baseline_s(1)) & (t_rel <= opts.spec_baseline_s(2));
if any(bmask)
    b = nanmean(M(bmask,:),1);
    M = M - repmat(b, size(M,1), 1);
end

% choose band based on freq range
if max(f_plot) > 20
    band_hz = opts.spec_band_gamma;
else
    band_hz = opts.spec_band_slow;
end
band_mask = (f_plot >= band_hz(1)) & (f_plot <= band_hz(2));
band_trace = nanmean(M(:,band_mask),2);

figure('Color','w', ...
       'Units','pixels', ...
       'Position',[100 50 400 900]); 
   annotation('textbox', [0.62 0.94 0.35 0.05], ...
       'String', sprintf('Trials used: %d / %d', nUsed, nTotal), ...
       'EdgeColor','none', ...
       'HorizontalAlignment','right', ...
       'FontSize',10);
subplot(3,1,1);
imagesc(t_rel, f_plot, M'); axis xy ; colormap('viridis')
ylabel('Frequency (Hz)')
hold on; plot([0 0],[min(f_plot) max(f_plot)],'k-','LineWidth',1); hold off
title(ttl, 'Interpreter','none'); box off

subplot(3,1,2);
plot(t_rel, band_trace, 'k', 'LineWidth',2); hold on
plot([0 0],[min(band_trace) max(band_trace)],'k-','LineWidth',1); hold off
ylabel(sprintf('[%g %g] Hz (baseline-corr log10)', band_hz(1), band_hz(2)))
xlim([t_rel(1) t_rel(end)]); box off

subplot(3,1,3);
histogram(nanFrac, 20);
xlabel('NaN fraction per event'); ylabel('Count'); box off

end
