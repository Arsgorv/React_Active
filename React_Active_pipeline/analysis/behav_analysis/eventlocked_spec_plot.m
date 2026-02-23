function eventlocked_spec_plot(spectrum_file, event_abs_s, twin_s, band_hz, freq_xlim)

S = load(spectrum_file,'Spectro');
Spectro = S.Spectro;

Praw = Spectro{1};
t = Spectro{2}(:);     % seconds (expected)
f = Spectro{3}(:);     % Hz

% --- P orientation check: make P be time x freq ---
% Common cases:
%   size(Praw) = [nT nF]  (time x freq)  -> OK
%   size(Praw) = [nF nT]  (freq x time)  -> transpose
nT = numel(t);
nF = numel(f);

if isequal(size(Praw), [nT nF])
    P = Praw;
elseif isequal(size(Praw), [nF nT])
    P = Praw.'; % -> time x freq
else
    % fallback: guess by matching one dimension to t
    if size(Praw,1) == nT
        P = Praw;         % rows are time
    elseif size(Praw,2) == nT
        P = Praw.';       % columns are time
    else
        error('Spectro dims mismatch: size(P)=[%d %d], numel(t)=%d, numel(f)=%d', ...
            size(Praw,1), size(Praw,2), nT, nF);
    end
end

P = log10(max(double(P), eps));  % time x freq, numeric

% --- Clean time vector t: finite, sorted, unique, and keep P aligned ---
good = isfinite(t);
t = t(good);
P = P(good,:);

% ensure monotonic increasing
[ts, ord] = sort(t);
P = P(ord,:);
t = ts;

% remove duplicates in t (interp1/griddedInterpolant can choke on duplicates)
% keep first occurrence after sort
[t, ia] = unique(t, 'stable');
P = P(ia,:);

% sanity check
if numel(t) ~= size(P,1)
    error('After cleanup: numel(t)=%d but size(P,1)=%d', numel(t), size(P,1));
end
if numel(t) < 2
    error('Not enough time points for interpolation in %s', spectrum_file);
end

dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    error('Bad dt from t in %s (dt=%g)', spectrum_file, dt);
end

t_rel = (twin_s(1):dt:twin_s(2))';
fmask = (f >= freq_xlim(1)) & (f <= freq_xlim(2));

if ~any(fmask)
    error('No frequencies in freq_xlim=[%g %g] Hz. f range=[%g %g] Hz (nF=%d).', ...
        freq_xlim(1), freq_xlim(2), min(f), max(f), numel(f));
end
f_plot = f(fmask);

event_abs_s = event_abs_s(:);
event_abs_s = event_abs_s(isfinite(event_abs_s));
if isempty(event_abs_s)
    warning('No events for %s', spectrum_file);
    return
end

X = nan(numel(event_abs_s), numel(t_rel), numel(f_plot));

V = P(:,fmask); % time x selected freq

for k = 1:numel(event_abs_s)
    tq = event_abs_s(k) + t_rel;
    % interp1 over time dimension (rows of V)
    X(k,:,:) = interp1(t, V, tq, 'linear', NaN);
end

M = squeeze(nanmean(X,1)); % nTrel x nFplot
band_mask = (f_plot >= band_hz(1)) & (f_plot <= band_hz(2));
if ~any(band_mask)
    error('No frequencies in band_hz=[%g %g] Hz inside selected f_plot range=[%g %g].', ...
        band_hz(1), band_hz(2), min(f_plot), max(f_plot));
end
band_trace = nanmean(M(:,band_mask),2);

figure('Color','w');
subplot(2,1,1);
imagesc(t_rel, f_plot, M'); axis xy
ylabel('Frequency (Hz)')
xlim([t_rel(1) t_rel(end)])
ylim(freq_xlim)
hold on; plot([0 0], freq_xlim, 'k-', 'LineWidth',1); hold off
box off

subplot(2,1,2);
plot(t_rel, band_trace, 'k', 'LineWidth',2); hold on
plot([0 0], [min(band_trace) max(band_trace)], 'k-', 'LineWidth',1)
hold off
xlabel('time (s)')
ylabel(sprintf('Band [%g %g] Hz power (log10, a.u.)', band_hz(1), band_hz(2)))
xlim([t_rel(1) t_rel(end)])
box off

end
