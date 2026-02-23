function mat = extract_matrix(t_s, y, onsets_abs_s, tvec, baseline_samples, smoothSamples)
% extract_matrix
% Interpolate continuous signal on each trial and build trial x time matrix.
% Handles duplicate timestamps and multi-column signals.
%
% If y is NxM with M>1, uses Euclidean norm across columns.

nT = numel(onsets_abs_s);
nS = numel(tvec);
mat = nan(nT, nS);

t_s = double(t_s(:));
y = double(y);
if isvector(y)
    y = y(:);
end

% Make time unique; average duplicates
[t_s, y] = ra_make_unique_xy(t_s, y);
if numel(t_s) < 2
    return
end

% Collapse multi-col signal to scalar (norm)
if size(y,2) > 1
    y = sqrt(sum(y.^2, 2));
end

tvec = double(tvec(:))';

for k = 1:nT
    t0 = onsets_abs_s(k);
    if ~isfinite(t0), continue; end

    tq = t0 + tvec(:);
    v  = interp1(t_s, y, tq, 'linear', NaN);

    if baseline_samples > 0 && baseline_samples <= numel(v)
        b = mean(v(1:baseline_samples), 'omitnan');
        v = v - b;
    end

    if smoothSamples > 1
        v = movmean(v, smoothSamples, 'omitnan');
    end

    mat(k,:) = v(:)';
end

end
