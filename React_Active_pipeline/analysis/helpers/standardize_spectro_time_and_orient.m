function [t_s, P_tf] = standardize_spectro_time_and_orient(t, f, Praw)

nT = numel(t);
nF = numel(f);

% orient
if isequal(size(Praw), [nT nF])
    P = Praw;
elseif isequal(size(Praw), [nF nT])
    P = Praw.';
else
    if size(Praw,1) == nT
        P = Praw;
    elseif size(Praw,2) == nT
        P = Praw.';
    else
        error('Spectro dims mismatch: size(P)=[%d %d], numel(t)=%d, numel(f)=%d', ...
            size(Praw,1), size(Praw,2), nT, nF);
    end
end

% finite t
good = isfinite(t);
t = t(good);
P = P(good,:);

% sort
[ts, ord] = sort(t);
P = P(ord,:);
t = ts;

% unique
[t, ia] = unique(t, 'stable');
P = P(ia,:);

% drop rows with any non-finite in P (rare but breaks griddedInterpolant logic downstream)
rowok = all(isfinite(P),2);
t = t(rowok);
P = P(rowok,:);

if numel(t) < 2
    error('Not enough points after cleanup.');
end

% infer units (seconds vs ts)
dt0 = median(diff(t));
if dt0 > 5
    t_s = t / 1e4;
else
    t_s = t;
end

P_tf = P;

end
