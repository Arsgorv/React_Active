function Null = itpc_null_from_phasebank(phaseBank, nEv, t_rel, nShuffle)
% Build a null distribution for ITPC by sampling random event times
% restricted so the full [t_rel] window fits in the recording.

if nargin < 4 || isempty(nShuffle), nShuffle = 200; end
fs = phaseBank.fs;
t_start = phaseBank.t_start;
dur = (phaseBank.nSamp - 1) / fs;
t_end = t_start + dur;

tmin = t_start - min(t_rel);
tmax = t_end   - max(t_rel);
if ~(isfinite(tmin) && isfinite(tmax) && tmax > tmin)
    error('itpc_null_from_phasebank: invalid valid-range for event sampling.');
end

nF = numel(phaseBank.fCenters);
nT = numel(t_rel);

acc = zeros(nF,nT);
acc2= zeros(nF,nT);
nOK = 0;

for s = 1:nShuffle
    e = tmin + (tmax - tmin) * rand(nEv,1);
    [itpc, nUsed] = itpc_compute_from_phasebank(phaseBank, e, t_rel, 1);
    if ~isfinite(nUsed(1)) || nUsed(1) < nEv*0.9
        % should not happen due to valid-range sampling, but keep safe
        continue
    end
    x = double(itpc);
    x(~isfinite(x)) = 0;
    acc = acc + x;
    acc2= acc2 + x.^2;
    nOK = nOK + 1;
end

if nOK < 5
    warning('itpc_null_from_phasebank: too few valid shuffles (%d).', nOK);
end

m = acc / max(1,nOK);
v = acc2 / max(1,nOK) - m.^2;
v(v<0) = 0;

Null = struct();
Null.nShuffle = nShuffle;
Null.nOK = nOK;
Null.mean = m;
Null.std  = sqrt(v);
end
