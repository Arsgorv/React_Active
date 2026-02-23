function Null = plv_null_from_phasebank(phaseBank, nEv, nShuffle)
% Null PLV by sampling random phase indices uniformly across the recording.

if nargin < 3 || isempty(nShuffle), nShuffle = 500; end

nF = numel(phaseBank.fCenters);
nSamp = phaseBank.nSamp;

acc = zeros(1,nF);
acc2= zeros(1,nF);
nOK = 0;

for s = 1:nShuffle
    idx = randi(nSamp, [nEv 1]);

    plv = nan(1,nF);
    for k = 1:nF
        ph = phaseBank.phCell{k};
        z = mean(exp(1i * ph(idx)));
        plv(k) = abs(z);
    end

    x = double(plv);
    x(~isfinite(x)) = 0;
    acc = acc + x;
    acc2= acc2 + x.^2;
    nOK = nOK + 1;
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
