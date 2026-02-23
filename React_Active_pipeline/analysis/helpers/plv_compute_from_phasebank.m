function [plv, pRay, nUse] = plv_compute_from_phasebank(phaseBank, event_abs_s, min_events)

event_abs_s = event_abs_s(:);
event_abs_s = event_abs_s(isfinite(event_abs_s));

fCenters = phaseBank.fCenters;
plv  = nan(size(fCenters));
pRay = nan(size(fCenters));
nUse = nan(size(fCenters));

fs = phaseBank.fs;
t_start = phaseBank.t_start;
nSamp = phaseBank.nSamp;

idx = round((event_abs_s - t_start) * fs) + 1;
ok  = (idx >= 1) & (idx <= nSamp);

idx = idx(ok);
nEv = numel(idx);

for k = 1:numel(fCenters)
    nUse(k) = nEv;
    if nEv < min_events
        continue
    end

    ph = phaseBank.phCell{k};
    ph_ev = ph(idx);

    z = mean(exp(1i*ph_ev));
    plv(k) = abs(z);

    n = nEv;
    R = n*plv(k);
    zRay = (R^2)/n;
    pRay(k) = exp(-zRay) * (1 + (2*zRay - zRay^2)/(4*n) - ...
        (24*zRay - 132*zRay^2 + 76*zRay^3 - 9*zRay^4)/(288*n^2));
end

end
