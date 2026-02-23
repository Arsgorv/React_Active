function [itpc, nUsed] = itpc_compute_from_phasebank(phaseBank, event_abs_s, t_rel, min_events)

event_abs_s = event_abs_s(:);
event_abs_s = event_abs_s(isfinite(event_abs_s));

fCenters = phaseBank.fCenters;
nF = numel(fCenters);
nT = numel(t_rel);

itpc = nan(nF, nT);
nUsed = nan(nF, 1);

fs = phaseBank.fs;
t_start = phaseBank.t_start;
nSamp = phaseBank.nSamp;

% precompute indices ONCE (independent of frequency)
TQ  = event_abs_s(:) + t_rel(:)';                 % [nEv x nT]
idx = round((TQ - t_start) * fs) + 1;             % [nEv x nT]
ok  = (idx >= 1) & (idx <= nSamp);

% events contributing at least one sample
used_events = all(ok,2);
nUsed_events = sum(used_events);
ok = ok(used_events,:);
idx = idx(used_events,:);

for k = 1:nF
    nUsed(k) = nUsed_events;
    if nUsed(k) < min_events
        continue
    end

    ph = phaseBank.phCell{k};                     % [nSamp x 1]
    Z = nan(size(idx));                           % complex
    Z(ok) = exp(1i * ph(idx(ok)));

    itpc(k,:) = abs(nanmean(Z,1));
end

end
