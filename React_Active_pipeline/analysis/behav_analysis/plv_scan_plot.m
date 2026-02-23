function plv_scan_plot(LFP, event_abs_s, fCenters, bw, tag)

event_abs_s = event_abs_s(:);
event_abs_s = event_abs_s(isfinite(event_abs_s));

plv  = nan(size(fCenters));
pRay = nan(size(fCenters));
nUse = nan(size(fCenters));

fs = 1024;

t0 = double(Range(LFP,'s'));
if numel(t0) < 2
    warning('plv_scan_plot: LFP too short');
    return
end
t_start = t0(1);

for k = 1:numel(fCenters)
    f1 = max(0.5, fCenters(k) - bw/2);
    f2 = fCenters(k) + bw/2;

    Fil = FilterLFP(LFP, [f1 f2], fs);
    ph  = angle(hilbert(double(Data(Fil))));
    ph  = ph(:);

    % event sample indices
    idx = round((event_abs_s - t_start) * fs) + 1;
    ok  = (idx >= 1) & (idx <= numel(ph));

    ph_ev = ph(idx(ok));
    nUse(k) = numel(ph_ev);

    if nUse(k) < 10
        continue
    end

    z = mean(exp(1i*ph_ev));
    plv(k) = abs(z);

    n = nUse(k);
    R = n*plv(k);
    zRay = (R^2)/n;
    pRay(k) = exp(-zRay) * (1 + (2*zRay - zRay^2)/(4*n) - ...
        (24*zRay - 132*zRay^2 + 76*zRay^3 - 9*zRay^4)/(288*n^2));
end

figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);

subplot(3,1,1);
plot(fCenters, plv, 'k', 'LineWidth',2);
xlabel('Frequency (Hz)'); ylabel('PLV');
title(sprintf('PLV scan (%s)', tag), 'Interpreter','none');
box off

subplot(3,1,2);
plot(fCenters, -log10(pRay), 'k', 'LineWidth',2);
xlabel('Frequency (Hz)'); ylabel('-log10(p) Rayleigh');
box off

subplot(3,1,3);
plot(fCenters, nUse, 'k', 'LineWidth',2);
xlabel('Frequency (Hz)'); ylabel('n events used');
box off

end
