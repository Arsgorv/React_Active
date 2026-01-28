function t0 = first_sustained_crossing(t, x, th, hold_n)
t0 = nan;
kk = find(x > th);
if isempty(kk), return; end
for j = 1:numel(kk)
    k0 = kk(j);
    if k0+hold_n-1 <= numel(x)
        if all(x(k0:k0+hold_n-1) > th)
            t0 = t(k0);
            return
        end
    end
end
end
