function y = scale_to_range(x, ymin, ymax)
x = double(x(:))';
if all(~isfinite(x))
    y = nan(size(x));
    return
end
xm = nanmin(x); xM = nanmax(x);
if ~isfinite(xm) || ~isfinite(xM) || (xM <= xm)
    y = ymin + 0*x;
else
    y = ymin + (x - xm) * (ymax - ymin) / (xM - xm);
end
end
