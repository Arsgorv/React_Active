function [t_s, y] = get_marker_xy(name, D, P)
t_s = []; y = [];
if isfield(P, name)
    sig = P.(name);
elseif isfield(D, name)
    sig = D.(name);
else
    return
end

try
    t_s = double(Range(sig,'s')); y = double(Data(sig));
    return
catch
end

if isstruct(sig)
    if isfield(sig,'t'), t = double(sig.t(:)); else, t = []; end
    if isfield(sig,'data'), y = double(sig.data(:));
    elseif isfield(sig,'Data'), y = double(sig.Data(:));
    else, y = []; end
    if ~isempty(t)
        if nanmedian(diff(t)) > 5
            t_s = t/1e4;
        else
            t_s = t;
        end
    end
end
end