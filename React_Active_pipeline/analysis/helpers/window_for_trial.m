function wsec = window_for_trial(wName, tr, stim_on_rel_s, stim_off_rel_s, arrival_rel_s, arr_anchor)

st = stim_on_rel_s(tr);
sp = stim_off_rel_s(tr);

% choose arrival time used for defining the window
arr = arrival_rel_s(tr);
if nargin >= 6 && isfinite(arr_anchor)
    arr = arr_anchor; % fixed, non-leaky anchor
end

switch wName
    case 'stimon_to_stimoff'
        wsec = [st sp];
    case 'stimon_to_arrival'
        wsec = [st arr];
    case 'stimoff_to_arrival'
        wsec = [sp arr];
    case 'arrival_m500_to_arrival'
        wsec = [arr-0.5 arr];
    otherwise
        wsec = [nan nan];
end
end
