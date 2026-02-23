function v = mean_in_window(row, wsec, fs, n_samples)
v = NaN;
if numel(wsec)~=2 || any(~isfinite(wsec)), return; end
if wsec(2) <= wsec(1), return; end
s0 = max(1, min(n_samples, round(wsec(1)*fs)+1));
s1 = max(1, min(n_samples, round(wsec(2)*fs)+1));
if s1 <= s0, return; end
v = mean(row(s0:s1), 2, 'omitnan');
end

function Baphy = load_baphy(datapath)
f1 = fullfile(datapath,'Master_sync.mat');
f2 = fullfile(datapath,'Baphy_RA.mat');
if exist(f1,'file')
    S = load(f1,'Baphy'); Baphy = S.Baphy;
elseif exist(f2,'file')
    S = load(f2,'Baphy'); Baphy = S.Baphy;
else
    error('Missing Master_sync.mat or Baphy_RA.mat');
end
end
