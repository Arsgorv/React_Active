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
