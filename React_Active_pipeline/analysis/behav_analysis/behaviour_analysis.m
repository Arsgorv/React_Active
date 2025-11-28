function behaviour_analysis(datapath, smoothing_win, varargin)
%{ 
behaviour_analysis(datapath, smoothing_win, Name,Value,...)
EXTRA NAME–VALUE PAIRS
  'TrialSet'   – 'all' (default) | 'NoSound' | 'NoMotor'
  'AntWin'     – struct with fields named exactly like marker variables.
                 Each field = [tStart tEnd] in seconds relative to trial onset.

EXAMPLE
  W.pupil_area_004 = [1.8 3.5];
  behaviour_analysis(dp, 0.1, 'TrialSet','NoSound', 'AntWin',W)


 Summarise behavioural, video and e-phys markers, locked to trial onset.
 Extract each trial into a 0 s -> reward+2 s matrix
 Baseline-correct and (optionally) smooth each trace
 Plot grand means ± SEM and single-trial overlays
%}

%% Initialization
smootime=.006;

% Parse input
p = inputParser;
addParameter(p,'TrialSet','all',@(s) any(strcmpi(s,{'all','nosound','nomotor'})));
% addParameter(p,'AntWin',struct(),@isstruct);
parse(p,varargin{:});
trialSet = lower(p.Results.TrialSet);
% antWin = p.Results.AntWin;

% Load behaviour data
load(fullfile(datapath, 'stim', 'trial_structure.mat'))
load(fullfile(datapath, 'video', 'DLC_data.mat'))
load(fullfile(datapath, 'behavResources.mat'), 'MovAcctsd')
time_video = time_face/1e4;

% Load ephys data
load(fullfile(datapath, 'SleepScoring_OBGamma.mat'), 'BrainPower')
load(fullfile(datapath, 'LFPData', 'LFP19.mat')); respiration_LFP = LFP;
load(fullfile(datapath, 'LFPData', 'LFP6.mat')); EMG_LFP = LFP;
SmoothGamma = BrainPower.Power{1};

%% Build trial masks %add it to baphy_preprocessing script

nTrials = numel(trial_structure.hit_rate.all);
idxTarget = trial_id_struct.Target(:);
idxReference = trial_id_struct.Reference(:);
idxGood = trial_id_struct.goodTrials(:);

idxNoSound = trial_id_struct.NoSound(:);
nsType = trial_id_struct.NoSoundType;

hasNoMotor = isfield(trial_id_struct,'NomotorExclusive');
if hasNoMotor
    idxNoMotor = trial_id_struct.NomotorExclusive(:);
    nmType = trial_id_struct.NomotorExclusiveType;
else
    idxNoMotor = []; nmType = {};
end

% onsets in seconds
tAll = nan(nTrials,1);
tAll(idxReference) = Range(trial_structure.trial_onset.Reference,'s');
tAll(idxTarget) = Range(trial_structure.trial_onset.Target,'s');
tAll(idxNoSound) = Range(trial_structure.trial_onset.NoSound,'s');
if hasNoMotor
    tAll(idxNoMotor) = Range(trial_structure.trial_onset.NomotorExclusive,'s');
end

% logical masks
M = struct(); % container for all masks
M.mskTarget = false(nTrials,1);  M.mskTarget(idxTarget) = true;
M.mskRef = false(nTrials,1); M.mskRef(idxReference) = true;
M.mskGood = false(nTrials,1); M.mskGood(idxGood) = true;

M.mskNS = false(nTrials,1); M.mskNS(idxNoSound) = true;
M.mskNS_T = false(nTrials,1);
M.mskNS_R = false(nTrials,1);
if ~isempty(idxNoSound)
    M.mskNS_T(idxNoSound(strcmp(nsType,'Target'))) = true;
    M.mskNS_R(idxNoSound(strcmp(nsType,'Reference'))) = true;
end

M.mskNM = false(nTrials,1);
M.mskNM_T = false(nTrials,1);
M.mskNM_R = false(nTrials,1);
if hasNoMotor
    M.mskNM(idxNoMotor) = true;
    M.mskNM_T(idxNoMotor(strcmp(nmType,'Target'))) = true;
    M.mskNM_R(idxNoMotor(strcmp(nmType,'Reference'))) = true;
end

M.mskReg = ~(M.mskNS | M.mskNM);

switch trialSet
    case 'nosound', keepMask = M.mskNS;
    case 'nomotor', keepMask = M.mskNM;
    otherwise,      keepMask = true(nTrials,1);
end
keepMask = keepMask & M.mskGood;

idxRef = find(M.mskRef & keepMask);
idxTar = find(M.mskTarget & keepMask);
ref_onsets = tAll(idxRef);
tar_onsets = tAll(idxTar);

%% Resample ephys data (REVIEW)
ds2videoFS = @(sig) tsd(time_video*1e4, interp1(Range(sig,'s'), Data(sig), time_video, 'linear', 'extrap'));

% Process OB (resample); sampled at 1250Hz
SmoothGamma_rs = ds2videoFS(SmoothGamma); 
load(fullfile(datapath, 'LFPData', 'LFP12.mat')); OB_LFP = LFP;
FilLFP_OB = FilterLFP(OB_LFP,[40 60],1024);
Phase_OB=tsd(Range(FilLFP_OB) , angle(hilbert(zscore(Data(FilLFP_OB))))*180/pi+180);
Phase_OB_rs = ds2videoFS(Phase_OB);

% Phase_OB_Above_350=thresholdIntervals(Phase_OB,350,'Direction','Above');
% OB_trig=ts((Stop(Phase_OB_Above_350)+Start(Phase_OB_Above_350))/2);
% OB_trig_rs = ds2videoFS(OB_trig);

%% Process EMG (filter, smooth, resample); sampled at 1250Hz
FilLFP_EMG = FilterLFP(EMG_LFP,[50 300],1024);
emgEnv = runmean(Data((FilLFP_EMG)).^2,ceil(smootime/median(diff(Range(FilLFP_EMG,'s')))));
EMGData = tsd(Range(FilLFP_EMG),emgEnv);
EMGData_rs = ds2videoFS(EMGData);
EMGData_log10 = tsd(Range(FilLFP_EMG),log10(emgEnv));
EMGData_log10_rs = ds2videoFS(EMGData_log10);
Phase_EMG=tsd(Range(FilLFP_EMG) , angle(hilbert(zscore(Data(FilLFP_EMG))))*180/pi+180);
Phase_EMG_rs = ds2videoFS(Phase_EMG);

%% Process accelerometer (filter, smooth, resample) ; sampled at 50Hz
MovAcctsd=tsd(Range(MovAcctsd),runmean(Data((MovAcctsd)).^2,ceil(smootime/median(diff(Range(MovAcctsd,'s'))))));
MovAcctsd_rs = ds2videoFS(MovAcctsd);
MovAcctsd_log10=tsd(Range(MovAcctsd),log10(runmean(Data((MovAcctsd)).^2,ceil(smootime/median(diff(Range(MovAcctsd,'s')))))));
MovAcctsd_log10_rs = ds2videoFS(MovAcctsd_log10);

%% Process respiration (filter, smooth, resample); sampled at 1250Hz

FilLFP_respi = FilterLFP(respiration_LFP,[0.1 1],1024);
respiENV = runmean(Data((FilLFP_respi)).^2,ceil(smootime/median(diff(Range(FilLFP_respi,'s')))));
BreathingPower = tsd(Range(FilLFP_respi),respiENV);
BreathingPower_rs = ds2videoFS(BreathingPower);
BreathingPower_log10 = tsd(Range(FilLFP_respi),log10(respiENV));
BreathingPower_log10_rs = ds2videoFS(BreathingPower_log10);
BreathingPower_var = tsd(Range(BreathingPower) , movstd(Data(BreathingPower) , ceil(30/median(diff(Range(FilLFP_respi,'s'))))));
BreathingPower_var_rs = ds2videoFS(BreathingPower_var);

Phase=tsd(Range(FilLFP_respi) , angle(hilbert(zscore(Data(FilLFP_respi))))*180/pi+180);
Phase_Above_350=thresholdIntervals(Phase,350,'Direction','Above');
Phase_rs = ds2videoFS(Phase);

Sniff=ts((Stop(Phase_Above_350)+Start(Phase_Above_350))/2);
Sniff_rs = ds2videoFS(Sniff);

% UltraLowSpectrumBM(datapath, 'LFP19','Piezzo');
% P = load('Piezzo_UltraLow_Spectrum.mat');
% Sptsd = tsd(P.Spectro{2}*1e4 , P.Spectro{1});
% Breathing = ConvertSpectrum_in_Frequencies_BM(P.Spectro{3} , Range(Sptsd) , Data(Sptsd) , 'frequency_band' , [.25 1]);
% Breathing_rs = ds2videoFS(Breathing);
% Breathing_var = tsd(Range(Breathing) , movstd(Data(Breathing) , ceil(30/median(diff(Range(Breathing,'s')))),'omitnan'));
% Breathing_var_rs = ds2videoFS(Breathing_var);

% Calculate respi frequency
% resp_z = detrend(Data(BreathingPower_rs),'linear');
% 
% % recompute zero-crossings on centred trace
% zeroX  = find(resp_z(1:end-1)<=0 & resp_z(2:end)>0);
% 
% if numel(zeroX) > 1
%     t    = Range(BreathingPower_rs,'s');          % seconds
%     sniffT    = diff(t(zeroX));                    % sec between sniffs
%     sniffRate = 1 ./ sniffT;                       % Hz
% 
%     % interpolate back to full 50 Hz grid
%     sniffRate50 = interp1(t(zeroX(2:end)), sniffRate, t, ...
%                           'linear','extrap');
%     respiration_rate = tsd(t*1e4, sniffRate50);
% else
%     warning('Not enough zero-crossings to compute respiration rate');
%     respiration_rate = tsd(Range(BreathingPower_rs,'s'), ...
%                            nan(size(resp_z)));
% end

% clf
% d = respiration_LFP;
% d1 = respiration_LFP_rs;
% hold all
% plot(Range(d, 's'), Data(d), '*');
% plot(Range(d1, 's'), Data(d1), '*');

%% Define body parts and associated markers ---
bodyparts = {
    'Pupil', {'pupil_area_004','pupil_center_004','pupil_center_004_mvt'};
    'Eye',   {'eye_area_004'};
    'Nostril',{'nostril_area','nostril_center','nostril_center_mvt'};
    'Nose',  {'nose_area','nose_center','nose_center_mvt'};
    'Cheek', {'cheek_center','cheek_center_mvt'};
    'Ear',   {'ear_center','ear_center_mvt'};
    'Jaw',   {'jaw_center','jaw_center_mvt'};
    'Tongue',{'tongue_center','tongue_center_mvt'};
    'Pupil_EyeCam', {'pupil_area_007','pupil_center_007','pupil_center_007_mvt'};
    'Eye_EyeCam',   {'eye_area_007'};
    'Spout', {'spout_likelihood'};
    'OB_gamma_power', {'SmoothGamma','SmoothGamma_rs', 'Phase_OB', 'Phase_OB_rs'};
    'Respiration', {'BreathingPower_rs', 'BreathingPower_log10_rs', 'Phase','Phase_rs'};
%     'Respiration_other', {'BreathingPower_var_rs', 'Breathing', 'Breathing_rs', 'Breathing_var_rs'};
    'EMG', {'EMGData_rs', 'EMGData_log10_rs', 'Phase_EMG', 'Phase_EMG_rs'};
    'Accelerometer', {'MovAcctsd_rs', 'MovAcctsd_log10_rs'};
    };

%% colour setup 
nbp          = size(bodyparts,1);         % number of body parts
baseColors   = lines(nbp);               % distinct, high-contrast roots
maxPerPart   = 5;                        % how many shades we pre-compute
shadeLevels  = linspace(0,0.55,maxPerPart+1);   % 0 = dark root, ? lighter
shadeLevels(1) = [];                     % remove 0 so root stays untouched

% helper: lighten a colour toward white by factor a (0–1)
lighten = @(c,a) c + (1-c).*a;           % same as imadjust(col,[],[],gamma)

fs = 50;
baseline_pre = 0.2; % s
stim_patch_len = 2.1; % s
epoch_post = 2.0;    % s after reward
reward_rel = nanmean(trial_structure.reward_onset.Target); % mean reward time rel to trial onset
trial_dur = reward_rel + epoch_post;
n_samples = round(trial_dur*fs)+1;
baseline_samples = round(baseline_pre*fs);
tvec = (0:n_samples-1)/fs;

%% Per-marker anticipatory window
delayMean = nanmedian(trial_structure.rel_spout_arrival); 
defaultWin = [3.7 delayMean]; % in seconds

% AW.nostril_area = defaultWin;
% AW.nostril_center = defaultWin;
% AW.nostril_center_mvt = defaultWin;

% AW.nose_area = [delayMean reward_rel]; % defaultWin
% AW.nose_center = [delayMean reward_rel]; % defaultWin
% AW.nose_center_mvt = [delayMean reward_rel]; % defaultWin
% 
% AW.cheek_center = [delayMean reward_rel]; % defaultWin
% AW.cheek_center_mvt = [delayMean reward_rel]; % defaultWin

% AW.ear_center = defaultWin;
% AW.ear_center_mvt = defaultWin;
% 
% AW.jaw_center = defaultWin;
% AW.jaw_center_mvt = defaultWin;

% AW.tongue_center = [delayMean reward_rel];
% AW.tongue_center_mvt = defaultWin;

% AW.pupil_area_007 = defaultWin;
% AW.pupil_center_007 = defaultWin;
% AW.pupil_center_007_mvt = defaultWin;

% AW.eye_area_007 = defaultWin;

% AW.SmoothGamma_rs = defaultWin;

% AW.respiration_LFP_rs = defaultWin;
% AW.respiration_rate =defaultWin;

% AW.EMGData_rs = defaultWin;
% AW.MovAcctsd_log10_rs = [0.2 2.3]; % [delayMean-0.1 delayMean]
% AW.MovAcctsd_rs = [0.2 2.3]; % [delayMean-0.1 delayMean]
% 
% antWin = AW;
           
%% Generate figures
metricNames = {'AshmanD','AUROC','CohensD','meanDiff'};
nMetrics    = numel(metricNames);
colsTotal   = 2 + nMetrics;
allMetrics  = struct();                % collect per bodypart

for b = 1:nbp
    bodypart_name = bodyparts{b,1};
    marker_names = bodyparts{b,2};
    baseCol   = baseColors(b,:);
    n_markers = numel(marker_names);
    
    f{b} = figure('Color', 'w', 'visible', 'off'); sgtitle(['Name ' bodypart_name '. Smoothing: ' num2str(smoothing_win) 's. Ref# ' num2str(numel(idxRef)) '. Tar# ' num2str(numel(idxTar))]); 
    set(f{b}, 'Position', [1 1 1980 n_markers*240]);
    set(gcf,'renderer','OpenGL');
    
    for m = 1:n_markers
        % pick colour: first marker = baseCol, next markers = lighter shades
        col = baseCol;
        if m > 1
            a   = shadeLevels(min(m-1,numel(shadeLevels)));
            col = lighten(baseCol,a);
        end

        marker_var = marker_names{m};
        if ~exist(marker_var, 'var')
            warning('%s does not exist in workspace, skipping...', marker_var);
            continue
        end
        marker = eval(marker_var);
        t_marker = Range(marker,'s');
        y_marker = Data(marker);
        
        % ------ marker-specific anticipatory window ------------------------
        if exist('antWin', 'var') && isfield(antWin, marker_var) && numel(antWin.(marker_var))==2
            WinSec = antWin.(marker_var);     % [start  end] in seconds
        else
            WinSec = defaultWin;              % global default
        end
        
        WinSamp = round(WinSec*fs);
        WinSamp  = max(1, min(n_samples, WinSamp));  % clamp to [1, n_samples]
        if diff(WinSamp) <= 0
            warning('%s: anticipatory window invalid/collapsed; skipping stats panels.', marker_var);
        end
        WinStart = WinSec(1);   WinEnd = WinSec(2);
        
        % extract trial matrices ----------------------------------------  
        smoothSamples = max( 1, round(smoothing_win * fs) ); % seconds ? samples
        ref_mat = extractTrialMatrix(t_marker, y_marker, ref_onsets, ...
                                     baseline_samples, n_samples, smoothSamples);
        tar_mat = extractTrialMatrix(t_marker, y_marker, tar_onsets, ...
                                     baseline_samples, n_samples, smoothSamples);
        
        if all(isnan(ref_mat(:))) && all(isnan(tar_mat(:)))
            warning('%s: no valid trials extracted. Skipping subplot.', marker_var);
            continue
        end
        
        % anticipatory means
        mRef = mean(ref_mat(:,WinSamp(1):WinSamp(2)),2,'omitnan');
        mTar = mean(tar_mat(:,WinSamp(1):WinSamp(2)),2,'omitnan');
        mean_ref = nanmean(ref_mat,1);
        mean_tar = nanmean(tar_mat,1);
        sem_ref  = nanstd(ref_mat,[],1)/sqrt(size(ref_mat,1));
        sem_tar  = nanstd(tar_mat,[],1)/sqrt(size(tar_mat,1));
        yl       = [min([mean_ref-sem_ref mean_tar-sem_tar]) , ...
            max([mean_ref+sem_ref mean_tar+sem_ref])];
        
        validRef = mRef(isfinite(mRef));
        validTar = mTar(isfinite(mTar));
                        
        % Metrics 
        metrics = struct('AshmanD',nan,'AUROC',nan,'CohensD',nan,'meanDiff',nan, ...
            'meanRef',nan,'meanTar',nan,'stdRef',nan,'stdTar',nan, ...
            'rocFPR',[],'rocTPR',[]);
        if numel(validRef)>3 && numel(validTar)>3
            metrics.meanRef = mean(validRef);
            metrics.meanTar = mean(validTar);
            metrics.stdRef  = std(validRef);
            metrics.stdTar  = std(validTar);
            % Ashman D
            vr = var(validRef); vt = var(validTar);
            metrics.AshmanD = abs(mean(validRef)-mean(validTar)) / ...
                sqrt(0.5*(vr+vt));
            % AUROC
            allv = [validRef; validTar];
            if numel(unique(allv)) > 1
                lab = [zeros(numel(validRef),1); ones(numel(validTar),1)];
                [rocFPR,rocTPR,~,AUC] = perfcurve(lab, allv, 1);
                metrics.AUROC = AUC; metrics.rocFPR = rocFPR; metrics.rocTPR = rocTPR;
            end
            % Cohen’s d
            sPool = sqrt(((numel(validRef)-1)*vr+(numel(validTar)-1)*vt)...
                /(numel(validRef)+numel(validTar)-2));
            metrics.CohensD = (mean(validTar)-mean(validRef))/sPool;
            % mean difference (Tar-Ref)
            metrics.meanDiff = mean(validTar)-mean(validRef);
        end
        allMetrics.(bodypart_name).(marker_var) = metrics;
        
        %% Time course plot 
        subplot(n_markers,colsTotal,(m-1)*colsTotal + 1); hold on      
        
        % Gray patch = stimulus
        hStim = patch([(stim_patch_len+0.2) 0.2 0.2 (stim_patch_len+0.2)], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.18);
        % Spout
        hSpout = plot([nanmean(trial_structure.rel_spout_arrival) nanmean(trial_structure.rel_spout_arrival)], yl, 'k--','LineWidth',1.2);
        % Reward
        hReward = plot([reward_rel reward_rel], yl, 'r--','LineWidth',1.2);
        % plot
        hRef = shadedErrorBar_BM(tvec, ref_mat, {'color',col,'LineWidth',2}, 1);
        set(hRef.patch,'FaceAlpha',0.15,'EdgeAlpha',0.25); 
        hTar = shadedErrorBar_BM(tvec, tar_mat, {'color',col*0.5,'LineWidth',2}, 1);
        set(hTar.patch,'FaceAlpha',0.15,'EdgeAlpha',0.25);
        hAnticip = patch([WinStart WinEnd WinEnd WinStart], ...
            yl([1 1 2 2]), ...
            [0 1 0],'FaceAlpha',0.08,'EdgeColor','none');
        % Labels
        ylabel(marker_var, 'Interpreter', 'none')
        if m == 1
            title(['Body Part: ' bodypart_name],'FontWeight','bold','FontSize',15)
        end
        if m == n_markers
            xlabel('Time from trial onset (s)')
        end
        lgd = legend([hStim, hAnticip, hRef.mainLine, hTar.mainLine, hSpout, hReward], ...
            {'Stimulus', 'Anticip window', 'Reference','Target', 'Spout', 'Reward'},'Location','westoutside');
        drawnow                  % make sure the legend has a position
        set(lgd,'Units','normalized')
        p = lgd.Position;        % [x y w h] in figure-normalized units
        lgd.Position = [p(1)-0.095  p(2)  p(3)  p(4)];
        xlim([tvec(1) tvec(end)]); ylim(yl); box off; hold off
        
        %% Distributions and Ashmans'D
        subplot(n_markers,colsTotal,(m-1)*colsTotal + 2); cla; hold on
        if ~isempty(validRef), histogram(validRef,30,'Normalization','pdf', ...
                'FaceColor',col,'FaceAlpha',0.4); end
        if ~isempty(validTar), histogram(validTar,30,'Normalization','pdf', ...
                'FaceColor',col*0.5,'FaceAlpha',0.4); end
        
        % KDE curves only if we have enough spread
        if numel(validRef) > 3 && numel(unique(validRef)) > 1
            [fR,xR] = ksdensity(validRef);  plot(xR,fR,'Color',[0.2 0.6 1],'LineWidth',1.2);
        end
        if numel(validTar) > 3 && numel(unique(validTar)) > 1
            [fT,xT] = ksdensity(validTar);  plot(xT,fT,'Color',[1 0.4 0.4],'LineWidth',1.2);
        end
        
        title(sprintf('Ashman D = %.2f', metrics.AshmanD), 'Interpreter','none','FontSize',8)
        axis tight; box off; set(gca,'YAxisLocation','right')
        if m==n_markers
            xlabel(sprintf('Marker value (window %.1f–%.1f s)',WinStart,WinEnd));
        end
        %%  other metrics
        for mc = 2:nMetrics
            subplot(n_markers,colsTotal,(m-1)*colsTotal+1+mc); cla
            MakeSpreadAndBoxPlot3_SB({mRef,mTar},{col,col*0.5},[1 2],...
                {'Ref','Tar'},'newfig',0,...
                'showpoints',1,'optiontest','ranksum','paired',0);
            title(metricNames{mc})
        end
        
        % ---- ROC curve panel -----------------------------------------------
        subplot(n_markers,colsTotal,(m-1)*colsTotal + 1 + 1 + nMetrics); cla
        plot(metrics.rocFPR, metrics.rocTPR, 'Color', col, 'LineWidth',1.5); hold on
        plot([0 1],[0 1],'k:')                    % chance diagonal
        axis square, box off
        xlabel('False-positive rate'), ylabel('True-positive rate')
        title(sprintf('ROC  (AUC = %.2f)', metrics.AUROC))
        
        hold off
    end
     %% Save per-bodypart outputs 
    outDir = fullfile(datapath,'analysis','behaviour',bodypart_name); [~,session_name] = fileparts(datapath);
    if ~exist(outDir,'dir'), mkdir(outDir), end
    saveas(f{b},fullfile(outDir,[session_name '_' bodypart_name '_' trialSet '_sw_' num2str(smoothing_win) '_aw_' num2str(round(WinSec(1))) '_' num2str(round(WinSec(2))) '.png']));
   
    outDir_1 = fullfile(fileparts(datapath), 'DS_figures', 'behaviour'); 
    if ~exist(outDir_1,'dir'), mkdir(outDir_1), end
    saveas(f{b},fullfile(outDir_1,[session_name '_' bodypart_name '_' trialSet '_sw_' num2str(smoothing_win) '_aw_' num2str(round(WinSec(1))) '_' num2str(round(WinSec(2))) '.png']));

    T = struct2table(allMetrics.(bodypart_name));
    writetable(T, fullfile(outDir,[bodypart_name '_' trialSet '_metrics.csv']));
    save(fullfile(outDir,[bodypart_name '_' trialSet '_metrics.mat']),'T')
end
% close all

%% Plot distribution of first licks
% lick_ref = trial_structure.lick_onset.Reference(find(~isnan(trial_structure.lick_onset.Reference)));
% lick_tar = trial_structure.lick_onset.Target(find(~isnan(trial_structure.lick_onset.Target)));
% figure('Name','First-lick latency'); clf
% MakeSpreadAndBoxPlot3_SB({lick_ref, lick_tar},{[1 0 0],[0 0 1]},[2,4],{'Reference','Target'});
% ylabel('Latency from trial onset (s)')
% title('First-lick latency distribution')

%% pieces of BM scripts physio
% % a bit irrelevant
% load(fullfile(datapath, 'LFPData', 'LFP12.mat'))
% FilGamma = FilterLFP(LFP,[40 60],1024);
% tEnveloppeGamma = tsd(Range(LFP), abs(hilbert(Data(FilGamma))) ); 
% smootime=.06;
% SmoothGamma = tsd(Range(tEnveloppeGamma), runmean(Data(tEnveloppeGamma), ...
%     ceil(smootime/median(diff(Range(tEnveloppeGamma,'s'))))));
% 
% 
% FilULow = FilterLFP(LFP,[.1 1],1024); 
% tEnveloppeULow = tsd(Range(LFP), abs(hilbert(Data(FilULow))) ); 
% smootime=.006;
% SmoothULow = tsd(Range(tEnveloppeULow), runmean(Data(tEnveloppeULow), ...
%     ceil(smootime/median(diff(Range(tEnveloppeULow,'s'))))));
% 
% 
% figure
% subplot(312)
% plot(Range(SmoothULow,'s') , Data(SmoothULow))
% hold on
% plot(Range(SmoothGamma,'s') , Data(SmoothGamma))
% % xlim([8249 8259])
% 
% 
% [c_all,lags] = xcorr(zscore(Data(SmoothULow)) , zscore(Data(Restrict(SmoothGamma,SmoothULow))) , 3750);
% 
% figure
% plot(linspace(-3,3,7501) , c_all , 'b' , 'LineWidth' , 2)
% vline(0,'--r')
% xlabel('lag (s)'), ylabel('Corr values (a.u.)'), xlim([-3 3])
% box off
% 
% % Breathing with piezo UL spectro
% UltraLowSpectrumBM(datapath, 'LFP19','Piezzo');
% UltraLowSpectrumBM(datapath, 'LFP12','B');
% 
% % Find a way to calculate Piezzo_ULow_Spectrum
% % Breathing Piezzo, frequency and variability
% P = load('Piezzo_UltraLow_Spectrum.mat');
% Sptsd = tsd(P.Spectro{2}*1e4 , P.Spectro{1});
% % [Sptsd_clean,~,EpochClean] = CleanSpectro(Sptsd , P.Spectro{3} , 8);
% 
% Piezzo  = Sptsd;
% 
% Spectrum_Frequency = ConvertSpectrum_in_Frequencies_BM(P.Spectro{3} , Range(Sptsd) , Data(Sptsd) , 'frequency_band' , [.25 1]);
% Breathing_var = tsd(Range(Spectrum_Frequency) , movstd(Data(Spectrum_Frequency) , ceil(30/median(diff(Range(Spectrum_Frequency,'s')))),'omitnan'));
% 
% Breathing  = Spectrum_Frequency;
% 
% Cols=[0 0 1];
% X = 1;
% Legends = {'All'};
% NoLegends = {'','','',''};
% 
% figure
% D{1} = Data(Breathing)*60;
% MakeSpreadAndBoxPlot3_SB(D,Cols,X,Legends,'showpoints',0,'paired',0)
% ylabel('Breathing / min')
% makepretty_BM2
% 
% figure
% imagesc(Range(Sptsd)/3.6e7 , P.Spectro{3} , runmean(runmean(log10(Data(Sptsd)'),5)',50)'), axis xy
% ylabel('Frequency (Hz)')
% colormap viridis
% caxis([-1.8 -.5])
% 
% figure
% plot(P.Spectro{3} , nanmean(Data(Piezzo)) , 'b')
% 
% figure
% [Y,X]=hist(Data(Breathing_var),1000);
% Y=Y/sum(Y);
% plot(X,runmean(Y,90),'b','LineWidth',1)
% xlabel('Breathing variability (a.u.)'), ylabel('PDF')
% box off,% xlim([1.5 6])
% legend('Wake','NREM','REM')
% 
% % Breathing power var
% load(fullfile(datapath, 'LFPData', 'LFP19.mat'))
% FilLFP = FilterLFP(LFP,[.1 1],1024);
% smootime=.006;
% 
% BreathingPower = tsd(Range(FilLFP),runmean(Data((FilLFP)).^2,ceil(smootime/median(diff(Range(FilLFP,'s'))))));
% BreathingPower_var = tsd(Range(BreathingPower) , movstd(Data(BreathingPower) , ceil(30/median(diff(Range(FilLFP,'s'))))));
% % figure
% % subplot(312)
% % plot(Range(BreathingPower_var,'s')/3600 , log10(Data(BreathingPower_var)))
% % 
% % figure
% % [Y,X]=hist(log10(Data(BreathingPower_var)),1000);
% % Y=Y/sum(Y);
% % plot(X,runmean(Y,90),'b','LineWidth',1)
% % xlabel('Breathing variability (a.u.)'), ylabel('PDF')
% % box off,% xlim([1.5 6])
% 
% % Sniff (nned B Middle spectro
% Phase=tsd(Range(FilLFP) , angle(hilbert(zscore(Data(FilLFP))))*180/pi+180);
% Phase_Above_350=thresholdIntervals(Phase,350,'Direction','Above');
% Sniff=ts((Stop(Phase_Above_350)+Start(Phase_Above_350))/2);
% 
% % % load(fullfile(datapath, 'B_Middle_Spectrum.mat'))
% % Stsd = tsd(Spectro{2}*1e4 , log10(Spectro{1}));
% % f = Spectro{3};
% % 
% % figure
% % subplot(3,4,[1 5])
% % [M,~,t]=AverageSpectrogram(Stsd,f,Sniff,50,250,0,.7,1);
% % imagesc(t/1E3,f,SmoothDec(M,2)), axis xy
% % % xlim([0 5]), ylim([30 100]), 
% % caxis([3.5 6.2]), ylabel('Frequency (Hz)')
% % title('Wake')
% % makepretty
% % 
% % subplot(349)
% % plot(t/1E3, nanmean(M) , 'k' , 'LineWidth' , 2), xlim([0 5])
% % xlabel('time (s)'), ylabel('OB gamma power (a.u.)')
% % makepretty
% % % ylim([4.1 4.8])
% % 
% % colormap jet
% % 
% % %
% % [M_sup,T] = PlotRipRaw(LFP, Range(Sniff)/1e4, 1e3, 0, 1,1);
% % figure
% % subplot(141)
% % plot(M_sup(:,1) , M_sup(:,2))
% % % hold on
% % % plot(M_sup_wake(:,1) , M_deep_wake(:,2))
% % % ylim([-600 600])
% % % 
% 
% % EMG
% EMGData=tsd(Range(FilLFP),runmean(Data((FilLFP)).^2,ceil(smootime/median(diff(Range(FilLFP,'s'))))));
% 
% subplot(6,6,[25 19 13 7 1])
% [Y,X] = hist(log10(Data(EMGData)),1000);
% a = area(X , runmean(Y,10)); a.FaceColor=[.8 .8 .8]; a.LineWidth=1.5; a.EdgeColor=[0 0 0];
% set(gca,'XDir','reverse'), camroll(270), box off
% v2=vline(3.5,'-r'); v2.LineWidth=3;
% xlabel('EMG power (log scale)'), xlim([2.5 6.5])
% 
% subplot(6,6,[2:6 8:12 14:18 20:24 26:30])
% X = log10(Data(SmoothGamma_intf)); Y = log10(Data(EMGData));
% plot(X(1:2e3:end) , Y(1:2e3:end) , '.k' , 'MarkerSize' , 3)
% axis square
% v1=vline(2.3,'-r'); v1.LineWidth=3;
% v2=hline(3.5,'-r'); v2.LineWidth=3;
% ylim([2.5 6.5])
% 
% % Accelero
% MovAcctsd=tsd(Range(MovAcctsd),runmean(Data((MovAcctsd)).^2,ceil(smootime/median(diff(Range(MovAcctsd,'s'))))));
% 
% figure
% plot(Range(Restrict(SmoothGamma, MovingEpoch),'s')/60 , Data(Restrict(SmoothGamma, MovingEpoch)))
% hold on
% plot(Range(Restrict(SmoothGamma, FreezeAccEpoch),'s')/60 , Data(Restrict(SmoothGamma, FreezeAccEpoch)))
% legend('Moving','Immobile')
% ylabel('Gamma power (a.u.)')
% xlabel('time (min)')
% title('Gamma power values across session')   
%     
%             
% MovAcctsd_NewRange = interp1(Range(NewMovAcctsd) , Data(NewMovAcctsd) , Range(SmoothGamma)); 
% MovAcctsd_NewRange_tsd = tsd(Range(SmoothGamma) , MovAcctsd_NewRange);
% 
% D1 = Data(MovAcctsd_NewRange_tsd);
% D2 = Data(SmoothGamma);
% D1_bis = D1(D1<1.3e7);
% D2_bis = D2(D1<1.3e7);
% 
% clf; bin=500;
% subplot(121)
% PlotCorrelations_BM(D1_bis(1:bin:end) , D2_bis(1:bin:end) , 5 , 1)
% ylim([50 150])
% xlabel('Accelerometer values (regular values)'); ylabel('Gamma values (regular values)')
% title('regular values')
% subplot(122)
% PlotCorrelations_BM(log10(D1_bis(1:bin:end)) , log10(D2_bis(1:bin:end)) , 5 , 1)
% xlabel('Accelerometer values (log values)'); ylabel('Gamma values (log values)')
% title('log values')
% 
% a=suptitle('OB gamma power = f(accelero)'); a.FontSize=20;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
end