function makePupilSummaryVideo(datapath , smoothing_win)
% -------------------------------------------------------------------------
%  Create a 4-panel MP4:
%     || Target frame | cum-plot Target || cum-plot Ref | Reference frame ||
%  for 20 random Target & 20 random Reference trials from the first half.
% -------------------------------------------------------------------------

% 1. Check sync between video and traces
% 2. Delete old moving xline
% 3. Leave only relevant traces in legend
% 4. Check what exactly do I calculate as pupil area. Try passing raw data

%% 0)  Parameters
nTrialsEach = 20;                % how many Target / Reference
fsVideo = 50;                % eye-cam FPS  (verified below)
soundDur = 2.1;               % s, stimulus period
rewardTail = 2;               % s, highlight after reward
baselineSec = 0.5;
outFile = fullfile(datapath,'video','pupil_summary_1.mp4');

% fadeAlpha = linspace(1,0.15,nTrialsEach);   % line-fade in cum-plots
colTar  = [1 0.4 0.4];   % red
colRef  = [0.25 0.6 1];  % blue
greyCol = 0.85*[1 1 1];  % stimulus patch

%% 1)  Load data
load(fullfile(datapath,'stim','trial_structure.mat'));
load(fullfile(datapath,'video','DLC_data.mat'),'pupil_area_007', 'time_eye');
timeVideo = time_eye / 1e4;
videoT0 = timeVideo(1);  % global time of the first video frame
pArea  = Data(pupil_area_007);

% eye-cam file
camFile = dir(fullfile(datapath,'video','*007*DLC*.mp4'));
assert(~isempty(camFile),'Eye-cam video not found');
vr = VideoReader(fullfile(camFile.folder,camFile.name));
assert(abs(vr.FrameRate-fsVideo)<0.2, ...
    'Video FPS %.2f ? %d',vr.FrameRate,fsVideo);
nFramesVid = ceil(vr.Duration*vr.FrameRate);

%% 2)  First-half trial pool
rng('shuffle')
good    = trial_structure.goodTrials;
halfN   = floor(numel(good)/2);
refPool = find(ismember(trial_id_struct.Reference,good(1:halfN)));
tarPool = find(ismember(trial_id_struct.Target   ,good(1:halfN)));
assert(numel(refPool)>=nTrialsEach && numel(tarPool)>=nTrialsEach, ...
    'Not enough good trials in first half');

selRef = refPool(randperm(numel(refPool),nTrialsEach));
selTar = tarPool(randperm(numel(tarPool),nTrialsEach));
fprintf('Ref trials: %s\n',mat2str(trial_id_struct.Reference(selRef).'))
fprintf('Tar trials: %s\n',mat2str(trial_id_struct.Target(selTar).'))

%% 3)  Prepare cumulative data  ------------------------------------------
stimOff   = 0.2;                                   % trial-to-stim offset
delayTar  = double(trial_structure.reward_onset.Target(selTar));
delayMean = mean(double(trial_structure.reward_onset.Target)); % for Ref

tRef        = double(Range(trial_structure.trial_onset.Reference,'s'));
tTar        = double(Range(trial_structure.trial_onset.Target ,'s'));

% window of interest here: [stim onset ; reward onset]
stimRefAbs = tRef(selRef) + stimOff;
stimTarAbs = tTar(selTar) + stimOff;
rewRefAbs  = tRef(selRef) + delayMean;
rewTarAbs  = tTar(selTar) + delayTar;

% include –baselineSec … +rewardTail around reward
winRef = arrayfun(@(s,e) find(timeVideo>=s-baselineSec & timeVideo<=e+rewardTail), ...
    stimRefAbs, rewRefAbs ,'uni',0);
winTar = arrayfun(@(s,e) find(timeVideo>=s-baselineSec & timeVideo<=e+rewardTail), ...
    stimTarAbs, rewTarAbs ,'uni',0);

% optional smoothing
if smoothing_win>1
    smoothFn = @(x) movmean(x,smoothing_win,'omitnan');
else
    smoothFn = @(x) x;
end

cellRef = cellfun(@(ix) smoothFn(pArea(ix).'), winRef,'uni',0);
cellTar = cellfun(@(ix) smoothFn(pArea(ix).'), winTar,'uni',0);

% subtract baseline (last baselineSec before stim = first baselineSec data)
nBase = round(baselineSec*fsVideo);
cellRef = cellfun(@(v) v - mean(v(1:nBase),'omitnan'),cellRef,'uni',0);
cellTar = cellfun(@(v) v - mean(v(1:nBase),'omitnan'),cellTar,'uni',0);

minLen  = min(cellfun(@numel,[cellRef;cellTar]));   % shortest trial length
cellRef = cellfun(@(v) v(1:minLen), cellRef,'uni',0);
cellTar = cellfun(@(v) v(1:minLen), cellTar,'uni',0);

matRef  = cell2mat(cellRef);
matTar  = cell2mat(cellTar);
tAxis   = (-baselineSec : 1/fsVideo : (-baselineSec)+(minLen-1)/fsVideo);

%% 4)  VideoWriter setup ---------------------------------------------------
vw = VideoWriter(outFile,'MPEG-4'); vw.FrameRate = fsVideo; open(vw);
fig = figure('Color','w','Units','pixels','Position',[100 100 840 720], 'visible', 'off');

axTar = subplot(2,2,1, 'Parent', fig); axis(axTar,'off');
axRef = subplot(2,2,2, 'Parent', fig); axis(axRef,'off');
axCum    = subplot(2,2,3, 'Parent', fig); hold(axCum,'on'); box(axCum,'off');
axTrial    = subplot(2,2,4, 'Parent', fig); hold(axTrial,'on'); box(axTrial,'off');
xlabel(axCum,'time (s)'); ylabel(axCum,'pupil area');
xlabel(axTrial,'time (s)'); ylabel(axTrial,'pupil area');

% grey stimulus patch (0-soundDur)
% patch(axCum,[0 soundDur soundDur 0],[-Inf Inf Inf -Inf], ...
%     greyCol,'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');

% hTar = shadedErrorBar_BM(tAxis,matTar,{'color',colTar,'LineWidth',2},1);
% hRef = shadedErrorBar_BM(tAxis,matRef,{'color',colRef,'LineWidth',2},1);
% legend(axCum,[hTar.mainLine, hRef.mainLine], ...
%     {'Target mean','Reference mean'},'Location','northwest');
% ylim(axCum,[min([matTar(:);matRef(:)]),max([matTar(:);matRef(:)])]+[-1e-3 1e-3])
% xlim(axCum,[tAxis(1) tAxis(end)])
% cursor = xline(axCum,0,'--k','LineWidth',1,'HandleVisibility','off');

%% ---------------- 5) MAIN LOOP  ----------------------------------------
cla(axCum); hold(axCum,'on'); box(axCum,'off')
hStim = patch(axCum, [0 soundDur soundDur 0], ...
    repmat([-1 -1 1 1]*1e9,1,1), ...
    greyCol,'EdgeColor','none','FaceAlpha',0.25, ...
    'HandleVisibility','off');

% mean ± SEM of *all* trials added so far
set(gcf, 'CurrentAxes', axCum)
hTar = shadedErrorBar_BM(tAxis, matTar, {'color',colTar,'LineWidth',2}, 1);
set(hTar.patch,'FaceAlpha',0.15,'EdgeAlpha',0.25);
hRef = shadedErrorBar_BM(tAxis, matRef, {'color',colRef,'LineWidth',2}, 1);
set(hRef.patch,'FaceAlpha',0.15,'EdgeAlpha',0.25);
hReward = xline([delayMean], 'r--','LineWidth',1.2);

legend([hStim, hRef.mainLine, hTar.mainLine, hReward], ...
    {'Stimulus','Reference','Target','Reward'},'Location','best')

ylim(axCum,[min([matTar(:);matRef(:)]) , ...
    max([matTar(:);matRef(:)])] + [-1e-3 1e-3]);
xlim(axCum,[tAxis(1) tAxis(end)]);
xlabel(axCum,'time (s)'); ylabel(axCum,'pupil area');
title(['number of trials = ' num2str(nTrialsEach)])

for k = 1:nTrialsEach
    disp(['Working on trial #' num2str(k)]);
    cla(axTrial); hold(axTrial,'on'); box(axTrial,'off')
    hStim = patch(axTrial, [0 soundDur soundDur 0], ...
        repmat([-1 -1 1 1]*1e9,1,1), ...
        greyCol,'EdgeColor','none','FaceAlpha',0.25, ...
        'HandleVisibility','off');
    
    curTar = matTar(k, :);
    curRef = matRef(k, :);
    set(gcf, 'CurrentAxes', axTrial)
    hTar = plot(axTrial,tAxis,curTar,'Color',colTar,'LineWidth',1);
    hRef = plot(axTrial,tAxis,curRef,'Color',colRef,'LineWidth',1);
    hReward = xline([delayMean], 'r--','LineWidth',1.2);
    legend([hStim, hRef, hTar, hReward], ...
        {'Stimulus','Reference','Target','Reward'},'Location','best')
        
    ylim(axTrial,[min([matTar(:);matRef(:)]) , ...
        max([matTar(:);matRef(:)])] + [-1e-3 1e-3]);
    xlim(axTrial,[tAxis(1) tAxis(end)]);
    xlabel(axTrial,'time (s)'); ylabel(axTrial,'pupil area');
    title(axTrial, ['Current Ref trial # ' num2str(trial_id_struct.Reference(selRef(k))) '. ' 'Current Tar trial # ' num2str(trial_id_struct.Target(selTar(k)))])
    
    % ------------------------------------------------ cumulative panel ---
    cursor = xline(axTrial,0,'--k','LineWidth',1,'HandleVisibility','off');
    
    % --------------- frame-by-frame playback for THIS trial -------------
    for r = 1:minLen
        tRel = tAxis(r);                           % -0.5 … (reward)
        absTar = stimTarAbs(k) + tRel;
        absRef = stimRefAbs(k) + tRel;
        
        frmTar = round((absTar - videoT0)*fsVideo)+1;
        frmRef = round((absRef - videoT0)*fsVideo)+1;
        
        % -------- TARGET panel (left) -----------------------------------
        if frmTar<=nFramesVid
            vr.CurrentTime = (frmTar-1)/vr.FrameRate;
            imshow(readFrame(vr),'Parent',axTar); axis(axTar,'off')
        end
        if      0<=tRel && tRel<=soundDur
            title(axTar,'Target Stimulus','Color','g');
            rectangle(axTar,'Position',[0 0 1 1],'EdgeColor','g','LineWidth',3);
        elseif  tRel>delayTar(k) && tRel<=delayTar(k)+rewardTail
            title(axTar,'Reward','Color','r');
            rectangle(axTar,'Position',[0 0 1 1],'EdgeColor','r','LineWidth',3);
        else,   title(axTar,'')
        end
        
        % -------- REFERENCE panel (right) -------------------------------
        if frmRef<=nFramesVid
            vr.CurrentTime = (frmRef-1)/vr.FrameRate;
            imshow(readFrame(vr),'Parent',axRef); axis(axRef,'off')
        end
        if      0<=tRel && tRel<=soundDur
            title(axRef,'Reference Stimulus','Color','g');
            rectangle(axRef,'Position',[0 0 1 1],'EdgeColor','g','LineWidth',3);
        elseif  tRel>delayMean && tRel<=delayMean+rewardTail
            title(axRef,'Sham','Color','r');
            rectangle(axRef,'Position',[0 0 1 1],'EdgeColor','r','LineWidth',3);
        else,   title(axRef,'')
        end
        
        % -------- move cursor & write frame -----------------------------
        set(cursor,'Value',tRel);
        drawnow limitrate nocallbacks
        writeVideo(vw,getframe(fig));
    end
    
    % -------- 1-s pause before next trial pair --------------------------
    pauseFrame = getframe(fig);
    for p = 1:fsVideo, writeVideo(vw,pauseFrame); end
end
close(vw)
end
