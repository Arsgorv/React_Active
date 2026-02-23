function behaviour_analysis(datapath, opts)
% behaviour_analysis (React_Active)
% Per-session: metrics + PSTH figures for each marker for several windows and trial subsets.
%
% Figures:
% For each bodypart/physio marker:
%   windows: stimon_to_arrival, stimoff_to_arrival, arrival_m500_to_arrival
%   trial subsets: regular, nosound, nomotor
%   comparison: Target vs Reference (within subset)
% Metrics: AshmanD, AUROC, CohensD, meanDiff
% Plot: PSTH + hist + meanDiff + CohenD + ROC

if nargin < 2, opts = struct(); end
if ~isfield(opts,'fs_video'), opts.fs_video = 50; end
if ~isfield(opts,'smoothing_win_s'), opts.smoothing_win_s = 0.50; end
if ~isfield(opts,'baseline_pre_s'), opts.baseline_pre_s = 0.20; end
if ~isfield(opts,'epoch_post_reward_s'), opts.epoch_post_reward_s = 2.0; end
if ~isfield(opts,'min_trials_per_group'), opts.min_trials_per_group = 3; end

fs = opts.fs_video;

% output folders ----------------
outBase = fullfile(datapath,'analysis','behaviour');
if ~exist(outBase,'dir'), mkdir(outBase); end
[~,sessname] = fileparts(datapath);

%% Load Baphy
Baphy = load_baphy(datapath);

%% Load DLC
D = load(fullfile(datapath,'video','DLC_data.mat'));
time_video_s = ra_time_to_seconds(double(D.time_face(:))); % robust seconds

%% Load physio
P = struct();
physio_file = fullfile(datapath,'ephys','physio_for_behaviour.mat');
if exist(physio_file,'file')
    P = load(physio_file);
end

%% trial timing
n = Baphy.n_trials;
trial_start_abs_s = Baphy.trial.abs_trialstart_s(:);
stim_on_abs_s     = Baphy.trial.abs_stim_start_s(:);
stim_off_abs_s    = Baphy.trial.abs_stim_stop_s(:);

arrival_abs_s = nan(n,1);
if isfield(Baphy.trial,'spout_arrival_abs_s')
    arrival_abs_s = Baphy.trial.spout_arrival_abs_s(:);
end

stim_on_rel_s  = stim_on_abs_s  - trial_start_abs_s;
stim_off_rel_s = stim_off_abs_s - trial_start_abs_s;
arrival_rel_s  = arrival_abs_s  - trial_start_abs_s;

% --- Fix for NoMotor: arrival is undefined -> use median arrival across trials that have it
arrival_imputed = ~isfinite(arrival_rel_s);
arr_med_all = nanmedian(arrival_rel_s(~arrival_imputed));
if isfinite(arr_med_all)
    arrival_rel_s(arrival_imputed) = arr_med_all;
end

%% reward relative for epoch sizing
reward_rel_s_pertrial = nan(n,1);
if isfield(Baphy.trial,'abs_reward_s')
    reward_rel_s_pertrial = Baphy.trial.abs_reward_s(:) - trial_start_abs_s;
elseif isfield(Baphy.trial,'abs_target_s')
    reward_rel_s_pertrial = Baphy.trial.abs_target_s(:) - trial_start_abs_s;
end
reward_rel_s = nanmedian(reward_rel_s_pertrial(isfinite(reward_rel_s_pertrial)));
if ~isfinite(reward_rel_s)
    reward_rel_s = nanmedian(stim_off_rel_s(isfinite(stim_off_rel_s))) + 2.0;
end

trial_dur_s = Baphy.trial.trial_dur_s;
n_samples = round(trial_dur_s*fs)+1;
tvec = (0:n_samples-1)/fs;
baseline_samples = max(1, round(opts.baseline_pre_s*fs));
smoothSamples = max(1, round(opts.smoothing_win_s*fs));

%% define windows
winNames = {'stimon_to_stimoff','stimon_to_arrival','stimoff_to_arrival','arrival_m500_to_arrival'};

%% trial masks
mskGood = false(n,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'goodTrials')
    mskGood(Baphy.idx.goodTrials(:)) = true;
else
    mskGood = isfinite(trial_start_abs_s);
end

mskTar = false(n,1); mskRef = false(n,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'Target')
    mskTar(Baphy.idx.Target(:)) = true;
    mskRef(Baphy.idx.Reference(:)) = true;
elseif isfield(Baphy.trial,'type')
    mskTar(strcmpi(Baphy.trial.type,'Target')) = true;
    mskRef(strcmpi(Baphy.trial.type,'Reference')) = true;
end

mskNS = false(n,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'NoSound')
    mskNS(Baphy.idx.NoSound(:)) = true;
elseif isfield(Baphy.trial,'is_nosound')
    mskNS = logical(Baphy.trial.is_nosound(:));
end

mskNM = false(n,1);
if isfield(Baphy,'idx') && isfield(Baphy.idx,'NomotorExclusive')
    mskNM(Baphy.idx.NomotorExclusive(:)) = true;
elseif isfield(Baphy.trial,'is_nomotor')
    mskNM = logical(Baphy.trial.is_nomotor(:));
end

mskReg = ~(mskNS | mskNM);

%% subset struct (Target vs Ref within each)
subsets = struct();
subsets(1).name = 'regular';
subsets(1).idxA = find(mskTar & mskGood & mskReg);
subsets(1).idxB = find(mskRef & mskGood & mskReg);

subsets(2).name = 'nosound';
subsets(2).idxA = find(mskTar & mskGood & mskNS);
subsets(2).idxB = find(mskRef & mskGood & mskNS);

subsets(3).name = 'nomotor';
subsets(3).idxA = find(mskTar & mskGood & mskNM);
subsets(3).idxB = find(mskRef & mskGood & mskNM);

%% Markers
bodyparts = {
    'Pupil', {'pupil_area_004','pupil_center_004','pupil_center_004_mvt'};
    'Eye',   {'eye_area_004'};
    'Nostril',{'nostril_area','nostril_center','nostril_center_mvt'};
    'Nose',  {'nose_area','nose_center','nose_center_mvt'};
    'Cheek', {'cheek_center','cheek_center_mvt'};
    'Ear',   {'ear_center','ear_center_mvt'};
    'Jaw',   {'jaw_center','jaw_center_mvt'};
    'Tongue',{'tongue_center','tongue_center_mvt'};
    'Pupil-EyeCam', {'pupil_area_007','pupil_center_007','pupil_center_007_mvt'};
    'Eye-EyeCam',   {'eye_area_007'};
    'Spout', {'spout_likelihood'};

    % Physio (from physio_for_behaviour.mat)
    'OB-delta-power', {'OB_delta_rs','OB_delta_log10_rs','OB_delta_z_rs'};
    'OB-gamma-power', {'OB_gamma_rs','OB_gamma_log10_rs','OB_gamma_z_rs'};
    'OB-gamma-fast',  {'OB_gamma_fast_rs','OB_gamma_fast_log10_rs'};
    'OB-delta-fast',  {'OB_delta_fast_rs','OB_delta_fast_log10_rs'};

    'Respiration', {'RespPower_rs','RespPower_log10_rs','RespPhase_rs'};
    'EMG', {'EMGPower_rs','EMGPower_log10_rs','EMGPhase_rs'};
    'Accelerometer', {'AccPower_rs','AccPower_log10_rs'};
    'Heart', {'HeartPower_rs','HeartPower_log10_rs','HeartPhase_rs'};
};

% Add imputed-arrival summary fields; keep ROC arrays in MAT (CSV drops them later)
varNames = {'bodypart','marker','subset','window', ...
    'nTar','nRef', ...
    'nArrivalImputed','fracArrivalImputed', ...
    'AshmanD','AUROC','CohensD','meanDiff', ...
    'meanRef','meanTar','stdRef','stdTar', ...
    'rocFPR','rocTPR'};

%% Colours
nbp = size(bodyparts,1);
baseColors = lines(nbp);
maxPerPart = 5;
shadeLevels = linspace(0,0.55,maxPerPart+1);
shadeLevels(1) = [];
lighten = @(c,a) c + (1-c).*a;

%% Main loop
for b = 1:size(bodyparts,1)
    bp = bodyparts{b,1};
    markers = bodyparts{b,2};
    outDir = fullfile(outBase, bp);
    if ~exist(outDir,'dir'), mkdir(outDir); end

    for ss = 1:numel(subsets)
        subsetName = subsets(ss).name;
        idxA = subsets(ss).idxA;
        idxB = subsets(ss).idxB;

        idxA = idxA(isfinite(trial_start_abs_s(idxA)));
        idxB = idxB(isfinite(trial_start_abs_s(idxB)));

        for w = 1:numel(winNames)
            wName = winNames{w};

            n_markers = numel(markers);
            colsTotal = 6; % PSTH | dist+kde | meanDiff | CohenD | AUROC | ROC

            f = figure('Color','w','Visible','on');
            set(f,'Position',[50 50 1900 max(420, 220*n_markers)]);
            sgtitle(sprintf('%s | %s | %s | sw=%.3fs | Tar=%d Ref=%d', ...
                bp, subsetName, wName, opts.smoothing_win_s, numel(idxA), numel(idxB)), ...
                'Interpreter','none');

            rows = cell(0, numel(varNames));
            allMet = struct('marker',{},'metrics',{});

            % fixed arrival anchor (avoid identity leakage)
            arrA = arrival_rel_s(idxA);
            arrB = arrival_rel_s(idxB);
            medA = nanmedian(arrA(isfinite(arrA)));
            medB = nanmedian(arrB(isfinite(arrB)));

            arr_anchor = nan;
            if isfinite(medA) && isfinite(medB)
                arr_anchor = min(medA, medB);
            elseif isfinite(medA)
                arr_anchor = medA;
            elseif isfinite(medB)
                arr_anchor = medB;
            end

            for m = 1:n_markers
                baseCol = baseColors(b,:);
                col = baseCol;
                if m > 1
                    a = shadeLevels(min(m-1, numel(shadeLevels)));
                    col = lighten(baseCol, a);
                end
                colRef = col;
                colTar = max(0, min(1, col*0.5));
                ColsRT = {colRef, colTar};
                LegRT  = {'Ref','Tar'};

                mk = markers{m};

                [t_s, y] = get_marker_xy(mk, D, P);
                if isempty(t_s)
                    continue
                end

                goodA = idxA(isfinite(trial_start_abs_s(idxA)));
                goodB = idxB(isfinite(trial_start_abs_s(idxB)));

                onA = trial_start_abs_s(goodA);
                onB = trial_start_abs_s(goodB);

                matA = extract_matrix(t_s, y, onA, tvec(:), baseline_samples, smoothSamples);
                matB = extract_matrix(t_s, y, onB, tvec(:), baseline_samples, smoothSamples);

                % per-trial window means
                xA = nan(size(matA,1),1);
                xB = nan(size(matB,1),1);

                for k = 1:size(matA,1)
                    tr = goodA(k);
                    wsec = window_for_trial(wName, tr, stim_on_rel_s, stim_off_rel_s, arrival_rel_s, arr_anchor);
                    xA(k) = mean_in_window(matA(k,:), wsec, fs, n_samples);
                end
                for k = 1:size(matB,1)
                    tr = goodB(k);
                    wsec = window_for_trial(wName, tr, stim_on_rel_s, stim_off_rel_s, arrival_rel_s, arr_anchor);
                    xB(k) = mean_in_window(matB(k,:), wsec, fs, n_samples);
                end

                dRef = xB(isfinite(xB));
                dTar = xA(isfinite(xA));
                nRef = numel(dRef);
                nTar = numel(dTar);

                doStats = (nRef >= opts.min_trials_per_group) && (nTar >= opts.min_trials_per_group);

                if doStats
                    met = effect_metrics(xA, xB);
                else
                    met = struct('AshmanD',nan,'AUROC',nan,'CohensD',nan,'meanDiff',nan, ...
                        'meanRef',nan,'meanTar',nan,'stdRef',nan,'stdTar',nan, ...
                        'rocFPR',[],'rocTPR',[]);
                end

                % arrival imputation summary for THIS subset (same for all markers, but we store per-row)
                mask_subset = false(n,1);
                mask_subset([goodA(:); goodB(:)]) = true;
                nImp = sum(arrival_imputed & mask_subset);
                nTot = sum(mask_subset);
                fracImp = nImp / max(nTot,1);

                % PSTH stats
                meanA = nanmean(matA,1);
                meanB = nanmean(matB,1);
                semA  = nanstd(matA,[],1) / max(1, sqrt(size(matA,1)));
                semB  = nanstd(matB,[],1) / max(1, sqrt(size(matB,1)));

                lo = min([meanA-semA, meanB-semB]); hi = max([meanA+semA, meanB+semB]);
                if ~isfinite(lo), lo = -1; end
                if ~isfinite(hi), hi =  1; end
                if lo == hi, lo = lo-1; hi = hi+1; end
                pad = 0.05*(hi-lo);
                yl = [lo-pad hi+pad];

                pooled = [goodA; goodB];
                stim_on_med  = nanmedian(stim_on_rel_s(pooled));
                stim_off_med = nanmedian(stim_off_rel_s(pooled));
                arr_med      = nanmedian(arrival_rel_s(pooled));
                rew_med      = nanmedian(reward_rel_s_pertrial(pooled));

                w0 = nan(numel(pooled),1); w1 = nan(numel(pooled),1);
                for kk = 1:numel(pooled)
                    tr = pooled(kk);
                    wsec = window_for_trial(wName, tr, stim_on_rel_s, stim_off_rel_s, arrival_rel_s, arr_anchor);
                    w0(kk) = wsec(1); w1(kk) = wsec(2);
                end
                w_med = [nanmedian(w0) nanmedian(w1)];

                % ----- plotting row -----
                row0 = (m-1)*colsTotal;

                hasTar = ~isempty(matA) && size(matA,1) >= 1;
                hasRef = ~isempty(matB) && size(matB,1) >= 1;

                %% 1) PSTH
                subplot(n_markers, colsTotal, row0+1); hold on

                hStim = [];
                if isfinite(stim_on_med) && isfinite(stim_off_med) && stim_off_med > stim_on_med
                    hStim = patch([stim_on_med stim_off_med stim_off_med stim_on_med], [yl(1) yl(1) yl(2) yl(2)], ...
                        [0.7 0.7 0.7], 'EdgeColor','none','FaceAlpha',0.18);
                end

                hSpout = [];
                hReward = [];
                if isfinite(arr_med), hSpout = plot([arr_med arr_med], yl, 'k--','LineWidth',1.2); end
                if isfinite(rew_med), hReward = plot([rew_med rew_med], yl, 'r--','LineWidth',1.2); end

                hAnt = [];
                if all(isfinite(w_med)) && w_med(2) > w_med(1)
                    hAnt = patch([w_med(1) w_med(2) w_med(2) w_med(1)], yl([1 1 2 2]), ...
                        [0 1 0], 'EdgeColor','none','FaceAlpha',0.08);
                end

                hTar = [];
                if hasTar
                    hTar = shadedErrorBar_BM(tvec, matA, {'color',colTar,'LineWidth',2}, 1);
                    set(hTar.patch,'FaceAlpha',0.15,'EdgeAlpha',0.25);
                end

                hRef = [];
                if hasRef
                    hRef = shadedErrorBar_BM(tvec, matB, {'color',colRef,'LineWidth',2}, 1);
                    set(hRef.patch,'FaceAlpha',0.15,'EdgeAlpha',0.25);
                end

                ylabel(mk,'Interpreter','none');
                if m==n_markers, xlabel('Time from trial onset (s)'); end
                xlim([tvec(1) tvec(end)]); ylim(yl);
                title(sprintf('%s (Tar=%d Ref=%d)', mk, nTar, nRef), 'Interpreter','none');
                box off

                if m==1
                    hh = []; ll = {};
                    if ~isempty(hStim), hh(end+1)=hStim; ll{end+1}='Stimulus'; end
                    if ~isempty(hAnt),  hh(end+1)=hAnt;  ll{end+1}='Anticip window'; end
                    if ~isempty(hRef),  hh(end+1)=hRef.mainLine; ll{end+1}='Reference'; end
                    if ~isempty(hTar),  hh(end+1)=hTar.mainLine; ll{end+1}='Target'; end
                    if ~isempty(hSpout),  hh(end+1)=hSpout;  ll{end+1}='Spout'; end
                    if ~isempty(hReward), hh(end+1)=hReward; ll{end+1}='Reward'; end
                    if ~isempty(hh)
                        lgd = legend(hh,ll,'Location','westoutside');
                        drawnow
                        set(lgd,'Units','normalized');
                        pL = lgd.Position;
                        lgd.Position = [pL(1)-0.08  pL(2)  pL(3)  pL(4)];
                    end
                end

                %% 2) Histogram + KDE
                subplot(n_markers, colsTotal, row0+2); hold on
                if nTar > 0
                    histogram(dTar,30,'Normalization','pdf','FaceColor',colTar,'FaceAlpha',0.35,'EdgeColor','none');
                end
                if nRef > 0
                    histogram(dRef,30,'Normalization','pdf','FaceColor',colRef,'FaceAlpha',0.35,'EdgeColor','none');
                end
                if nTar >= 7 && numel(unique(dTar)) > 1
                    [f1,xi1] = ksdensity(dTar); plot(xi1,f1,'Color',colTar,'LineWidth',1.2);
                end
                if nRef >= 7 && numel(unique(dRef)) > 1
                    [f2,xi2] = ksdensity(dRef); plot(xi2,f2,'Color',colRef,'LineWidth',1.2);
                end
                title(sprintf('AshmanD=%.2f', met.AshmanD), 'Interpreter','none');
                box off; set(gca,'YAxisLocation','right');

                %% 3) meanDiff
                subplot(n_markers, colsTotal, row0+3); cla; hold on
                if doStats
                    MakeSpreadAndBoxPlot3_SB({dRef, dTar}, ColsRT, [1 2], LegRT, ...
                        'newfig',0,'paired',0,'optiontest','ranksum','showpoints',1,'ShowSigstar','sig');
                    title(sprintf('meanDiff (Tar-Ref) = %.3g', met.meanDiff), 'Interpreter','none');
                else
                    axis off
                    title(sprintf('meanDiff skipped (nRef=%d nTar=%d)', nRef, nTar), 'Interpreter','none');
                end
                box off

                %% 4) CohenD
                subplot(n_markers, colsTotal, row0+4); cla; hold on
                if doStats
                    MakeSpreadAndBoxPlot3_SB({dRef, dTar}, ColsRT, [1 2], LegRT, ...
                        'newfig',0,'paired',0,'optiontest','ranksum','showpoints',1,'ShowSigstar','sig');
                    title(sprintf('Cohen''s d = %.2f', met.CohensD), 'Interpreter','none');
                else
                    axis off
                    title(sprintf('Cohen''s d skipped (nRef=%d nTar=%d)', nRef, nTar), 'Interpreter','none');
                end
                box off

                %% 5) AUROC
                subplot(n_markers, colsTotal, row0+5); cla; hold on
                if doStats
                    MakeSpreadAndBoxPlot3_SB({dRef, dTar}, ColsRT, [1 2], LegRT, ...
                        'newfig',0,'paired',0,'optiontest','ranksum','showpoints',1,'ShowSigstar','sig');
                    title(sprintf('AUROC = %.2f', met.AUROC), 'Interpreter','none');
                else
                    axis off
                    title(sprintf('AUROC skipped (nRef=%d nTar=%d)', nRef, nTar), 'Interpreter','none');
                end
                box off

                %% 6) ROC
                subplot(n_markers, colsTotal, row0+6); cla; hold on
                if doStats && ~isempty(met.rocFPR)
                    plot(met.rocFPR, met.rocTPR, 'k', 'LineWidth',1.5);
                    plot([0 1],[0 1],'k:');
                    xlim([0 1]); ylim([0 1]); axis square
                    title(sprintf('AUC=%.2f', met.AUROC), 'Interpreter','none');
                else
                    plot([0 1],[0 1],'k:'); xlim([0 1]); ylim([0 1]); axis square
                    title(sprintf('ROC skipped'), 'Interpreter','none');
                end
                xlabel('FPR'); ylabel('TPR'); box off

                %% store row
                rows(end+1,:) = {bp, mk, subsetName, wName, ...
                    nTar, nRef, ...
                    nImp, fracImp, ...
                    met.AshmanD, met.AUROC, met.CohensD, met.meanDiff, ...
                    met.meanRef, met.meanTar, met.stdRef, met.stdTar, ...
                    {met.rocFPR}, {met.rocTPR}};

                allMet(end+1).marker = mk;
                allMet(end).metrics = met;
            end

            if isempty(rows)
                close(f);
                continue
            end

            T = cell2table(rows, 'VariableNames', varNames);

            % save figure + per-(bp,subset,window) metrics
            baseName = sprintf('%s_%s_%s_%s_sw%.3f', sessname, bp, subsetName, wName, opts.smoothing_win_s);

            saveas(f, fullfile(outDir, [baseName '.png']));
            saveas(f, fullfile(outDir, [baseName '.svg']));

            close(f);

            % Save full table (keeps ROC arrays) in this condition-specific file
            S = struct();
            S.table = T;
            S.session = sessname;
            S.bodypart = bp;
            S.subset = subsetName;
            S.window = wName;
            S.allMet = allMet;
            save(fullfile(outDir, [baseName '_metrics.mat']), '-struct','S');

            % Also export CSV without ROC arrays (CSV-safe)
            T_full = T;
            vars_to_drop = intersect({'rocFPR','rocTPR'}, T_full.Properties.VariableNames);
            T_csv = T_full;
            if ~isempty(vars_to_drop)
                T_csv = removevars(T_csv, vars_to_drop);
            end
            writetable(T_csv, fullfile(outDir, [baseName '_metrics.csv']));

            % Session-level “latest” pointers (for across-session collection):
            % keep one full MAT and one CSV per session per bodypart by overwriting.
            save(fullfile(outDir, sprintf('%s_metrics_full.mat', sessname)), 'T_full');
            writetable(T_csv, fullfile(outDir, sprintf('%s_metrics.csv', sessname)));
        end
    end
end

end

function t_s = ra_time_to_seconds(t)
t = t(:);
if numel(t) < 2
    t_s = t;
    return
end
dt = nanmedian(diff(t));
if ~isfinite(dt)
    t_s = t;
    return
end
if dt > 1
    t_s = t / 1e4;
else
    t_s = t;
end
end


%% outdated: Plot distribution of first licks
% lick_ref = trial_structure.lick_onset.Reference(find(~isnan(trial_structure.lick_onset.Reference)));
% lick_tar = trial_structure.lick_onset.Target(find(~isnan(trial_structure.lick_onset.Target)));
% figure('Name','First-lick latency'); clf
% MakeSpreadAndBoxPlot3_SB({lick_ref, lick_tar},{[1 0 0],[0 0 1]},[2,4],{'Reference','Target'});
% ylabel('Latency from trial onset (s)')
% title('First-lick latency distribution')

%% outdated: pieces of BM scripts physio
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






