function decoding_trial_identity(datapath)
%% Initialize
load(fullfile(datapath, 'stim', 'trial_structure.mat'))
load(fullfile(datapath, 'video', 'DLC_data.mat'))

%% Extract marker trajectories for each trial in decoding window
% bodyparts = {
%     'Pupil', {'pupil_area_004','pupil_center_004','pupil_center_004_mvt','pupil_center_004_velocity'};
%     'Eye',   {'eye_area_004'};
%     'Nostril',{'nostril_area','nostril_center','nostril_center_mvt','nostril_center_velocity'};
%     'Nose',  {'nose_area','nose_center','nose_center_mvt','nose_center_velocity'};
%     'Cheek', {'cheek_center','cheek_center_mvt','cheek_center_velocity'};
%     'Ear',   {'ear_center','ear_center_mvt','ear_center_velocity'};
%     'Jaw',   {'jaw_center','jaw_center_mvt','jaw_center_velocity'};
%     'Tongue',{'tongue_center','tongue_center_mvt','tongue_center_velocity'};
%     'Pupil (EyeCam)', {'pupil_area_007','pupil_center_007','pupil_center_007_mvt','pupil_center_007_velocity'};
%     'Eye (EyeCam)',   {'eye_area_007'};
% };

marker_names = {
    'pupil_area_004','pupil_center_004','pupil_center_004_mvt','pupil_center_004_velocity',...
    'eye_area_004','nostril_area','nostril_center','nostril_center_mvt','nostril_center_velocity',...
    'nose_area','nose_center','nose_center_mvt','nose_center_velocity',...
    'cheek_center','cheek_center_mvt','cheek_center_velocity',...
    'ear_center','ear_center_mvt','ear_center_velocity',...
    'jaw_center','jaw_center_mvt','jaw_center_velocity',...
    'tongue_center','tongue_center_mvt','tongue_center_velocity',...
    'pupil_area_007','pupil_center_007','pupil_center_007_mvt','pupil_center_007_velocity',...
    'eye_area_007'
    };

%% Define decoding window
n_trials = min(length(trial_structure.trial_onset.Target), length(trial_structure.trial_onset.Reference));

stim_tar = trial_structure.stim_onset.Target(1:n_trials);
reward_tar = trial_structure.reward_onset.Target(1:n_trials);
stim_ref = trial_structure.stim_onset.Reference(1:n_trials);

trial_tar = Range(trial_structure.trial_onset.Target, 's');
trial_ref = Range(trial_structure.trial_onset.Reference, 's');

win_starts_tar = stim_tar + 2.1;   
win_ends_tar   = reward_tar;       
win_starts_ref = stim_ref + 2.1;
win_ends_ref   = reward_tar;

baseline_pre = 0.5; % seconds

%% Build a feature matrix
n_feat_per_marker = 8;
feature_tar = nan(n_trials, n_markers * n_feat_per_marker);
feature_ref = nan(n_trials, n_markers * n_feat_per_marker);
for m = 1:n_markers
    marker = eval(marker_names{m});
    t_marker = Range(marker, 's');
    y_marker = Data(marker);
    for i = 1:n_trials
        % TARGET
        idx_base = t_marker >= (trial_tar(i)+win_starts_tar(i)-baseline_pre) & t_marker < (trial_tar(i)+win_starts_tar(i));
        idx_win  = t_marker >= (trial_tar(i)+win_starts_tar(i)) & t_marker < (trial_tar(i)+win_ends_tar(i));
        window_vals = y_marker(idx_win) - nanmean(y_marker(idx_base));
        if isempty(window_vals), continue; end
        
        % Features
        f_idx = (m-1)*n_feat_per_marker;
        feature_tar(i, f_idx+1) = nanmean(window_vals);    % Mean
        feature_tar(i, f_idx+2) = nanstd(window_vals);     % Std
        feature_tar(i, f_idx+3) = nanmax(window_vals);     % Max
        feature_tar(i, f_idx+4) = nanmin(window_vals);     % Min
        feature_tar(i, f_idx+5) = (window_vals(end)-window_vals(1))/(win_ends_tar(i)-win_starts_tar(i)); % Slope
        feature_tar(i, f_idx+6) = window_vals(1);          % Start value
        feature_tar(i, f_idx+7) = window_vals(end);        % End value
        feature_tar(i, f_idx+8) = nansum(window_vals);   % AUC (optional)

        % repeat for feature_ref
        idx_base = t_marker >= (trial_ref(i)+win_starts_ref(i)-baseline_pre) & t_marker < (trial_ref(i)+win_starts_ref(i));
        idx_win  = t_marker >= (trial_ref(i)+win_starts_ref(i)) & t_marker < (trial_ref(i)+win_ends_ref(i));
        window_vals = y_marker(idx_win) - nanmean(y_marker(idx_base));
        if isempty(window_vals), continue; end
        feature_ref(i, f_idx+1) = nanmean(window_vals);
        feature_ref(i, f_idx+2) = nanstd(window_vals);
        feature_ref(i, f_idx+3) = nanmax(window_vals);
        feature_ref(i, f_idx+4) = nanmin(window_vals);
        feature_ref(i, f_idx+5) = (window_vals(end)-window_vals(1))/(win_ends_ref(i)-win_starts_ref(i));
        feature_ref(i, f_idx+6) = window_vals(1);
        feature_ref(i, f_idx+7) = window_vals(end);
        feature_ref(i, f_idx+8) = nansum(window_vals);
    end
end

% Remove trials with any NaN
valid_tar = all(~isnan(feature_tar),2);
valid_ref = all(~isnan(feature_ref),2);
feature_tar = feature_tar(valid_tar,:);
feature_ref = feature_ref(valid_ref,:);

%% Unsupervised Dimensionality reduction. Either PCA or t-SNE or UMAP
X = [feature_ref; feature_tar];
y = [zeros(size(feature_ref,1),1); ones(size(feature_tar,1),1)];

[~, frvals, frvecs, trnsfrmd] = pca(X, 'raw', 95, 1);

figure;
gscatter(trnsfrmd(:,1), trnsfrmd(:,2), y, 'kb', 'ox');
xlabel('PC1'); ylabel('PC2');
legend({'Reference','Target'}); title('PCA of Marker Features'); grid on;


Y = tsne(X, 'NumDimensions',2, 'Standardize',true);
figure;
gscatter(Y(:,1), Y(:,2), y, 'kb', 'ox');
xlabel('t-SNE1'); ylabel('t-SNE2');
legend({'Reference','Target'}); title('t-SNE of Marker Features'); grid on;

%% Supervised Decoding
% LDA ; Logistic regression ; random forest ; support vector machine ;
% cross-validation

%% Find the most informative marker
%{ 
Look at classifier feature weights (for LDA/logistic regression)
Or, use feature importance from a Random Forest
Ablation: re-run decoding with each marker left out and see if performance drops
Alternatively, try univariate ROC/AUC: For each marker alone, how well does it distinguish Target vs Reference?
%}

%% Visualize
%{
Plot mean value of the most important marker(s) per trial type (as you’ve done).
Plot PCA/UMAP projections colored by trial type.
Report cross-validated decoding accuracy.
Plot classifier feature weights for each marker. 
%}





end