function MetITPC = itpc_scan_plot(LFP, ev, G, opts, tag, phaseBank)

if ~isfield(opts,'min_events'),         opts.min_events = 10; end
if ~isfield(opts,'itpc_dt_s'),          opts.itpc_dt_s = 1/1024; end
if ~isfield(opts,'itpc_twin_s'),        opts.itpc_twin_s = [-0.5 7]; end
if ~isfield(opts,'itpc_f_centers'),     opts.itpc_f_centers = 2:2:120; end
if ~isfield(opts,'itpc_bw_hz'),         opts.itpc_bw_hz = 2; end

if ~isfield(opts,'ob_band_delta'),      opts.ob_band_delta = [0.5 4]; end
if ~isfield(opts,'ob_band_theta'),      opts.ob_band_theta = [4 6]; end

if ~isfield(opts,'itpc_plot_f_xlim'),   opts.itpc_plot_f_xlim = []; end
if ~isfield(opts,'itpc_band_hz'),       opts.itpc_band_hz  = opts.ob_band_delta; end
if ~isfield(opts,'itpc_band2_hz'),      opts.itpc_band2_hz = []; end

% null correction options
if ~isfield(opts,'itpc_use_null'),      opts.itpc_use_null = true; end
if ~isfield(opts,'itpc_null_n'),        opts.itpc_null_n = 200; end

fCenters = opts.itpc_f_centers(:)';   % 1 x nF
bw       = opts.itpc_bw_hz;
twin_s   = opts.itpc_twin_s;

if isempty(opts.itpc_plot_f_xlim)
    fShowMask = true(size(fCenters));
else
    fShowMask = (fCenters >= opts.itpc_plot_f_xlim(1)) & (fCenters <= opts.itpc_plot_f_xlim(2));
end
if ~any(fShowMask), fShowMask = true(size(fCenters)); end

dt = opts.itpc_dt_s;
t_rel = (twin_s(1):dt:twin_s(2))';

evList = {ev.all, ev.tar, ev.ref, ev.nosound, ev.nomotor};
colNames = {'all','tar','ref','nosound','nomotor'};

% Build phase bank if needed
if nargin < 6 || isempty(phaseBank)
    phaseBank = prepare_phase_bank(LFP, opts);
end
if ~isequal(phaseBank.fCenters, fCenters) || phaseBank.bw ~= bw
    phaseBank = prepare_phase_bank(LFP, opts);
end

b1Mask = (fCenters >= opts.itpc_band_hz(1)) & (fCenters <= opts.itpc_band_hz(2));
b2Mask = false(size(fCenters));
if ~isempty(opts.itpc_band2_hz)
    b2Mask = (fCenters >= opts.itpc_band2_hz(1)) & (fCenters <= opts.itpc_band2_hz(2));
end

IT_obs = cell(1,5);
IT_null = cell(1,5); % each has mean/std/nOK
NU = cell(1,5);
nEv = nan(1,5);

clim = [inf -inf];
yl_mean = [inf -inf];

for c = 1:5
    e = evList{c};
    e = e(:); e = e(isfinite(e));
    nEv(c) = numel(e);
    
    if numel(e) < opts.min_events
        IT_obs{c} = [];
        IT_null{c} = [];
        NU{c} = [];
        continue
    end
    
    [itpc_obs, nUsed] = itpc_compute_from_phasebank(phaseBank, e, t_rel, opts.min_events);
    IT_obs{c} = itpc_obs;
    NU{c} = nUsed;
    
    if opts.itpc_use_null
        Null = itpc_null_from_phasebank(phaseBank, numel(e), t_rel, opts.itpc_null_n);
        IT_null{c} = Null;
        
        itpc_corr = itpc_obs - Null.mean;
        clim(1) = min(clim(1), nanmin(itpc_corr(:)));
        clim(2) = max(clim(2), nanmax(itpc_corr(:)));
        
        m = nanmean(itpc_corr(fShowMask,:), 1);
        yl_mean(1) = min(yl_mean(1), nanmin(m));
        yl_mean(2) = max(yl_mean(2), nanmax(m));
    else
        clim(1) = min(clim(1), nanmin(itpc_obs(:)));
        clim(2) = max(clim(2), nanmax(itpc_obs(:)));
        
        m = nanmean(itpc_obs(fShowMask,:), 1);
        yl_mean(1) = min(yl_mean(1), nanmin(m));
        yl_mean(2) = max(yl_mean(2), nanmax(m));
    end
end

if ~all(isfinite(clim)) || clim(2) <= clim(1), clim = []; end
if ~all(isfinite(yl_mean)) || yl_mean(2) <= yl_mean(1), yl_mean = []; end

% ---- output struct
MetITPC = struct();
MetITPC.tag = tag;
MetITPC.t_rel = t_rel;
MetITPC.fCenters = fCenters;
MetITPC.opts_snapshot = opts;
MetITPC.group = struct([]);

% ---- plot: heatmap row + mean trace row (per group)
figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);

for c = 1:5
    MetITPC.group(c).name = colNames{c};
    MetITPC.group(c).nEv = nEv(c);
    MetITPC.group(c).itpc_obs = IT_obs{c};
    MetITPC.group(c).nUsed = NU{c};
    MetITPC.group(c).null = IT_null{c};
    
    % heatmap
    ax1 = subplot(2,5,c);
    if isempty(IT_obs{c})
        axis(ax1,'off');
        title(ax1, sprintf('%s (n=%d)', colNames{c}, nEv(c)), 'Interpreter','none');
        continue
    end
    
    if opts.itpc_use_null && ~isempty(IT_null{c})
        X = IT_obs{c} - IT_null{c}.mean;
        imagesc(ax1, t_rel, fCenters(fShowMask), X(fShowMask,:));
        title(ax1, sprintf('%s (n=%d) ITPC-NULL', colNames{c}, nEv(c)), 'Interpreter','none');
    else
        imagesc(ax1, t_rel, fCenters(fShowMask), IT_obs{c}(fShowMask,:));
        title(ax1, sprintf('%s (n=%d) ITPC', colNames{c}, nEv(c)), 'Interpreter','none');
    end
    axis(ax1,'xy'); box(ax1,'off');
    xline(ax1,0,'k-');
    ylabel(ax1,'Hz'); xlabel(ax1,'t (s)');
    if ~isempty(clim), caxis(ax1, clim); end
    
    % mean trace
    ax2 = subplot(2,5,5+c);
    if opts.itpc_use_null && ~isempty(IT_null{c})
        X = IT_obs{c} - IT_null{c}.mean;
        m = nanmean(X(fShowMask,:), 1);
        plot(ax2, t_rel, m, 'LineWidth', 1); hold(ax2,'on');
        
        % also show band means (optional)
        if any(b1Mask)
            d1 = nanmean(X(b1Mask,:),1);
            plot(ax2, t_rel, d1, 'LineWidth', 1);
        end
        if any(b2Mask)
            d2 = nanmean(X(b2Mask,:),1);
            plot(ax2, t_rel, d2, 'LineWidth', 1);
        end
        hold(ax2,'off');
        legend(ax2, ra_legend_names(any(b1Mask),any(b2Mask)), 'Box','off', 'Location','northwest');
        title(ax2,'mean (corrected)','Interpreter','none');
    else
        m = nanmean(IT_obs{c}(fShowMask,:), 1);
        plot(ax2, t_rel, m, 'LineWidth', 1);
        title(ax2,'mean (obs)','Interpreter','none');
    end
    xline(ax2,0,'k-'); yline(ax2,0,'k:');
    box(ax2,'off'); xlabel(ax2,'t (s)');
    if ~isempty(yl_mean), ylim(ax2, yl_mean); end
    
    % --- summary metrics for across-session usage
    postMask = (t_rel >= 0) & (t_rel <= 1.0);
    
    if opts.itpc_use_null && ~isempty(IT_null{c})
        Xsum = IT_obs{c} - IT_null{c}.mean;
    else
        Xsum = IT_obs{c};
    end
    
    mAll = nanmean(Xsum(fShowMask, postMask), 'all');
    
    MetITPC.group(c).mean_post_0_1s = mAll;
    
    if any(b1Mask)
        MetITPC.group(c).band1_post_0_1s = nanmean(Xsum(b1Mask, postMask), 'all');
        MetITPC.group(c).band1_peak = nanmax(nanmean(Xsum(b1Mask,:),1));
    end
    if any(b2Mask)
        MetITPC.group(c).band2_post_0_1s = nanmean(Xsum(b2Mask, postMask), 'all');
        MetITPC.group(c).band2_peak = nanmax(nanmean(Xsum(b2Mask,:),1));
    end
    
    MetITPC.group(c).nEvents = nEv(c);
end

sgtitle(sprintf('ITPC scan: %s', tag), 'Interpreter','none');

end

function L = ra_legend_names(hasB1, hasB2)
L = {'mean'};
if hasB1, L{end+1} = 'band1'; end %#ok<AGROW>
if hasB2, L{end+1} = 'band2'; end %#ok<AGROW>
end
