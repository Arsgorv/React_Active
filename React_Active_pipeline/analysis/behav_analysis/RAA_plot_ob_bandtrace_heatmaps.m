function RAA_plot_ob_bandtrace_heatmaps(RES, opts)
if nargin < 2, opts = struct(); end
if ~isfield(opts,'caxis_prct'), opts.caxis_prct = [5 95]; end
if ~isfield(opts,'groups'), opts.groups = {'all','tar','ref','nosound','nomotor'}; end

sessLabels = RES.sessions;

for b = 1:numel(RES.band)
    bname = RES.band(b).name;
    
    % collect matrices for requested groups and compute shared clim
    Ms = cell(numel(opts.groups),1);
    t_rel0 = [];
    for gi = 1:numel(opts.groups)
        gname = string(opts.groups{gi});
        gIdx = find(strcmpi(string({RES.band(b).group.name}), gname),1,'first');
        if isempty(gIdx), Ms{gi} = []; continue, end
        Ms{gi} = RES.band(b).group(gIdx).M;
        if isempty(t_rel0), t_rel0 = RES.band(b).group(gIdx).t_rel; end
    end
    
    v = [];
    for gi = 1:numel(Ms)
        if isempty(Ms{gi}), continue, end
        vv = Ms{gi}(isfinite(Ms{gi}));
        v = [v; vv(:)]; %#ok<AGROW>
    end
    clim = [];
    if ~isempty(v)
        lo = prctile(v, opts.caxis_prct(1));
        hi = prctile(v, opts.caxis_prct(2));
        if isfinite(lo) && isfinite(hi) && hi > lo, clim = [lo hi]; end
    end
    
    figure('Color','w','Units','pixels','Position',[50 50 1700 700], ...
        'Name', sprintf('OB_%s_bandtrace_heat', char(bname)));
    
    for gi = 1:numel(opts.groups)
        subplot(1,numel(opts.groups),gi);
        M = Ms{gi};
        if isempty(M) || all(~isfinite(M(:)))
            axis off; title(opts.groups{gi},'Interpreter','none'); continue
        end
        
        imagesc(t_rel0, 1:size(M,1), M);
        axis xy
        xlabel('t from stimOn (s)');
        if gi==1
            ylabel('session');
            yticks(1:numel(sessLabels));
            yticklabels(cellstr(sessLabels));
        else
            yticks([]);
        end
        title(sprintf('%s | %s', char(bname), opts.groups{gi}), 'Interpreter','none');
        try colormap(gca,'viridis'); catch, colormap(gca,parula); end
        if ~isempty(clim), caxis(clim); end
        colorbar; box off
        xline(0,'k-');
    end
    
    % ---- after plotting existing group heatmaps for this band ----
    % Find tar/ref groups
    gNames = strings(numel(RES.band(b).group),1);
    for gg = 1:numel(RES.band(b).group)
        gNames(gg) = string(RES.band(b).group(gg).name);
    end
    gTar = find(strcmpi(gNames,'tar'), 1);
    gRef = find(strcmpi(gNames,'ref'), 1);
    
    if ~isempty(gTar) && ~isempty(gRef)
        t_rel = RES.band(b).group(gTar).t_rel(:);
        
        Mtar = RES.band(b).group(gTar).M;
        Mref = RES.band(b).group(gRef).M;
        Mdiff = Mtar - Mref; % Tar - Ref
        
        if any(isfinite(Mdiff(:)))
            figure('Color','w','Units','pixels','Position',[50 50 900 750], ...
                'Name', sprintf('OB_%s_TarMinusRef_bandtrace_heat', char(bname)));
            
            imagesc(t_rel, 1:size(Mdiff,1), Mdiff);
            axis xy
            xlabel(sprintf('t from %s (s)', char(RES.alignName)));
            ylabel('session');
            title(sprintf('%s | Tar-Ref | band-trace (BL-corr log10)', char(bname)), 'Interpreter','none');
            yticks(1:numel(sessLabels));
            yticklabels(cellstr(sessLabels));
            colorbar
            box off
            xline(0,'k-');
            
            colormap(gca, ra_redblue(256));
            
            % symmetric caxis (robust)
            v = Mdiff(isfinite(Mdiff));
            if ~isempty(v)
                lim = prctile(abs(v), 95);
                if isfinite(lim) && lim > 0
                    caxis([-lim lim]);
                end
            end
        end
    end
    
end



end

% ---- helper (place at end of file) ----
function cmap = ra_redblue(n)
% Simple diverging red-white-blue colormap (R2018b safe)
if nargin < 1, n = 256; end
n2 = floor(n/2);
t  = linspace(0,1,n2)';

% blue -> white
c1 = [t t ones(n2,1)];           % (0,0,1) to (1,1,1)
% white -> red
t2 = linspace(1,0,n-n2)';        % start at white
c2 = [ones(n-n2,1) t2 t2];       % (1,1,1) to (1,0,0)

cmap = [c1; c2];
end