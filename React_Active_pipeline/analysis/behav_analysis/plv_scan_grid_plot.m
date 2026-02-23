function MetPLV = plv_scan_grid_plot(LFP, ev, G, opts, tag, phaseBank)

if ~isfield(opts,'min_events'),    opts.min_events = 10; end
if ~isfield(opts,'plv_f_centers'), opts.plv_f_centers = 2:2:120; end
if ~isfield(opts,'plv_bw_hz'),     opts.plv_bw_hz = 2; end

% null correction options
if ~isfield(opts,'plv_use_null'),  opts.plv_use_null = true; end
if ~isfield(opts,'plv_null_n'),    opts.plv_null_n = 500; end

fCenters = opts.plv_f_centers(:)'; bw = opts.plv_bw_hz;

evList   = {ev.all, ev.tar, ev.ref, ev.nosound, ev.nomotor};
colNames = {'all','tar','ref','nosound','nomotor'};

% phase bank from prepare_phase_bank (reuse ITPC bank interface)
if nargin < 6 || isempty(phaseBank)
    tmp = opts;
    tmp.itpc_f_centers = fCenters;
    tmp.itpc_bw_hz = bw;
    phaseBank = prepare_phase_bank(LFP, tmp);
end
if ~isequal(phaseBank.fCenters, fCenters) || phaseBank.bw ~= bw
    tmp = opts; tmp.itpc_f_centers = fCenters; tmp.itpc_bw_hz = bw;
    phaseBank = prepare_phase_bank(LFP, tmp);
end

PLV  = cell(1,5);
PRAY = cell(1,5);
NUSE = cell(1,5);
NULL = cell(1,5);

yl_plv = [inf -inf];
yl_z   = [inf -inf];
yl_n   = [inf -inf];

for c = 1:5
    e = evList{c}; e = e(:); e = e(isfinite(e));
    if numel(e) < opts.min_events
        PLV{c} = []; PRAY{c} = []; NUSE{c} = []; NULL{c} = [];
        continue
    end

    [plv, pRay, nUse] = plv_compute_from_phasebank(phaseBank, e, opts.min_events);
    PLV{c} = plv;
    PRAY{c} = pRay;
    NUSE{c} = nUse;

    if opts.plv_use_null
        Null = plv_null_from_phasebank(phaseBank, numel(e), opts.plv_null_n);
        NULL{c} = Null;
        z = (plv - Null.mean) ./ max(Null.std, 1e-9);
        yl_z(1) = min(yl_z(1), nanmin(z));
        yl_z(2) = max(yl_z(2), nanmax(z));
    end

    yl_plv(1) = min(yl_plv(1), nanmin(plv));
    yl_plv(2) = max(yl_plv(2), nanmax(plv));
    yl_n(1) = min(yl_n(1), nanmin(nUse));
    yl_n(2) = max(yl_n(2), nanmax(nUse));
end

if ~all(isfinite(yl_plv)) || yl_plv(2)<=yl_plv(1), yl_plv=[]; end
if ~all(isfinite(yl_z))   || yl_z(2)<=yl_z(1),     yl_z=[]; end
if ~all(isfinite(yl_n))   || yl_n(2)<=yl_n(1),     yl_n=[]; end

MetPLV = struct();
MetPLV.tag = tag;
MetPLV.fCenters = fCenters;
MetPLV.opts_snapshot = opts;
MetPLV.group = struct([]);

figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);

for c = 1:5
    MetPLV.group(c).name = colNames{c};
    MetPLV.group(c).plv = PLV{c};
    MetPLV.group(c).pRay = PRAY{c};
    MetPLV.group(c).nUse = NUSE{c};
    MetPLV.group(c).null = NULL{c};

    % row 1: PLV obs
    ax1 = subplot(3,5,c);
    if isempty(PLV{c})
        axis(ax1,'off');
        title(ax1, sprintf('%s', colNames{c}), 'Interpreter','none');
    else
        plot(ax1, fCenters, PLV{c}, 'LineWidth', 1);
        xlabel(ax1,'Hz'); ylabel(ax1,'PLV');
        title(ax1, sprintf('%s', colNames{c}), 'Interpreter','none');
        box(ax1,'off');
        if ~isempty(yl_plv), ylim(ax1, yl_plv); end
    end

    % row 2: zscore vs null (preferred) OR -log10(pRay)
    ax2 = subplot(3,5,5+c);
    if isempty(PLV{c})
        axis(ax2,'off');
    else
        if opts.plv_use_null && ~isempty(NULL{c})
            z = (PLV{c} - NULL{c}.mean) ./ max(NULL{c}.std, 1e-9);
            plot(ax2, fCenters, z, 'LineWidth', 1);
            xlabel(ax2,'Hz'); ylabel(ax2,'z(PLV-null)');
            yline(ax2,0,'k:');
            box(ax2,'off');
            if ~isempty(yl_z), ylim(ax2, yl_z); end
        else
            plot(ax2, fCenters, -log10(PRAY{c}), 'LineWidth', 1);
            xlabel(ax2,'Hz'); ylabel(ax2,'-log10(pRay)');
            box(ax2,'off');
        end
    end

    % row 3: n used
    ax3 = subplot(3,5,10+c);
    if isempty(PLV{c})
        axis(ax3,'off');
    else
        plot(ax3, fCenters, NUSE{c}, 'LineWidth', 1);
        xlabel(ax3,'Hz'); ylabel(ax3,'n');
        box(ax3,'off');
        if ~isempty(yl_n), ylim(ax3, yl_n); end
    end
end

sgtitle(sprintf('PLV scan: %s', tag), 'Interpreter','none');

end
