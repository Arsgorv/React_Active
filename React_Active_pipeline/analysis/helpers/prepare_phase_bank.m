function phaseBank = prepare_phase_bank(LFP, opts)
% prepare_phase_bank
% Precompute analytic phase time series for all requested frequency bands ONCE.
% Reuse in ITPC/PLV to avoid repeated FilterLFP+hilbert.
%
% phaseBank fields:
%   .t_start
%   .fs
%   .nSamp
%   .fCenters
%   .bw
%   .phCell   {1 x nF} each [nSamp x 1] phase (double)
%
% Notes:
% - Uses FilterLFP(LFP,[f1 f2],fs) as you did.
% - Keeps your inferred fs to match Range(LFP,'s') indexing.

if ~isfield(opts,'min_events'), opts.min_events = 10; end 
if ~isfield(opts,'itpc_f_centers'), opts.itpc_f_centers = 2:2:120; end
if ~isfield(opts,'itpc_bw_hz'), opts.itpc_bw_hz = 2; end

t0 = double(Range(LFP,'s'));
if numel(t0) < 2
    error('prepare_phase_bank: LFP too short');
end

t_start = t0(1);
fs = round(1/median(diff(t0)));
if ~isfinite(fs) || fs <= 0
    fs = 1024;
end

fCenters = opts.itpc_f_centers(:)';
bw = opts.itpc_bw_hz;

% compute once
nF = numel(fCenters);
phCell = cell(1,nF);

for k = 1:nF
    f1 = max(0.5, fCenters(k) - bw/2);
    f2 = fCenters(k) + bw/2;

    Fil = FilterLFP(LFP, [f1 f2], fs);
    x = double(Data(Fil));
    x = x(:);

    % analytic phase
    phCell{k} = angle(hilbert(x));
end

phaseBank = struct();
phaseBank.t_start = t_start;
phaseBank.fs = fs;
phaseBank.nSamp = numel(phCell{1});
phaseBank.fCenters = fCenters;
phaseBank.bw = bw;
phaseBank.phCell = phCell;

end
