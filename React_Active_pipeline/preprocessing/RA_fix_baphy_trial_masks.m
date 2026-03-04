function B = RA_fix_baphy_trial_masks(B)
% RA_fix_baphy_trial_masks
% Ensures per-trial logical masks match B.n_trials after concatenation/merge.

if ~isfield(B,'trial') || ~isstruct(B.trial)
    return
end

% ground truth n
n = [];
if isfield(B,'n_trials') && ~isempty(B.n_trials)
    n = double(B.n_trials);
end
if isempty(n) || ~isfinite(n) || n<=0
    if isfield(B.trial,'type'), n = numel(B.trial.type);
    else, return
    end
end
B.n_trials = n;
B.trial_id = (1:n)';

% ---- is_nosound ----
need = (~isfield(B.trial,'is_nosound')) || (numel(B.trial.is_nosound) ~= n);
if need
    m = false(n,1);
    if isfield(B,'idx') && isfield(B.idx,'NoSound') && ~isempty(B.idx.NoSound)
        ix = B.idx.NoSound(:);
        ix = ix(isfinite(ix) & ix>=1 & ix<=n);
        m(ix) = true;
    end
    B.trial.is_nosound = m;
end

% ---- is_nomotor ----
need = (~isfield(B.trial,'is_nomotor')) || (numel(B.trial.is_nomotor) ~= n);
if need
    m = false(n,1);
    if isfield(B,'idx')
        if isfield(B.idx,'NoMotorExclusive') && ~isempty(B.idx.NoMotorExclusive)
            ix = B.idx.NoMotorExclusive(:);
            ix = ix(isfinite(ix) & ix>=1 & ix<=n);
            m(ix) = true;
        elseif isfield(B.idx,'NoMotor') && ~isempty(B.idx.NoMotor)
            ix = B.idx.NoMotor(:);
            ix = ix(isfinite(ix) & ix>=1 & ix<=n);
            m(ix) = true;
        end
    end
    B.trial.is_nomotor = m;
end

end