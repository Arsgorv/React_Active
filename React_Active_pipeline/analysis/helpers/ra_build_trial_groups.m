function G = ra_build_trial_groups(Baphy)

n = Baphy.n_trials;

G = struct();
G.names = {'all','tar','ref','nosound','nomotor'};
G.mask  = false(n, numel(G.names));

% all
G.mask(:,1) = true;

% tar/ref from Baphy.trial.type
isTar = false(n,1);
isRef = false(n,1);

if isfield(Baphy,'trial') && isfield(Baphy.trial,'type') && ~isempty(Baphy.trial.type)
    tt = Baphy.trial.type;
    if iscell(tt)
        tt = lower(string(tt(:)));
    else
        tt = lower(string(tt));
    end
    isTar = contains(tt,'target');
    isRef = contains(tt,'reference');
end

G.mask(:,2) = isTar;
G.mask(:,3) = isRef;

% nosound / nomotor directly present
isNoSound = false(n,1);
isNoMotor = false(n,1);

if isfield(Baphy,'trial')
    if isfield(Baphy.trial,'is_nosound')
        isNoSound = logical(Baphy.trial.is_nosound(:));
    end
    if isfield(Baphy.trial,'is_nomotor')
        isNoMotor = logical(Baphy.trial.is_nomotor(:));
    end
end

G.mask(:,4) = isNoSound;
G.mask(:,5) = isNoMotor;

end
