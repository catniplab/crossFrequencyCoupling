function estimators = CFCestimatorFactory(names)

if nargin < 1 || isempty(names)
    estimators = {'MCS', 'PLV', 'ESC', 'R2R', 'EAA'};
    return
end

if ischar(names)
    if strcmpi(names, 'all')
        names = {'MCS', 'PLV', 'ESC', 'R2R', 'EAA'};
    else
        estimators = resolveOne(names);
        return
    end
end

for k = 1:numel(names)
    estimators(k) = resolveOne(names{k});
end

end % CFCestimatorFactory

function estimator = resolveOne(name)

switch(upper(name))
    case {'MCS', 'CANOLTY2006'}
        estimator.ID = 'MCS';
        estimator.ref = 'Canolty2006';
        estimator.desc = 'magnitude of mean compound signal';
        estimator.handle = @estimateCFC_MCS_Canolty2006;
    case {'PLV', 'LACHAUX1999'}
        estimator.ID = 'PLV';
        estimator.ref = 'Lachaux1999';
        estimator.desc = 'phase locking value';
        estimator.handle = @estimateCFC_PLV_Lachaux1999;
    case {'ESC', 'BRUNS2004'}
        estimator.ID = 'ESC';
        estimator.ref = 'Bruns2004';
        estimator.desc = 'envelope-signal-correlation';
        estimator.handle = @estimateCFC_ESC_Bruns2004;
    case {'R2R', 'Penny2008'}
        estimator.ID = 'R2R';
        estimator.ref = 'Penny2008';
        estimator.desc = 'r^2 of phase-amplitude-regression';
        estimator.handle = @estimateCFC_R2R_Penny2008;
    case {'EAA', 'Tort2008'}
        estimator.ID = 'EAA';
        estimator.ref = 'Tort2008';
        estimator.desc = 'Entropy of average amplitude per phase';
        estimator.handle = @estimateCFC_EAA_Tort2008;
    otherwise
        error('No such estimator [%s]', name);
end

end % resolveOne
