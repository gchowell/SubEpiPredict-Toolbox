function varargout = subepidemicFitSignature(saved)
%SUBEPIDEMICFITSIGNATURE What made a fit: the fit revision and the settings that shape it.
%   signature = subepidemicFitSignature()
%   [same,differences] = subepidemicFitSignature(saved)
%
%   Without input: the signature of a fit made now. Run_Fit_subepidemicFramework.m
%   stores it as the variable fitSignature in every ranked file.
%     signature.fitRevision  changes whenever a revision of the toolbox can
%                            change the numbers of a fit (see FITREVISION below)
%     signature.settings     the settings of subepidemicFitSettings.m that
%                            shape the fit: scales, bounds, likelihood
%                            option, optimizer options, attempts, parallel
%                            use, polish, start points, random seed
%   Settings that only change reports are left out (PIIscale,
%   evidenceRatioLimit, notIdentifiedRatio, baselineWindow,
%   maxSkippedRefitFraction): changing one of them is no reason to fit again.
%
%   With the signature SAVED of a stored fit: SAME is true if a fit made now
%   would have the same signature; DIFFERENCES names what differs
%   ('fit revision' or the names of settings). A fit stored before revision
%   m1 has no signature; pass [] for it: SAME is false.
%
%   Used by subepidemicRollingForecast.m, which fits a date again unless its
%   stored fit has the current signature (and the same B, numstartpoints and
%   topmodelsx, which it compares itself). The method, distribution, kernel,
%   smoothing and calibration period are part of the file name of a fit and
%   need no check.
%
%   Added in revision m1 (2026-11-18).

% FITREVISION: to be changed in every revision that can change a fit's numbers.
%   'm1'  first signature; the fit code is that of revision m (extra start
%         points at the first threshold, ordered labels).
fitRevision = 'm1';

reportOnly = {'PIIscale','evidenceRatioLimit','notIdentifiedRatio','baselineWindow','maxSkippedRefitFraction'};
settings = subepidemicFitSettings();
settings = rmfield(settings, reportOnly(isfield(settings,reportOnly)));
signature = struct('fitRevision',fitRevision,'settings',settings);

if nargin==0
    varargout{1} = signature;
    return
end

differences = {};
if ~(isstruct(saved) && isscalar(saved) && isfield(saved,'fitRevision') && isfield(saved,'settings') && isstruct(saved.settings))
    differences = {'no signature (fit made before revision m1)'};
else
    if ~isequal(saved.fitRevision,signature.fitRevision)
        differences{end+1} = sprintf('fit revision (%s then, %s now)',char(saved.fitRevision),signature.fitRevision);
    end
    names = union(fieldnames(saved.settings),fieldnames(signature.settings));
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(saved.settings,name) || ~isfield(signature.settings,name) || ...
                ~isequaln(saved.settings.(name),signature.settings.(name))
            differences{end+1} = name; %#ok<AGROW>
        end
    end
end
varargout{1} = isempty(differences);
varargout{2} = differences;
end
