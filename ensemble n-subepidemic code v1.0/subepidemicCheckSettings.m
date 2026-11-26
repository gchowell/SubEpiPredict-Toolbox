function subepidemicCheckSettings(S)
%SUBEPIDEMICCHECKSETTINGS Check the settings of a fit before any fitting starts.
%   subepidemicCheckSettings(S) returns silently when the settings can run and
%   otherwise stops with error SubEpiPredict:InvalidSettings, which lists EVERY
%   problem found, each with the name of the setting in options.m.
%
%   S is a struct with the fields
%     method1, dist1, flag1, npatches_fixed, onset_fixed, smoothfactor1,
%     calibrationperiod1, numstartpoints, topmodelsx, B, DT, cumulative1,
%     outbreakx
%
%   Estimation method and error structure must match:
%     method1 = 0 (least squares)        dist1 = 0, 1 or 2
%     method1 = 1 (Poisson likelihood)   dist1 = 1
%     method1 = 3, 4, 5 (negative binomial likelihoods)   dist1 = method1
%     method1 = 6 (least absolute deviations)   dist1 = 6, or 1
%   options.m sets dist1 by itself for methods 1, 3, 4 and 5 but not for 0, so
%   changing only method1 to 0 left dist1 at 3: the bootstrap then had no noise
%   scale and the run stopped after the whole screen with an unrelated message.
%   method1 = 2 appears in older comments; no objective implements it.
%
%   The growth model flag1 is one integer from 0 to 5 for all sub-epidemics
%   (0 GGM, 1 generalized logistic, 2 generalized Richards, 3 logistic,
%   4 Richards, 5 Gompertz); a vector is not supported.
%
%   The comparison of the number of parameters with the length of the
%   calibration period stays in the fit driver.
%
%   See also Run_Fit_subepidemicFramework, options

if ~(isstruct(S) && isscalar(S))
    error('SubEpiPredict:InvalidSettings','The settings must be passed as one struct.');
end
names = {'method1','dist1','flag1','npatches_fixed','onset_fixed','smoothfactor1', ...
    'calibrationperiod1','numstartpoints','topmodelsx','B','DT','cumulative1','outbreakx'};
problems = {};
for k = 1:numel(names)
    if ~isfield(S,names{k})
        problems{end+1} = sprintf('%s is missing.',names{k}); %#ok<AGROW>
    end
end
if ~isempty(problems)
    error('SubEpiPredict:InvalidSettings','%s',localJoin(problems));
end

isInt = @(v) (isnumeric(v) || islogical(v)) && isreal(v) && isscalar(v) && isfinite(v) && v==fix(v);

% ---- estimation method and error structure ----
methodOK = isInt(S.method1) && any(S.method1==[0 1 3 4 5 6]);
if ~methodOK
    if isInt(S.method1) && S.method1==2
        problems{end+1} = 'method1 = 2 (Pearson chi-squared) is not implemented. Use 0, 1, 3, 4, 5 or 6.';
    else
        problems{end+1} = sprintf('method1 = %s is not valid. Use 0, 1, 3, 4, 5 or 6.',localShow(S.method1));
    end
end
distOK = isInt(S.dist1) && any(S.dist1==0:6);
if ~distOK
    problems{end+1} = sprintf('dist1 = %s is not valid. Use an integer from 0 to 6.',localShow(S.dist1));
end
if methodOK && distOK
    switch S.method1
        case 0
            allowed = [0 1 2];
            hint = ['Least squares (method1 = 0) needs dist1 = 0 (normal), 1 (Poisson) or 2 ', ...
                '(negative binomial with VAR = factor*mean). options.m does not set dist1 for method1 = 0.'];
        case 1
            allowed = 1;
            hint = 'The Poisson likelihood (method1 = 1) needs dist1 = 1.';
        case {3,4,5}
            allowed = S.method1;
            hint = sprintf('The negative binomial likelihood method1 = %d needs dist1 = %d.',S.method1,S.method1);
        case 6
            allowed = [6 1];
            hint = 'Least absolute deviations (method1 = 6) needs dist1 = 6 (Laplace) or 1 (Poisson).';
    end
    if ~any(S.dist1==allowed)
        problems{end+1} = sprintf('method1 = %d with dist1 = %d cannot run. %s',S.method1,S.dist1,hint);
    end
end

% ---- growth model ----
if ~(isInt(S.flag1) && S.flag1>=0 && S.flag1<=5)
    problems{end+1} = sprintf(['flag1 = %s is not valid. Use one integer from 0 to 5 for all ', ...
        'sub-epidemics; a vector is not supported.'],localShow(S.flag1));
end

% ---- sizes and switches ----
if ~(isInt(S.npatches_fixed) && S.npatches_fixed>=1)
    problems{end+1} = sprintf('npatches_fixed = %s must be a positive integer.',localShow(S.npatches_fixed));
end
if ~(isInt(S.onset_fixed) && any(S.onset_fixed==[0 1]))
    problems{end+1} = sprintf('onset_fixed = %s must be 0 or 1.',localShow(S.onset_fixed));
end
if ~(isInt(S.smoothfactor1) && S.smoothfactor1>=1)
    problems{end+1} = sprintf('smoothfactor1 = %s must be a positive integer (1 = no smoothing).',localShow(S.smoothfactor1));
end
if ~(isInt(S.calibrationperiod1) && S.calibrationperiod1>=2)
    problems{end+1} = sprintf('calibrationperiod1 = %s must be an integer of at least 2.',localShow(S.calibrationperiod1));
end
if ~(isInt(S.numstartpoints) && S.numstartpoints>=1)
    problems{end+1} = sprintf('numstartpoints = %s must be a positive integer.',localShow(S.numstartpoints));
end
if ~(isInt(S.topmodelsx) && S.topmodelsx>=1)
    problems{end+1} = sprintf('topmodelsx = %s must be a positive integer.',localShow(S.topmodelsx));
end
if ~(isInt(S.B) && S.B>=1)
    problems{end+1} = sprintf('B = %s (bootstrap realizations) must be a positive integer.',localShow(S.B));
end
if ~(isInt(S.DT) && any(S.DT==[1 7 365]))
    problems{end+1} = sprintf('DT = %s must be 1 (daily), 7 (weekly) or 365 (yearly).',localShow(S.DT));
end
if ~(isInt(S.cumulative1) && any(S.cumulative1==[0 1]))
    problems{end+1} = sprintf('cumulative1 = %s must be 0 (incidence file) or 1 (cumulative file).',localShow(S.cumulative1));
end
v = S.outbreakx;
if ~(isnumeric(v) && isreal(v) && isvector(v) && ~isempty(v) && all(isfinite(v(:))) && ...
        all(v(:)==fix(v(:))) && all(v(:)>=1))
    problems{end+1} = sprintf('outbreakx = %s must hold positive integers (columns of the data file).',localShow(v));
end

if ~isempty(problems)
    error('SubEpiPredict:InvalidSettings','%s',localJoin(problems));
end
end

% -----------------------------------------------------------------------------
function s = localShow(v)
if isnumeric(v) || islogical(v)
    if isempty(v)
        s = '[]';
    elseif numel(v)<=6
        s = mat2str(double(full(v)),6);
    else
        s = sprintf('[%d values]',numel(v));
    end
elseif ischar(v)
    s = ['''' v ''''];
else
    s = ['a ' class(v)];
end
end

function s = localJoin(problems)
if numel(problems)==1
    s = ['The settings cannot run. ' problems{1}];
else
    s = sprintf('The settings cannot run (%d problems):',numel(problems));
    for k = 1:numel(problems)
        s = sprintf('%s\n  %d. %s',s,k,problems{k});
    end
end
end
