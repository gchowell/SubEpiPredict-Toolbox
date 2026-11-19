function info = subepidemicSeed(label)
%SUBEPIDEMICSEED Apply the random seed setting and record the generator state.
%   info = subepidemicSeed(label) is called at the start of the fit driver
%   (once per area) and of the forecast, fit-plot and Rt scripts.
%
%   Setting rngSeed (options_fit_settings.m; default 1 since revision l):
%     an integer the generator is reset with rng(rngSeed), so the same call
%                on the same files gives the same numbers again.
%     empty      the generator is left as it is (the default until revision
%                k). Results differ from run to run.
%
%   info.seed   the seed used, or [] when none was set
%   info.state  the generator state at this point (the output of rng). The
%               fit driver saves info with its results, so a fit made
%               without a seed can still be repeated: load a ranked file
%               and call rng(rngInfo.state) before running the driver again.
%
%   The start points of the search, the bootstrap data and the forecast noise
%   are drawn on the MATLAB client, so one seed fixes them with and without a
%   parallel pool. LABEL is printed with the seed.

if nargin<1 || isempty(label)
    label = 'SubEpiPredict';
end
settings = subepidemicFitSettings();
info = struct('seed',[],'state',[]);
info.seed = settings.rngSeed;
if ~isempty(info.seed)
    rng(info.seed);
    fprintf('%s: random seed %d (setting rngSeed).\n',label,info.seed);
end
info.state = rng;
end
