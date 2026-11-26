function subepidemicSkippedRefits(label,rank1,nSkipped,nRefits,firstMessage)
%SUBEPIDEMICSKIPPEDREFITS One policy for refits that could not be simulated.
%   subepidemicSkippedRefits(label,rank1,nSkipped,nRefits,firstMessage) is
%   called by the forecast and Rt scripts after the loop over the refits of
%   one ranked model.
%
%   nSkipped of the nRefits could not be simulated over the forecast period
%   and were left out (subepidemicSimulateRefit.m). Then:
%     none skipped                         nothing happens
%     all skipped                          error SubEpiPredict:NoForecastCurves
%     share above maxSkippedRefitFraction  error SubEpiPredict:TooManyRefitsSkipped
%     otherwise                            warning SubEpiPredict:ForecastCurvesSkipped
%
%   maxSkippedRefitFraction is a setting of options_fit_settings.m, default
%   0.05. With 0 the script stops at the first refit that fails (the behaviour
%   before revision g); with 1 it stops only when no refit is left. The refits
%   that fail tend to be the extreme parameter sets, so leaving many out would
%   narrow the interval; that is why there is a limit.

if nSkipped<=0
    return
end
settings = subepidemicFitSettings();
limit = settings.maxSkippedRefitFraction;
nUsed = nRefits-nSkipped;
if nUsed<=0
    error('SubEpiPredict:NoForecastCurves', ...
        '%s, ranked model %d: none of the %d refits could be simulated over the forecast period. First failure: %s', ...
        label,rank1,nRefits,firstMessage);
end
if nSkipped/nRefits>limit
    error('SubEpiPredict:TooManyRefitsSkipped', ...
        ['%s, ranked model %d: %d of %d refits could not be simulated over the forecast period. ', ...
         'That is more than the share allowed (maxSkippedRefitFraction = %g). First failure: %s'], ...
        label,rank1,nSkipped,nRefits,limit,firstMessage);
end
subepidemicWarning('SubEpiPredict:ForecastCurvesSkipped', ...
    ['%s, ranked model %d: %d of %d refits could not be simulated over the forecast period and were left out. ', ...
     'The result rests on %d refits. First failure: %s'], ...
    label,rank1,nSkipped,nRefits,nUsed,firstMessage);
end
