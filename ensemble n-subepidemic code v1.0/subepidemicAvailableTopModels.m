function topmodels = subepidemicAvailableTopModels(topmodels,npatches_fixed,onset_fixed,smoothfactor1,cadfilename2,flag1,method1,dist1,calibrationperiod1)
%SUBEPIDEMICAVAILABLETOPMODELS Limit the ranked models to those the saved fit holds.
%   topmodels = subepidemicAvailableTopModels(topmodels, ...) takes the ranks a
%   script is about to process (1:topmodelsx from options.m) and returns the
%   leading part for which a ranked file of the fit exists in ./output. The
%   other arguments are the pieces of the file name, as the scripts use them.
%
%   options.m can ask for more top models than the fit saved: the screen had
%   fewer candidates, or topmodelsx was raised after the fit. The scripts then
%   stopped with a missing-file or index error. Now they process the ranks
%   that exist and warn (SubEpiPredict:TopModelsCapped).
%
%   If no ranked file exists at all, the list is returned unchanged and the
%   script stops at its own load, as before.

prefix = strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed), ...
    '-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2, ...
    '-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1), ...
    '-calibrationperiod-',num2str(calibrationperiod1),'-rank-');
nAvailable = 0;
for k = 1:numel(topmodels)
    if exist(strcat(prefix,num2str(topmodels(k)),'.mat'),'file')==2
        nAvailable = nAvailable+1;
    else
        break
    end
end
if nAvailable>0 && nAvailable<numel(topmodels)
    subepidemicWarning('SubEpiPredict:TopModelsCapped', ...
        ['options.m asks for %d top-ranked models, but the saved fit has %d. ', ...
         'Ranks 1 to %d are used.'],numel(topmodels),nAvailable,nAvailable);
    topmodels = topmodels(1:nAvailable);
end
end
