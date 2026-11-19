function R = subepidemicRollingForecast(forecastDates,forecastingperiod,outbreakx,weight_type1,varargin)
%SUBEPIDEMICROLLINGFORECAST Fit and forecast from a list of dates and score every horizon.
%   R = subepidemicRollingForecast(forecastDates)
%   R = subepidemicRollingForecast(forecastDates,forecastingperiod,outbreakx,weight_type1)
%   R = subepidemicRollingForecast(...,'refit',true,'closeFigures',false)
%
%   One forecast from one date says little about a model. This function
%   repeats the fit (Run_Fit_subepidemicFramework.m) and the scored forecast
%   (plotForecast_subepidemicFramework.m) from every date of a list, with the
%   settings of options.m and options_forecast.m, and averages the scores of
%   each horizon over the dates, next to the flat baseline
%   (subepidemicBaselineForecast.m).
%
%   forecastDates      cell array of dates 'mm-dd-yyyy' (the last date of the
%                      data used for each fit), or a vector of datenum values.
%                      Example, every 2 weeks:
%                      cellstr(datestr(datenum(2020,4,20):14:datenum(2020,8,24),'mm-dd-yyyy'))
%   forecastingperiod  horizon in time units ([] = the value in options_forecast.m)
%   outbreakx          one area, a column of the data file ([] = options.m)
%   weight_type1       ensemble weights ([] = options_forecast.m)
%   'refit'            false (default): a date whose fit is already in ./output,
%                      made with the same B, numstartpoints and topmodelsx,
%                      by the same fit revision with the same fit settings
%                      (subepidemicFitSignature.m) and complete, is not
%                      fitted again, so an interrupted run can
%                      be continued. true: every date is fitted again.
%   'closeFigures'     true (default): the figures a date opens are closed
%                      before the next date.
%   'quiet'            false (default): each date prints what a fit and a
%                      forecast print. true: that text goes to a log file per
%                      date in ./output (rolling-log-...txt) and each date
%                      prints a few lines: data, fit, warnings, the score of
%                      ranked model 1 and the time. (Revision m.)
%   'keepReplaced'     true (default): before a date is fitted again, the fit
%                      files of that date already in ./output are copied to
%                      ./output/replaced-by-rolling/, unless a copy is there
%                      already. So the fit you made of a date with other
%                      settings is not lost. (Revision m.)
%
%   DATA. The data of every date come from one file: the file of the most
%   recent date (datevecend1 in options.m), which also supplies the later
%   observations the forecasts are scored against. For a date without its own
%   file in ./input, the file is written from the first rows of that most
%   recent file. Such a file holds the data as they are known today, not as
%   they were known on that date: later corrections of the counts are
%   included, which favours every forecast a little. A file that already
%   exists is used as it is, and a difference from the most recent file is
%   reported.
%
%   SCORES. As in the file performance-forecasting-byHorizon-...csv, the score
%   "over h" is the average over the first h time units of a forecast; the
%   score "at step h" is that of time unit h alone. Both are averaged over
%   the dates. Three comparisons with the baseline are given for the WIS:
%     ratio of means   mean WIS of the model / mean WIS of the baseline; dates
%                      with large counts weigh more
%     typical ratio    geometric mean over the dates of (WIS of the model /
%                      WIS of the baseline); every date weighs the same
%     share better     the share of dates on which the model's WIS is below
%                      the baseline's
%   Values below 1 (share above 0.5) favour the model. Forecasts from nearby
%   dates overlap and are not independent, so a difference seen on a handful
%   of dates is weak evidence.
%
%   OUTPUT. R.dates, R.status (one row per date), R.scores (one row per date,
%   model and horizon), R.summary (one row per model and horizon), the column
%   names in R.scoresColumns and R.summaryColumns, and the names of the two
%   files written to ./output (rolling-forecast-scores-... and
%   rolling-forecast-summary-...). Model codes: 0 = baseline, k = ranked k,
%   100+k = Ensemble(k).
%
%   Added in revision j (2026-11-17).

% ---------------------------------------------------------------- inputs
if nargin<1 || isempty(forecastDates)
    error('SubEpiPredict:RollingInput','Give the forecast dates: a cell array of ''mm-dd-yyyy'' or a vector of datenum values.');
end
if nargin<2, forecastingperiod=[]; end
if nargin<3, outbreakx=[]; end
if nargin<4, weight_type1=[]; end
refit=false;
closeFigures=true;
quiet=false;
keepReplaced=true;
if mod(numel(varargin),2)~=0
    error('SubEpiPredict:RollingInput','Options are given as name and value pairs, e.g. ''refit'',true.');
end
for k=1:2:numel(varargin)
    switch lower(char(varargin{k}))
        case 'refit'
            refit=logical(varargin{k+1});
        case 'closefigures'
            closeFigures=logical(varargin{k+1});
        case 'quiet'
            quiet=logical(varargin{k+1});
        case 'keepreplaced'
            keepReplaced=logical(varargin{k+1});
        otherwise
            error('SubEpiPredict:RollingInput','Unknown option ''%s''. The options are ''refit'', ''closeFigures'', ''quiet'' and ''keepReplaced''.',char(varargin{k}));
    end
end

global method1 npatches_fixed onset_fixed dist1 smoothfactor1 calibrationperiod1 %#ok<GVMIS>

[cumulative1,outbreakx_INP,~,cadregion,caddisease,datatype,DT,datevecfirst1,datevecend1,numstartpoints,topmodelsSetting,B,flag1]=options;
[getperformance,~,forecastingperiod_INP,weight_type1_INP]=options_forecast;

% The fit driver clears these global variables when it starts, so their values
% are kept here in ordinary variables.
methodFit=method1; distFit=dist1; nSubepidemics=npatches_fixed; onsetFixed=onset_fixed;
smoothing=smoothfactor1; calibrationLength=calibrationperiod1;
if nSubepidemics==1
    topmodelsSetting=1; % as in the fit driver
end

if isempty(forecastingperiod), forecastingperiod=forecastingperiod_INP; end
if isempty(outbreakx), outbreakx=outbreakx_INP; end
if isempty(weight_type1), weight_type1=weight_type1_INP; end
if ~isscalar(outbreakx)
    error('SubEpiPredict:RollingInput','Give one area (outbreakx is a single column of the data file); run the function once per area.');
end
if ~(isscalar(forecastingperiod) && forecastingperiod>=1 && forecastingperiod==round(forecastingperiod))
    error('SubEpiPredict:RollingInput','The forecast horizon must be a whole number of at least 1 (it is %s).',mat2str(forecastingperiod));
end
if getperformance~=1
    error('SubEpiPredict:RollingInput','options_forecast.m has getperformance = %d. Set it to 1: without it no forecast is scored.',getperformance);
end

switch DT
    case 1
        cadtemporal='daily';
    case 7
        cadtemporal='weekly';
    case 365
        cadtemporal='yearly';
    otherwise
        error('SubEpiPredict:RollingInput','DT = %s in options.m; the toolbox reads daily (1), weekly (7) and yearly (365) data.',mat2str(DT));
end
if cumulative1==1
    filePrefix=strcat('cumulative-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-');
else
    filePrefix=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-');
end

% dates as numbers, in increasing order, without repeats
if iscell(forecastDates) || ischar(forecastDates)
    forecastDates=cellstr(forecastDates);
    dateNumbers=zeros(numel(forecastDates),1);
    for k=1:numel(forecastDates)
        d=strtrim(forecastDates{k});
        if numel(d)~=10 || d(3)~='-' || d(6)~='-' || any(~ismember(d([1 2 4 5 7:10]),'0123456789'))
            error('SubEpiPredict:RollingInput','Forecast date %d is ''%s''. Dates are written mm-dd-yyyy, e.g. 05-11-2020.',k,d);
        end
        monthDayYear=[str2double(d(1:2)) str2double(d(4:5)) str2double(d(7:10))];
        dateNumbers(k)=datenum([monthDayYear(3) monthDayYear(1) monthDayYear(2)]);
        backAgain=datevec(dateNumbers(k));
        if ~isequal(backAgain(1:3),monthDayYear([3 1 2]))
            error('SubEpiPredict:RollingInput','Forecast date %d, ''%s'', is not a calendar date (dates are written mm-dd-yyyy).',k,d);
        end
    end
else
    dateNumbers=double(forecastDates(:));
end
dateNumbers=unique(dateNumbers);
nDates=numel(dateNumbers);

% ---------------------------------------------- the most recent data file
latestName=strcat(filePrefix,datestr(datenum(datevecend1),'mm-dd-yyyy'),'.txt');
if exist(strcat('./input/',latestName),'file')~=2
    error('SubEpiPredict:RollingInput', ...
        ['The file ./input/%s was not found. It is named by datevecend1 in options.m and holds the whole series: ', ...
         'the data of every forecast date and the later observations the forecasts are scored against.'],latestName);
end
latest=load(strcat('./input/',latestName));
if outbreakx<1 || outbreakx>size(latest,2) || outbreakx~=round(outbreakx)
    error('SubEpiPredict:RollingInput','outbreakx = %s, but the data file %s has %d column(s).',mat2str(outbreakx),latestName,size(latest,2));
end
nLatest=size(latest,1);

% number of data rows up to each date; later data needed: rows+1 ... rows+horizon
if DT==365
    dv=datevec(dateNumbers);
    rowsAtDate=dv(:,1)-datevecfirst1(1)+1;
    onGrid=true(nDates,1);
else
    steps=(dateNumbers-datenum(datevecfirst1))/DT;
    onGrid=abs(steps-round(steps))<1e-9;
    rowsAtDate=round(steps)+1;
end
if any(~onGrid)
    k=find(~onGrid,1);
    error('SubEpiPredict:RollingInput', ...
        ['Forecast date %s is not one of the dates of the series: the series starts on %s and has one value every %d days. ', ...
         'The nearest dates of the series are %s and %s.'], ...
        datestr(dateNumbers(k),'mm-dd-yyyy'),datestr(datenum(datevecfirst1),'mm-dd-yyyy'),DT, ...
        datestr(datenum(datevecfirst1)+floor(steps(k))*DT,'mm-dd-yyyy'),datestr(datenum(datevecfirst1)+ceil(steps(k))*DT,'mm-dd-yyyy'));
end
if any(rowsAtDate<1)
    k=find(rowsAtDate<1,1);
    error('SubEpiPredict:RollingInput','Forecast date %s is before the first date of the series (%s).', ...
        datestr(dateNumbers(k),'mm-dd-yyyy'),datestr(datenum(datevecfirst1),'mm-dd-yyyy'));
end

% status columns: 1 date number, 2 data rows, 3 input file (1 written, 2 existing and equal to the
% most recent file, 3 existing and different), 4 fit (1 new, 2 reused), 5 scored (1/0), 6 warnings, 7 seconds
status=[dateNumbers rowsAtDate nan(nDates,2) zeros(nDates,2) nan(nDates,1)];
notes=cell(nDates,1);
for k=1:nDates, notes{k}=''; end

fprintf('\nRolling forecasts: %d date(s) from %s to %s, horizon %d, area %d, %s %s %s (%s).\n', ...
    nDates,datestr(dateNumbers(1),'mm-dd-yyyy'),datestr(dateNumbers(end),'mm-dd-yyyy'),forecastingperiod,outbreakx,cadtemporal,caddisease,datatype,cadregion);

scores=zeros(0,8);
for k=1:nDates

    dateString=datestr(dateNumbers(k),'mm-dd-yyyy');
    timeStart=tic;
    fprintf('\n===== Rolling forecast %d of %d: data to %s =====\n',k,nDates,dateString);

    % ---- later data
    if rowsAtDate(k)+forecastingperiod>nLatest
        notes{k}=sprintf('skipped: the forecast needs data up to row %d and %s has %d rows',rowsAtDate(k)+forecastingperiod,latestName,nLatest);
        warning('SubEpiPredict:RollingDateSkipped','Forecast date %s %s.',dateString,notes{k});
        continue
    end

    % ---- input file of this date
    dateName=strcat(filePrefix,dateString,'.txt');
    if exist(strcat('./input/',dateName),'file')==2
        existing=load(strcat('./input/',dateName));
        if size(existing,1)~=rowsAtDate(k) || size(existing,2)~=size(latest,2)
            notes{k}=sprintf(['skipped: the file ./input/%s has %d row(s) and %d column(s); a series that starts on %s needs ', ...
                '%d row(s) on that date and %d column(s)'],dateName,size(existing,1),size(existing,2), ...
                datestr(datenum(datevecfirst1),'mm-dd-yyyy'),rowsAtDate(k),size(latest,2));
            warning('SubEpiPredict:RollingDateSkipped','Forecast date %s %s.',dateString,notes{k});
            continue
        end
        if isequal(existing(:,outbreakx),latest(1:rowsAtDate(k),outbreakx))
            status(k,3)=2;
        else
            status(k,3)=3;
            nDifferent=sum(existing(:,outbreakx)~=latest(1:rowsAtDate(k),outbreakx));
            notes{k}=sprintf('the existing file of this date differs from the most recent file in %d of %d values (the fit uses the existing file)',nDifferent,rowsAtDate(k));
            fprintf('Data: ./input/%s exists and %s.\n',dateName,notes{k});
        end
    else
        fid=fopen(strcat('./input/',dateName),'w');
        if fid<0
            error('SubEpiPredict:RollingInput','The file ./input/%s could not be written.',dateName);
        end
        fprintf(fid,[repmat('%.15g\t',1,size(latest,2)-1) '%.15g\n'],latest(1:rowsAtDate(k),:)');
        fclose(fid);
        status(k,3)=1;
        fprintf('Data: ./input/%s written (the first %d rows of %s).\n',dateName,rowsAtDate(k),latestName);
    end

    figuresBefore=findall(0,'type','figure');
    try

        % ---- fit, unless a fit of this date with the same B and start points is saved
        fitName=strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(nSubepidemics), ...
            '-onsetfixed-',num2str(onsetFixed),'-smoothing-',num2str(smoothing),'-', ...
            cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',dateString, ...
            '-flag1-',num2str(flag1(1)),'-method-',num2str(methodFit),'-dist-',num2str(distFit), ...
            '-calibrationperiod-',num2str(calibrationLength),'-rank-1.mat');
        reuse=false;
        if ~refit && exist(fitName,'file')==2
            % The file of rank 1 holds the driver's variables at the time it was
            % saved. The fit is used again if it was made with the same B, start
            % points and number of top models, and if it is complete: the file of
            % its last rank exists (an interrupted fit leaves the first ranks only).
            inFile=whos('-file',fitName);
            wanted={'M','numstartpoints','topmodelsx','topmodelsRequested','fitSignature'};
            wanted=wanted(ismember(wanted,{inFile.name}));
            saved=load(fitName,wanted{:});
            sameSettings=isfield(saved,'M') && isfield(saved,'numstartpoints') && isequal(saved.M,B) && isequal(saved.numstartpoints,numstartpoints);
            if sameSettings && isfield(saved,'topmodelsRequested')
                sameSettings=isequal(saved.topmodelsRequested,topmodelsSetting);
            end
            if isfield(saved,'topmodelsx')
                lastRank=saved.topmodelsx;
            else
                lastRank=topmodelsSetting;
            end
            complete=exist(strcat(fitName(1:end-5),num2str(lastRank),'.mat'),'file')==2;
            % Revision m1: the fit must also have been made by the current fit
            % revision with the same fit settings (subepidemicFitSignature.m).
            if isfield(saved,'fitSignature')
                [sameSignature,signatureDifferences]=subepidemicFitSignature(saved.fitSignature);
            else
                [sameSignature,signatureDifferences]=subepidemicFitSignature([]);
            end
            clear saved
            if sameSettings && complete && sameSignature
                reuse=true;
            elseif sameSettings && complete
                fprintf('Fit: the saved fit of this date is not reused (%s); the date is fitted again.\n',strjoin(signatureDifferences,', '));
            elseif ~sameSettings
                fprintf('Fit: a saved fit of this date was made with other values of B, numstartpoints or topmodelsx; the date is fitted again.\n');
            else
                fprintf('Fit: the saved fit of this date is incomplete (the file of rank %d is missing); the date is fitted again.\n',lastRank);
            end
        end
        nWarnings=0;
        consoleText='';
        if reuse
            status(k,4)=2;
            fprintf('Fit: the saved fit of this date is used (give ''refit'',true to fit every date again).\n');
        else
            % ---- revision m: keep a copy of the fit files this fit replaces
            if keepReplaced
                screenName=strrep(strrep(fitName,'modifiedLogisticPatch-ensem-','ABC-ensem-'),'-rank-1.mat','.mat');
                existingFiles=[dir(strcat(fitName(1:end-5),'*.mat'));dir(screenName)];
                nCopied=0;
                for q=1:numel(existingFiles)
                    if exist('./output/replaced-by-rolling','dir')~=7
                        mkdir('./output/replaced-by-rolling');
                    end
                    copyName=fullfile('./output/replaced-by-rolling',existingFiles(q).name);
                    if exist(copyName,'file')~=2
                        copyfile(fullfile('./output',existingFiles(q).name),copyName);
                        nCopied=nCopied+1;
                    end
                end
                if nCopied>0
                    fprintf('Fit: %d fit file(s) of this date already in ./output were copied to ./output/replaced-by-rolling/ before the new fit.\n',nCopied);
                end
            end
            if quiet
                consoleText=evalc('Run_Fit_subepidemicFramework(outbreakx,dateString);');
            else
                Run_Fit_subepidemicFramework(outbreakx,dateString);
            end
            status(k,4)=1;
            nWarnings=nWarnings+numel(subepidemicWarningLog('get'));
            if quiet
                fprintf('Fit: new (%d kind(s) of warning).\n',nWarnings);
            end
        end

        % ---- forecast and scores
        byHorizon=[];
        if quiet
            consoleText=[consoleText evalc('[~,~,~,byHorizon]=plotForecast_subepidemicFramework(outbreakx,dateString,forecastingperiod,weight_type1);')];
            logName=strcat('./output/rolling-log-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',dateString,'.txt');
            fid=fopen(logName,'w');
            if fid>0
                fwrite(fid,consoleText);
                fclose(fid);
            end
        else
            [~,~,~,byHorizon]=plotForecast_subepidemicFramework(outbreakx,dateString,forecastingperiod,weight_type1);
        end
        nWarningsForecast=numel(subepidemicWarningLog('get'));
        nWarnings=nWarnings+nWarningsForecast;
        status(k,6)=nWarnings;
        if quiet
            scoreText='no scores';
            if ~isempty(byHorizon)
                rowBase=byHorizon(byHorizon(:,1)==0 & byHorizon(:,2)==forecastingperiod,:);
                rowOne=byHorizon(byHorizon(:,1)==1 & byHorizon(:,2)==forecastingperiod,:);
                if ~isempty(rowBase) && ~isempty(rowOne)
                    scoreText=sprintf('WIS over %d time unit(s): ranked 1 %.2f, baseline %.2f, ratio %.2f',forecastingperiod,rowOne(1,6),rowBase(1,6),rowOne(1,6)/rowBase(1,6));
                end
            end
            fprintf('Forecast: %s (%d kind(s) of warning). %.0f s for this date.\nLog: %s\n',scoreText,nWarningsForecast,toc(timeStart),logName);
        end
        if isempty(byHorizon)
            notes{k}=strtrim([notes{k} ' no scores were returned']);
        else
            % WIS of the baseline at the same horizon, next to every row
            baseRows=byHorizon(byHorizon(:,1)==0,:);
            baseWIS=nan(forecastingperiod,1);
            baseWIS(baseRows(:,2))=baseRows(:,6);
            scores=[scores;[dateNumbers(k)+zeros(size(byHorizon,1),1) byHorizon(:,1:6) baseWIS(byHorizon(:,2))]]; %#ok<AGROW>
            status(k,5)=1;
        end

    catch problem
        notes{k}=sprintf('failed: %s',problem.message);
        warning('SubEpiPredict:RollingDateFailed','Forecast date %s failed and is left out: %s',dateString,problem.message);
    end

    if closeFigures
        figuresNow=findall(0,'type','figure');
        if ~isempty(figuresNow)
            if isempty(figuresBefore)
                close(figuresNow);
            else
                close(figuresNow(~ismember(figuresNow,figuresBefore)));
            end
        end
    end
    status(k,7)=toc(timeStart);

end

% ----------------------------------------------------------- summary
[summary,summaryColumns]=subepidemicRollingSummary(scores,forecastingperiod);

R=struct();
R.dates=cellstr(datestr(dateNumbers,'mm-dd-yyyy'));
R.status=status;
R.statusColumns={'date number','data rows','input file (1 written, 2 existing, 3 existing and different)','fit (1 new, 2 reused)','scored','kinds of warning','seconds'};
R.notes=notes;
R.scores=scores;
R.scoresColumns={'date number','model','horizon','MAE','MSE','Coverage 95%PI','WIS','WIS of the baseline'};
R.summary=summary;
R.summaryColumns=summaryColumns;
R.horizon=forecastingperiod;
R.scoresFile='';
R.summaryFile='';

nScored=sum(status(:,5)==1);
if nScored==0
    fprintf('\nRolling forecasts: no date was scored.\n');
    for k=1:nDates
        fprintf('  %s  %s\n',R.dates{k},notes{k});
    end
    return
end

% files, named by the first and the last scored date
fileTag=strcat('onsetfixed-',num2str(onsetFixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(methodFit),'-dist-',num2str(distFit), ...
    '-horizon-',num2str(forecastingperiod),'-weight_type-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion, ...
    '-area-',num2str(outbreakx),'-',datestr(min(scores(:,1)),'mm-dd-yyyy'),'-to-',datestr(max(scores(:,1)),'mm-dd-yyyy'),'.csv');
dv=datevec(scores(:,1));
T=array2table([dv(:,1:3) scores(:,2:8)]);
T.Properties.VariableNames(1:10)={'year','month','day','model','horizon','MAE','MSE','Coverage 95%PI','WIS','WIS of the baseline'};
R.scoresFile=strcat('./output/rolling-forecast-scores-',fileTag);
writetable(T,R.scoresFile);
T=array2table(summary);
T.Properties.VariableNames(1:numel(summaryColumns))=summaryColumns;
R.summaryFile=strcat('./output/rolling-forecast-summary-',fileTag);
writetable(T,R.summaryFile);

% table on the screen
if forecastingperiod<=6
    shown=1:forecastingperiod;
elseif DT==1 && forecastingperiod>=14
    shown=unique([7:7:forecastingperiod forecastingperiod]);
else
    shown=unique(round(linspace(1,forecastingperiod,5)));
end
models=unique(summary(:,1))';
fprintf('\nRolling forecasts: %d of %d date(s) scored, %s to %s, horizon %d.\n',nScored,nDates,R.dates{1},R.dates{end},forecastingperiod);
blocks={'WIS over the first h time units, relative to the baseline (ratio of means; below 1 = better than the baseline)',7; ...
        'The same, typical ratio (geometric mean over the dates)',8; ...
        'Share of dates with a WIS over the first h time units below the baseline''s',9; ...
        'WIS at step h alone, relative to the baseline (ratio of means)',12; ...
        'Coverage of the 95% prediction interval over the first h time units (%; nominal 95)',5};
for b=1:size(blocks,1)
    fprintf('\n%s\n',blocks{b,1});
    fprintf('  %-14s %6s',' ','dates');
    for h=shown
        fprintf(' %8s',sprintf('h=%d',h));
    end
    fprintf('\n');
    for m=models
        if m==0 && any(blocks{b,2}==[7 8 9 12])
            continue
        end
        if m==0
            modelLabel='baseline';
        elseif m>100
            modelLabel=sprintf('Ensemble(%d)',m-100);
        else
            modelLabel=sprintf('ranked %d',m);
        end
        fprintf('  %-14s',modelLabel);
        rowsModel=summary(summary(:,1)==m,:);
        fprintf(' %6d',max(rowsModel(:,3)));
        for h=shown
            value=rowsModel(rowsModel(:,2)==h,blocks{b,2});
            if isempty(value), value=NaN; end
            if blocks{b,2}==5
                fprintf(' %8.1f',value);
            else
                fprintf(' %8.2f',value);
            end
        end
        fprintf('\n');
    end
end
fprintf('\nFiles: %s\n       %s\n',R.scoresFile,R.summaryFile);
if any(status(:,5)~=1) || any(status(:,3)==3)
    fprintf('Notes:\n');
    for k=1:nDates
        if ~isempty(notes{k})
            fprintf('  %s  %s\n',R.dates{k},notes{k});
        end
    end
end
if any(status(:,6)>0)
    fprintf('Dates with toolbox warnings (see the list printed after each run above): %s\n',strjoin(R.dates(status(:,6)>0)',', '));
end
scoredDates=dateNumbers(status(:,5)==1);
if DT==1 && nScored>1 && all(mod(scoredDates,7)==mod(scoredDates(1),7))
    fprintf(['All scored dates fall on the same day of the week: with a weekly reporting cycle in the data, the scores ', ...
        'at a single step h mix the horizon with the day of the week. Compare steps 7 days apart, or use the scores over the first 7, 14, ... days.\n']);
end
if nScored<10
    fprintf('Only %d date(s): differences between the models and the baseline are weak evidence at this number.\n',nScored);
end

end
