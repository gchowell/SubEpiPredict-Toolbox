function [summary,columns] = subepidemicRollingSummary(scores,horizon)
%SUBEPIDEMICROLLINGSUMMARY Average forecast scores over dates, by model and horizon.
%   [summary,columns] = subepidemicRollingSummary(scores,horizon)
%
%   scores   one row per date, model and horizon, as collected by
%            subepidemicRollingForecast.m (R.scores): date number, model,
%            horizon, MAE, MSE, coverage of the 95% PI (%), WIS, WIS of the
%            baseline. Model 0 is the baseline. The scores of horizon h are
%            averages over the first h time units of a forecast.
%   horizon  the largest horizon
%
%   summary  one row per model and horizon:
%     1  model
%     2  horizon h
%     3  number of dates with a score
%     4  MAE over the first h time units, mean over the dates
%     5  coverage of the 95% PI over the first h time units (%), mean over the dates
%     6  WIS over the first h time units, mean over the dates
%     7  WIS relative to the baseline, ratio of means: mean WIS of the model /
%        mean WIS of the baseline, over the dates where both have a value
%     8  WIS relative to the baseline, typical ratio: geometric mean over the
%        dates of (WIS of the model / WIS of the baseline)
%     9  share of dates on which the model's WIS is below the baseline's (a tie counts half)
%     10 MAE relative to the baseline, ratio of means
%     11 WIS at step h alone, mean over the dates
%     12 WIS at step h alone relative to the baseline, ratio of means
%     13 coverage of the 95% PI at step h alone (%), the share of dates on
%        which the observation at step h was inside the interval
%
%   The score at step h alone is recovered from the averages:
%   h*score(h) - (h-1)*score(h-1). Columns 7 and 12 weigh dates with large
%   counts more (the WIS grows with the counts); columns 8 and 9 weigh every
%   date the same. Added in revision j (2026-11-17).

columns={'model','horizon','dates','MAE','Coverage 95%PI','WIS','WIS relative to baseline (ratio of means)', ...
    'WIS relative to baseline (typical ratio)','share of dates better than baseline','MAE relative to baseline (ratio of means)', ...
    'WIS at this step','WIS at this step relative to baseline','Coverage 95%PI at this step'};
summary=zeros(0,13);
if isempty(scores)
    return
end

dates=unique(scores(:,1));
models=unique(scores(:,2));
nDates=numel(dates);

% scores of each model as date-by-horizon arrays (NaN where there is no row)
maeOf=cell(numel(models),1); covOf=maeOf; wisOf=maeOf;
for m=1:numel(models)
    maeOf{m}=nan(nDates,horizon); covOf{m}=nan(nDates,horizon); wisOf{m}=nan(nDates,horizon);
    rowsModel=scores(scores(:,2)==models(m),:);
    [~,dateIndex]=ismember(rowsModel(:,1),dates);
    for k=1:size(rowsModel,1)
        h=rowsModel(k,3);
        if h>=1 && h<=horizon
            maeOf{m}(dateIndex(k),h)=rowsModel(k,4);
            covOf{m}(dateIndex(k),h)=rowsModel(k,6);
            wisOf{m}(dateIndex(k),h)=rowsModel(k,7);
        end
    end
end
base=find(models==0,1);
steps=repmat(1:horizon,nDates,1);

for m=1:numel(models)
    maeM=maeOf{m}; covM=covOf{m}; wisM=wisOf{m};
    % scores at step h alone
    wisStep=steps.*wisM-(steps-1).*[zeros(nDates,1) wisM(:,1:end-1)];
    covStep=steps.*covM-(steps-1).*[zeros(nDates,1) covM(:,1:end-1)];
    if isempty(base)
        maeB=nan(nDates,horizon); wisB=maeB; wisStepB=maeB;
    else
        maeB=maeOf{base}; wisB=wisOf{base};
        wisStepB=steps.*wisB-(steps-1).*[zeros(nDates,1) wisB(:,1:end-1)];
    end
    for h=1:horizon
        row=nan(1,13);
        row(1)=models(m); row(2)=h;
        have=isfinite(maeM(:,h));
        row(3)=sum(have);
        if row(3)==0
            continue
        end
        row(4)=mean(maeM(have,h));
        haveCov=isfinite(covM(:,h));
        if any(haveCov), row(5)=mean(covM(haveCov,h)); end
        haveWIS=isfinite(wisM(:,h));
        if any(haveWIS), row(6)=mean(wisM(haveWIS,h)); end
        both=haveWIS & isfinite(wisB(:,h));
        if any(both)
            if sum(wisB(both,h))>0
                row(7)=sum(wisM(both,h))/sum(wisB(both,h));
            end
            positive=both & wisM(:,h)>0 & wisB(:,h)>0;
            if any(positive)
                row(8)=exp(mean(log(wisM(positive,h)./wisB(positive,h))));
            end
            if models(m)~=0
                row(9)=(sum(wisM(both,h)<wisB(both,h))+0.5*sum(wisM(both,h)==wisB(both,h)))/sum(both);
            end
        end
        bothMAE=have & isfinite(maeB(:,h));
        if any(bothMAE) && sum(maeB(bothMAE,h))>0
            row(10)=sum(maeM(bothMAE,h))/sum(maeB(bothMAE,h));
        end
        haveStep=isfinite(wisStep(:,h));
        if any(haveStep), row(11)=mean(wisStep(haveStep,h)); end
        bothStep=haveStep & isfinite(wisStepB(:,h));
        if any(bothStep) && sum(wisStepB(bothStep,h))>0
            row(12)=sum(wisStep(bothStep,h))/sum(wisStepB(bothStep,h));
        end
        haveCovStep=isfinite(covStep(:,h));
        if any(haveCovStep), row(13)=mean(covStep(haveCovStep,h)); end
        summary=[summary;row]; %#ok<AGROW>
    end
end

end
