

function [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1 quantilesc quantilesf]=getensemblesubepidemics(cumulative1,cadfilename1,datevecfirst1,npatches_fixed,onset_fixed,smoothfactor1,outbreakx,cadregion,caddate1,caddisease,datatype,flag1,method1,dist1,calibrationperiod1,topmodels1,forecastingperiod,getperformance,weight_type1,WISC_hash,WISF_hash,printscreen1,ensemblePanelGrid)

load(strcat('./output/ABC-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename1,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')


% remove any repeated rows
% [npatches onset_thr AICc];

[RMSES,index1]=unique(RMSES,'rows','stable');
PS=PS(index1,:);

%RMSES(:,1) %npatches
%RMSES(:,2) %onset_thr
%RMSES(:,3) % AICc


% Weights of the member models (subepidemicEnsembleWeights.m): -1 equal,
% 0 and 1 Akaike weights, 2 inverse WISC, 3 inverse WISF. Code 0 used to be
% (1/AICc)/sum(1/AICc), which is not a weight; it now gives Akaike weights and
% warns. The scores are read only for the code that uses them.
scoresC=[];
scoresF=[];
if weight_type1==2
    scoresC=WISC_hash(topmodels1);
elseif weight_type1==3
    scoresF=WISF_hash(topmodels1);
end
weights1=subepidemicEnsembleWeights(weight_type1,RMSES(topmodels1,3),scoresC,scoresF);

% weight_type1
% weights1
% sum(weights1)
% pause

%relativelik_i

%min1=min(AICc_best);
%weights1=exp(-(AICc_best-min1)/2)./(sum(exp(-(AICc_best-min1)/2)));


%%

% Revision of 2026-11-24. The ensemble is the mixture of its members. Member k
% with weight weights1(k) gives ALL its curves, and each of its M curves
% carries the weight weights1(k)/M; medians, intervals, quantiles and scores
% of the ensemble are weighted quantiles (subepidemicWeightedQuantile.m).
% Before, member k gave round(M*weights1(k)) curves drawn at random, so
% Ensemble(4) with equal weights used a quarter of the curves of each member,
% the rounding changed the weights (a member below 1/(2M) gave no curve), and
% the ensemble changed with the random seed. With equal weights and the same
% number of curves per member, the ensemble is now simply all curves together.
ensCurves1=[];   % model trajectories of the members, side by side
ensCurves2=[];   % predictive curves (trajectory plus noise) of the members
ensWeights1=[];  % one weight per column of ensCurves1
ensWeights2=[];  % one weight per column of ensCurves2

'ensemble';
topmodels1;

for rank1=topmodels1

    % Resetting the quantile array
    combinedQuantiles = [];

    load(strcat('./output/Forecast-modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename1,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-forecastingperiod-',num2str(forecastingperiod),'-rank-',num2str(rank1),'.mat'))

    M1=length(curvesforecasts1(1,:));

    M2=length(curvesforecasts2(1,:));

    if weights1(rank1)>0 && M1>0 && M2>0
        ensCurves1=[ensCurves1 curvesforecasts1];
        ensWeights1=[ensWeights1 repmat(weights1(rank1)/M1,1,M1)];

        ensCurves2=[ensCurves2 curvesforecasts2];
        ensWeights2=[ensWeights2 repmat(weights1(rank1)/M2,1,M2)];
    end

end

if isempty(ensCurves1) || isempty(ensCurves2)
    error('SubEpiPredict:EnsembleEmpty', ...
        'Ensemble of %d models: no member with a positive weight has forecast curves.',numel(topmodels1));
end

curvesforecasts1=ensCurves1;

curvesforecasts2=ensCurves2;

% the weights sum to 1 (they already do unless a member was left out above)
ensWeights1=ensWeights1/sum(ensWeights1);
ensWeights2=ensWeights2/sum(ensWeights2);

% Median and 95% limits of the ensemble. They are computed here, not in the
% plot section, because the forecast files below need them with or without plots.
ensLimits=subepidemicWeightedQuantile(curvesforecasts2,[0.025 0.5 0.975],ensWeights2);
ensMedian=ensLimits(:,2);

LB1=ensLimits(:,1)';
LB1=(LB1>=0).*LB1;

UB1=ensLimits(:,3)';
UB1=(UB1>=0).*UB1;

% store forecast curves
forecast1=[ensMedian ensLimits(:,3) ensLimits(:,1)];

[quantilesc,quantilesf]=computeQuantiles(data1,curvesforecasts2,forecastingperiod,ensWeights2);

% Names for quantile forecast table
quantNamesRanked = {'Q_0.010', 'Q_0.025', 'Q_0.050', 'Q_0.100', 'Q_0.150', 'Q_0.200', 'Q_0.250', 'Q_0.300', 'Q_0.350', 'Q_0.400', 'Q_0.450', 'Q_0.500', 'Q_0.550', 'Q_0.600', 'Q_0.650', 'Q_0.700', 'Q_0.750', 'Q_0.800', 'Q_0.850', 'Q_0.900', 'Q_0.950', 'Q_0.975', 'Q_0.990'};

% Quantile forecast array
combinedQuantiles = [quantilesc; quantilesf];
combinedQuantilesTable = array2table(combinedQuantiles, 'VariableNames', quantNamesRanked);

% Exporting the quantile forecasts
writetable(combinedQuantilesTable,strcat('./output/quantileTimes-Ensemble(',num2str(topmodels1(end)),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-weighttype-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


if printscreen1

    figure(500)

    % One panel per ensemble, Ensemble(2) in panel 1. The grid was fixed at
    % 1 by 3, so the ensemble of five models asked for a fourth panel and the
    % run stopped. The forecast script now passes the grid for all ensembles;
    % without it, three panels as before, or as many as this ensemble needs.
    if nargin<23 || isempty(ensemblePanelGrid)
        ensemblePanelGrid=[1 max(3,topmodels1(end)-1)];
    end
    if ensemblePanelGrid(1)*ensemblePanelGrid(2)<topmodels1(end)-1
        ensemblePanelGrid=[ceil((topmodels1(end)-1)/3) 3];
    end
    subplot(ensemblePanelGrid(1),ensemblePanelGrid(2),topmodels1(end)-1)

    datenum1=datenum([str2num(caddate1(7:10)) str2num(caddate1(1:2)) str2num(caddate1(4:5))]);

    datevec1=datevec(datenum1+forecastingperiod*DT);


    wave=[datevecfirst1 datevec1(1:3)];


    hold on

    size(LB1);
    size(timevect2);

    h=area(timevect2',[LB1' UB1'-LB1']);
    hold on

    h(1).FaceColor = [1 1 1];
    h(2).FaceColor = [0.8 0.8 0.8];

    %line1=plot(timevect2,quantile(curvesforecasts2',0.5),'r-')

    line1=plot(timevect2,ensMedian,'r-');

    set(line1,'LineWidth',2)

    if  1
        line1=plot(timevect2,LB1,'k--');
        set(line1,'LineWidth',2)

        line1=plot(timevect2,UB1,'k--');
        set(line1,'LineWidth',2)
    end

    gray1=gray(10);

    % plot time series datalatest
    line1=plot(data1(:,1),data1(:,2),'ko');
    set(line1,'LineWidth',2)


    axis([0 length(timevect2)-1 0 max(UB1)*1.2])

    line2=[timevect(end) 0;timevect(end) max(UB1)*1.20];

    box on

    line1=plot(line2(:,1),line2(:,2),'k--');
    set(line1,'LineWidth',2)


    % plot dates in x axis
    'day=';
    datenum1=datenum(wave(1:3))+timelags*DT; % start of fall wave (reference date)
    datestr(datenum1);

    datenumIni=datenum1;
    datenumEnd=datenum(wave(4:6));

    dates1=datestr(datenumIni:DT:datenumEnd,'mm-dd');

    if DT==1

        set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
        set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));

    elseif DT==7

        set(gca, 'XTick', 0:2:length(dates1(:,1))-1);
        set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:2:end,:)));

    elseif DT==365

        years1=wave(1)+timelags:wave(1)+timelags+length(dates1(:,1))-1;

        set(gca,'XTick',0:1:length(years1)-1);
        set(gca, 'XTickLabel', strcat('\fontsize{14}',num2str(years1')));

    end

    xticklabel_rotate;

    ylabel(strcat(caddisease,{' '},datatype))

    %title(strcat('Ensemble Model Forecast -',{' '},getUSstateName(outbreakx),{' '},'- Reported by',{' '},caddate1))
    title(strcat('Ensemble(',num2str(topmodels1(end)),')'))

    set(gca,'FontSize',GetAdjustedFontSize)
    set(gcf,'color','white')

end


% compute doubling times

% Revision 2026-11-16: meandoublingtime starts empty. It used to start as M zeros
% and the M values were appended after them, which halved the mean and put the
% lower limit at 0. Curves with no doubling give NaN and are left out.
% Revision of 2026-11-24: every curve of every member enters with its weight
% (ensWeights1), so the mean, the limits and the probability of the i-th
% doubling are weighted. The arrays are sized after the number of curves of
% the ensemble; they were sized after the last member loaded.
M1=length(curvesforecasts1(1,:));

meandoublingtime=zeros(M1,1)+NaN;

doublingtimess=zeros(30,M1)+NaN;

maxd=1;

for j=1:M1

    [tds,C0data,curve,doublingtimes]=getDoublingTimeCurve(max(curvesforecasts1(:,j),0),DT,0);

    doublingtimess(1:length(doublingtimes),j)=doublingtimes;

    if maxd<length(doublingtimes)
        maxd=length(doublingtimes);
    end

    meandoublingtime(j)=mean(doublingtimes);
end

doublingtimess=doublingtimess(1:maxd,1:M1);

seq_doublingtimes=[];

for j=1:maxd

    index1=find(~isnan(doublingtimess(j,:)));

    if isempty(index1)
        seq_doublingtimes=[seq_doublingtimes;[j NaN NaN NaN 0]];
    else
        weightsj=ensWeights1(index1);
        limitsj=subepidemicWeightedQuantile(doublingtimess(j,index1),[0.025 0.975],weightsj);
        seq_doublingtimes=[seq_doublingtimes;[j sum(weightsj.*doublingtimess(j,index1))/sum(weightsj) limitsj(1) limitsj(2) sum(weightsj)/sum(ensWeights1)]];
    end

end

seq_doublingtimes; % [ith doubling, mean, 95%CI LB, 95%CI UB, prob. i_th doubling]

% Mean doubling times
index1=find(~isnan(meandoublingtime));
if isempty(index1) % no curve doubled
    dmean=NaN;
    dLB=NaN;
    dUB=NaN;
else
    weightsj=ensWeights1(index1);
    dmean=sum(weightsj.*meandoublingtime(index1)')/sum(weightsj);
    limitsj=subepidemicWeightedQuantile(meandoublingtime(index1)',[0.025 0.975],weightsj);
    dLB=limitsj(1);
    dUB=limitsj(2);
end

param_doubling=[dmean dLB dUB]

% <=============================================================================================>
% <============================== Save file with doubling time estimates =======================>
% <=============================================================================================>

T = array2table(seq_doublingtimes);
T.Properties.VariableNames(1:5) = {'i_th doubling','db mean','db 95%CI LB','db 95% CI UB','prob. i_th doubling'};
writetable(T,strcat('./output/doublingTimes-Ensemble(',num2str(topmodels1(end)),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-weighttype-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


if getperformance

    % plot most recent data

    datenum1=datenum([str2num(caddate1(7:10)) str2num(caddate1(1:2)) str2num(caddate1(4:5))]);

    if (DT~=365)

        datenum1=datenum1+DT;

    end

    if getperformance
        data2=getData(cumulative1,cadtemporal,caddisease,datatype,cadregion,DT,datevecfirst1,datevecend1,datevec(datenum1),outbreakx,forecastingperiod*getperformance);

        timevect2=(data1(end,1)+1:(data1(end,1)+1+getperformance*forecastingperiod-1));

        if printscreen1
            line2=plot(timevect2,data2,'ro');
            set(line2,'LineWidth',2)
        end

        datalatest2=[data1;[timevect2' data2]];

        % <=============================================================================================>
        % <============================== Save file with forecast ======================================>
        % <=============================================================================================>
        forecastdata=[str2num(datestr((datenumIni:DT:datenumEnd)','yyyy')) str2num(datestr((datenumIni:DT:datenumEnd)','mm')) str2num(datestr((datenumIni:DT:datenumEnd)','dd')) [data1(:,2);data2] ensMedian LB1' UB1'];

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:7) = {'year','month','day','data','median','LB','UB'};

        writetable(T,strcat('./output/Ensemble(',num2str(topmodels1(end)),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-weighttype-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))
    else
        datalatest2=data1;

    end


    %%  compute performance metrics

    [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1]=computeforecastperformance(data1,datalatest2,curvesforecasts1,curvesforecasts2,forecastingperiod*getperformance,ensWeights1,ensWeights2);

    [WISC,WISFS]=computeWIS(data1,datalatest2,curvesforecasts2,forecastingperiod*getperformance,ensWeights2);

else

    % plot most recent data

    datenum1=datenum([str2num(caddate1(7:10)) str2num(caddate1(1:2)) str2num(caddate1(4:5))]);

    if (DT~=365)

        datenum1=datenum1+DT;

    end

    RMSECS_model1=NaN; MSECS_model1=NaN; MAECS_model1=NaN;  PICS_model1=NaN; MISCS_model1=NaN; RMSEFS_model1=NaN; MSEFS_model1=NaN; MAEFS_model1=NaN; PIFS_model1=NaN; MISFS_model1=NaN;

    WISC=NaN;

    WISFS=NaN;


    % <=============================================================================================>
    % <============================== Save file with forecast ======================================>
    % <=============================================================================================>
    forecastdata=[str2num(datestr((datenumIni:DT:datenumEnd)','yyyy')) str2num(datestr((datenumIni:DT:datenumEnd)','mm')) str2num(datestr((datenumIni:DT:datenumEnd)','dd')) [data1(:,2);zeros(forecastingperiod,1)+NaN] ensMedian LB1' UB1'];


    T = array2table(forecastdata);
    T.Properties.VariableNames(1:7) = {'year','month','day','data','median','LB','UB'};

    writetable(T,strcat('./output/Ensemble(',num2str(topmodels1(end)),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-weighttype-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


end
