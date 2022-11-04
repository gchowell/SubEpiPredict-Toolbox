
function [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1 quantilesc quantilesf]=getensemblesubepidemics(cadfilename1,datevecfirst1,npatches_fixed,onset_fixed,smoothfactor1,outbreakx,cadregion,caddate1,flag1,method1,dist1,calibrationperiod1,topmodels1,forecastingperiod,getperformance,weight_type1,WISC_hash,WISF_hash,printscreen1)

load(strcat('./output/ABC-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename1,'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')

% remove repeated rows
[RMSES,index1]=unique(RMSES,'rows','stable');
PS=PS(index1,:);

%topmodels1=1:2;

switch weight_type1
    
    case -1 % equally weighted from top models models
        
        weights1=ones(length(topmodels1),1)./length(topmodels1);
        
    case 0 % based on AICc
        
        AICc_best=RMSES(topmodels1,3);
        
        weights1=(1./AICc_best)/(sum(1./AICc_best)); % weights based on AICc

    case 1 %based on relative likelihood or Akaike weights

        AICmin=RMSES(1,3);

        relativelik_i=exp((AICmin-RMSES(topmodels1,3))/2);

        weights1=relativelik_i./sum(relativelik_i);  % weights based on relative likelihood

    case 2 % based on WISC during calibration of the models

        weights1=(1./WISC_hash(topmodels1))/sum(1./WISC_hash(topmodels1));
 
    case 3 % based on WISF based on forecasting performance during previous week (time period)

        weights1=(1./WISF_hash(topmodels1))/sum(1./WISF_hash(topmodels1));

end

% weight_type1
% weights1
% sum(weights1)
% pause

%relativelik_i

%min1=min(AICc_best);
%weights1=exp(-(AICc_best-min1)/2)./(sum(exp(-(AICc_best-min1)/2)));


%%

curvesforecasts1ens=[];
curvesforecasts2ens=[];

'ensemble'
topmodels1

for rank1=topmodels1
    
    
    load(strcat('./output/Forecast-modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',cadfilename1,'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-forecastingperiod-',num2str(forecastingperiod),'-rank-',num2str(rank1),'.mat'))
    
    M1=length(curvesforecasts1(1,:));
    
    index1=datasample(1:M1,round(M1*weights1(rank1)),'Replace',false);
    
    if length(index1)>0
        curvesforecasts1ens=[curvesforecasts1ens curvesforecasts1(:,index1)];
    end
    
    %
    
    M2=length(curvesforecasts2(1,:));
    
    index2=datasample(1:M2,round(M2*weights1(rank1)),'Replace',false);
    
    if length(index2)>0
        curvesforecasts2ens=[curvesforecasts2ens curvesforecasts2(:,index2)];
    end
    
end


curvesforecasts1=curvesforecasts1ens;

curvesforecasts2=curvesforecasts2ens;

% store forecast curves
forecast1=[median(curvesforecasts2,2) quantile(curvesforecasts2',0.975)' quantile(curvesforecasts2',0.025)'];

[quantilesc,quantilesf]=computeQuantiles(data1,curvesforecasts2,forecastingperiod);

if printscreen1
    
    figure(400)
    subplot(1,3,topmodels1(end)-1)
    
    datenum1=datenum([str2num(caddate1(7:8))+2000 str2num(caddate1(1:2)) str2num(caddate1(4:5))]);
    
    datevec1=datevec(datenum1+forecastingperiod);
    
    
    wave=[datevecfirst1 datevec1(1:3)];
    
    
    hold on
    
    quantile(curvesforecasts2',0.025)
    
    LB1=quantile(curvesforecasts2',0.025);
    LB1=(LB1>=0).*LB1;

    UB1=quantile(curvesforecasts2',0.975);
    UB1=(UB1>=0).*UB1;

    size(LB1)
    size(timevect2)
    
    h=area(timevect2',[LB1' UB1'-LB1'])
    hold on
    
    h(1).FaceColor = [1 1 1];
    h(2).FaceColor = [0.8 0.8 0.8];
    
    %line1=plot(timevect2,quantile(curvesforecasts2',0.5),'r-')
    
    line1=plot(timevect2,median(curvesforecasts2,2),'r-')
    
    set(line1,'LineWidth',2)
    
    if  1
        line1=plot(timevect2,LB1,'k--')
        set(line1,'LineWidth',2)
        
        line1=plot(timevect2,UB1,'k--')
        set(line1,'LineWidth',2)
    end
    
    gray1=gray(10);
    
    % plot time series datalatest
    line1=plot(data1(:,1)*DT,data1(:,2),'ko')
    set(line1,'LineWidth',2)
    
    
    axis([0 length(timevect2)-1 0 max(quantile(curvesforecasts2',0.975))*1.2])
    
    line2=[timevect(end) 0;timevect(end) max(quantile(curvesforecasts2',0.975))*1.20];
    
    box on
    
    line1=plot(line2(:,1),line2(:,2),'k--')
    set(line1,'LineWidth',2)
    
    
    % plot dates in x axis
    'day='
    datenum1=datenum(wave(1:3))+timelags; % start of fall wave (reference date)
    datestr(datenum1)
    
    datenumIni=datenum1;
    datenumEnd=datenum(wave(4:6))
    
    dates1=datestr(datenumIni:1:datenumEnd,'mm-dd');
    
    set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
    set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));
    xticklabel_rotate;
    
    ylabel('COVID-19 deaths')
    
    %title(strcat('Ensemble Model Forecast -',{' '},getUSstateName(outbreakx),{' '},'- Reported by',{' '},caddate1))
    title(strcat('Ensemble(',num2str(topmodels1(end)),')'))

    set(gca,'FontSize',24)
    set(gcf,'color','white')

end

% <=============================================================================================>
% <============================== Save file with forecast ======================================>
% <=============================================================================================>
forecastdata=[str2num(datestr((datenumIni:1:datenumEnd)','mm')) str2num(datestr((datenumIni:1:datenumEnd)','dd')) [data1(:,2);zeros(forecastingperiod,1)+NaN] median(curvesforecasts2,2) LB1' UB1'];

T = array2table(forecastdata);
T.Properties.VariableNames(1:6) = {'month','day','data','median','LB','UB'};

writetable(T,strcat('Ensemble(',num2str(topmodels1(end)),')-',cadregion,'-',caddate1,'.csv'))


if getperformance
    
    % plot most recent data
    
    
    %forecastingperiod=20;
    
    datenum1=datenum([str2num(caddate1(7:8))+2000 str2num(caddate1(1:2)) str2num(caddate1(4:5))]);
    
    datenum1=datenum1+1;
    
    
    data2=getCOVIDData_USA_Deaths(cadtemporal,datenum1,outbreakx,forecastingperiod)
    
    timevect2=(data1(end,1)+1:(data1(end,1)+1+forecastingperiod-1))*DT;
    
    if printscreen1
        line2=plot(timevect2,data2,'ro')
        set(line2,'LineWidth',2)
    end
    
    datalatest2=[data1;[timevect2' data2]];
    
    %%  compute performance metrics
    
    [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1]=computeforecastperformance(data1,datalatest2,curvesforecasts1,curvesforecasts2,forecastingperiod);
    
    [WISC,WISFS]=computeWIS(data1,datalatest2,curvesforecasts2,forecastingperiod);

else
    
    RMSECS_model1=NaN; MSECS_model1=NaN; MAECS_model1=NaN;  PICS_model1=NaN; MISCS_model1=NaN; RMSEFS_model1=NaN; MSEFS_model1=NaN; MAEFS_model1=NaN; PIFS_model1=NaN; MISFS_model1=NaN;
    
    WISC=NaN;

    WISFS=NaN;


end
