%% generate short-term forecasts using best fitting models and derive ensemble model

clear
close all

% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global invasions
global timeinvasions
global Cinvasions
global npatches_fixed
global onset_fixed

global method1 dist1 factor1

global smoothfactor1

global calibrationperiod1

% <============================================================================>
% <================== Load the parameter values ===============================>
% <============================================================================>

[outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datafilename1_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP]=options

[getperformance_INP, deletetempfiles_INP, forecastingperiod_INP, printscreen1_INP, weight_type1_INP]=options_forecast


% <============================================================================>
% <================================ Dataset ===================================>
% <============================================================================>

outbreakx=outbreakx_INP;

caddate1=caddate1_INP;

cadregion=cadregion_INP;

caddisease=caddisease_INP;

datatype=datatype_INP;

datevecfirst1=datevecfirst1_INP;

datevecend1=datevecend1_INP;


DT=DT_INP; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
end

% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

%smoothfactor1=7; % <smoothfactor1>-day rolling average smoothing of the case series

%calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)

% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

%method1=0;  % Type of estimation method: 0 = LSQ

%dist1=0; % Normnal distribution to model error structure

M=M_INP; %number of bootstrap realizations to generate uncertainty estimates

% <==============================================================================>
% <========================= Growth model ==========================================>
% <==============================================================================>

%npatches_fixed=2; % maximum number of subepidemics considered in epidemic trajectory fit

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flagx=flag1_INP; % Sequence of subepidemic growth models considered in epidemic trajectory

%flagx=[RICH RICH]

% <===============================================================================================>
% <============= Number of best fitting models used to generate ensemble model ===================>
% <===============================================================================================>

topmodels1=1:topmodelsx_INP;

if npatches_fixed==1
    topmodels1=1;
end

factors=factor(length(topmodels1));
if length(factors)==1
    rows=factors;
    cols=1;
    
elseif length(factors)==3
    rows=factors(1)*factors(2);
    cols=factors(3);
else
    rows=factors(1);
    cols=factors(2);
end

            
% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

deletetempfiles=deletetempfiles_INP; %flag or indicator variable (1/0) to delete Forecast..mat files after use

forecastingperiod=forecastingperiod_INP; %forecast horizon (number of data points ahead)

%caddatex=[2020 04 20];
%ndays=2;
printscreen1=printscreen1_INP;  % print plots with the results

% <==============================================================================>
% <====================== weighting scheme for ensemble model ============================>
% <==============================================================================>

weight_type1=weight_type1_INP; % -1= equally weighted from the top models, 0=based on AICc, 1= based on relative likelihood (Akaike weights), 2=based on WISC during calibration, 3=based on WISF during forecasting performance at previous time period (week)

WISC_hash=zeros(length(topmodels1),1); % vector that saves the WISC based on calibration to be used with weight_type1=2

WISF_hash=zeros(length(topmodels1),200); % vector that saves the WISF based on prior forecasting performance to be used with weight_type1=2


% <=================================================================================>
% <========== Initialize variables to store forecast metrics across models =============================>
% <=================================================================================>

RMSECSS=[];
MSECSS=[];
MAECSS=[];
PICSS=[];
MISCSS=[];
RMSEFSS=[];
MSEFSS=[];
MAEFSS=[];
PIFSS=[];
MISFSS=[];

WISCSS=[];
WISFSS=[];

quantilescs=[];
quantilesfs=[];

% <==================================================================================================>
% <========== plot short-term forecasts for top-ranking subepidemic models and the ensemble models =========================>
% <==================================================================================================>

forecasts_best=[];
forecasts_ENS2=[];
forecasts_ENS3=[];
forecasts_ENS4=[];

%for run_id=0:1:52*ndays-1
for run_id=-1
%for run_id=0:1:97
    %for run_id=86
    %for run_id=0
    
    cc1=1;
    
    close all
    
    %i=(run_id)*30+1;
    %ARIMA_mean1=ARIMAforecasts(i:1:i+forecastingperiod-1,1);
    %ARIMA_lb1=ARIMAforecasts(i:1:i+forecastingperiod-1,11);
    %ARIMA_ub1=ARIMAforecasts(i:1:i+forecastingperiod-1,end-1);
    
    run_id
    
    if run_id==-1
        %outbreakx=52;
        
        run_id=0;
        
        %caddate1='11-15-21';
        %caddate1='06-22-20';
        %caddate1='05-31-20';
        
        %caddate1='04-20-20';
        %caddate1='05-11-20'; % for paper practical use

        %caddate1='06-29-20';
        %caddate1='07-20-20';

        %caddate1='09-28-20';
        %caddate1='03-22-21';
        
        %caddate1='01-24-22';
        
        cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);
       
        
    else
        
        %
        
        datenum1=datenum(caddatex);
        
        %state_id = rem(run_id,52)+1;
        state_id = 52;
        
        nm=fix(run_id/52);
        %date_id = datetime(caddatex) + caldays(nm);
        date_id = datetime(caddatex) + run_id*7;
        
        date=datestr(date_id,'mm-dd-yy');
        
        outbreakx=state_id;
        
        caddate1=date;
        
        cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);
                
        
        %
    end
    
    
    for rankx=topmodels1
        
        rankx
        
        npatches=npatches_fixed;
        
        %caddate1='04-05-21';
        
        % <========================================================================================>
        % <================================ Load model results ==========================================>
        % <========================================================================================>
        
        
        load (strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flagx(1)),'-flag1-',num2str(flagx(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rankx),'.mat'))
        
        
        d_hat=1;
        
        cc3=1;
        
        % <========================================================================================>
        % <================================ Compute short-term forecast ===================================>
        % <========================================================================================>
        
        timevect=(data1(:,1))*DT;
        
        timevect2=(0:t_window(end)-1+forecastingperiod)*DT;
        
        % vector to store forecast mean curves
        curvesforecasts1=[];
        
        % vector to store forecast prediction curves
        curvesforecasts2=[];
        
        color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];
        
        
        if printscreen1
            figure(10)
            %subplot(1,length(topmodels1),cc1)
            
            subplot(rows,cols,cc1)
            
        end
        
        % generate forecast curves from each bootstrap realization
        for realization=1:M
            
            rs_hat=Phatss(realization,1:npatches);
            ps_hat=Phatss(realization,npatches+1:2*npatches);
            as_hat=Phatss(realization,2*npatches+1:3*npatches);
            Ks_hat=Phatss(realization,3*npatches+1:4*npatches);
            
            
            IC=zeros(npatches,1);
            
            
            if method1>=3
                factor1=Phatss(realization,end-1);
            end
            
            if method1==5
                d_hat=Phatss(realization,end);
            end
            
            
            if onset_thr>0
                IC(1,1)=data1(1,2);
                IC(2:end,1)=1;
                
                invasions=zeros(npatches,1);
                timeinvasions=zeros(npatches,1);
                Cinvasions=zeros(npatches,1);
                
                invasions(1)=1;
                timeinvasions(1)=0;
                Cinvasions(1)=0;
            else
                IC(1:end,1)=data1(1,2)./length(IC(1:end,1));
                
                invasions=zeros(npatches,1);
                timeinvasions=zeros(npatches,1);
                Cinvasions=zeros(npatches,1);
                
                invasions(1:end)=1;
                timeinvasions(1:end)=0;
                Cinvasions(1:end)=0;
            end
            
            
            [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect2,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1);
            
            
            for j=1:npatches
                
                incidence1=[x(1,j);diff(x(:,j))];
                
                if printscreen1
                    plot(timevect2,incidence1,color1(j,:))
                    hold on
                end
                
            end
            
            y=sum(x,2);
            
            totinc=[y(1,1);diff(y(:,1))];
            
            totinc(1)=totinc(1)-(npatches-1);
            
            bestfit=totinc;
            
            gray1=gray(10);
            
            if printscreen1
                plot(timevect2,totinc,'color',gray1(7,:))
                
            end
            
            
            curvesforecasts1=[curvesforecasts1 totinc];
            
            forecasts2=AddPoissonError(cumsum(totinc),20,dist1,factor1,d_hat);
            
            curvesforecasts2=[curvesforecasts2 forecasts2];

        end


        [quantilesc,quantilesf]=computeQuantiles(data1,curvesforecasts2,forecastingperiod);

        quantilescs=[quantilescs;quantilesc];

        quantilesfs=[quantilesfs;quantilesf];


        % <==============================================================================================>
        % <================================ Plot short-term forecast ============================================>
        % <==============================================================================================>
        
        
        if printscreen1
            
            title(strcat(num2ordinal(rank1),' Ranked Model'))
            
            line1=plot(data1(:,1)*DT,data1(:,2),'ko')
            set(line1,'LineWidth',2)
            
            
            %title(getUSstateName(outbreak2))
            
            
            axis([0 length(timevect2)-1 0 max(data1(:,2))*2])
            
            
            line2=[timevect(end) 0;timevect(end) max(data1(:,2))*2];
            
            
            line1=plot(line2(:,1),line2(:,2),'k--')
            set(line1,'LineWidth',2)
            
            
            %wave=[2020 2 27 2020 9 08];
            
            %caddate1=caddate1(6:end);
            
            datenum1=datenum([str2num(caddate1(7:8))+2000 str2num(caddate1(1:2)) str2num(caddate1(4:5))]);
            
            datevec1=datevec(datenum1+forecastingperiod);
            
            wave=[datevecfirst1 datevec1(1:3)];
            
            
            
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
            
            
            line1=plot(line2(:,1),line2(:,2),'k--')
            set(line1,'LineWidth',2)
            
            ylabel(strcat(caddisease,{' '},datatype))
            
            %title(strcat('Sub-epidemic Model Forecast',{' '},getUSstateName(outbreakx),{' '},'- Reported by',{' '},caddate1))

            set(gca,'FontSize',24)
            set(gcf,'color','white')

        end


        LB1=quantile(curvesforecasts2',0.025);
        LB1=(LB1>=0).*LB1;

        UB1=quantile(curvesforecasts2',0.975);
        UB1=(UB1>=0).*UB1;

        if rank1==1 %top model
            %store forecast for best model
            forecasts_best=[forecasts_best;[median(curvesforecasts2(end-forecastingperiod+1:end,:),2) LB1(end-forecastingperiod+1:end)' UB1(end-forecastingperiod+1:end)']];
        end


        % <=============================================================================================>
        % <============================== Save file with forecast ======================================>
        % <=============================================================================================>

        forecastdata=[str2num(datestr((datenumIni:1:datenumEnd)','mm')) str2num(datestr((datenumIni:1:datenumEnd)','dd')) [data1(:,2);zeros(forecastingperiod,1)+NaN] median(curvesforecasts2,2) LB1' UB1'];

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:6) = {'month','day','data','median','LB','UB'};
        writetable(T,strcat('ranked-',cadregion,'-',caddate1,'-',num2str(rank1),'.csv'))


    
        if rank1==1 %top model
            %store forecast for top-ranking subepidemic model
            forecasts_best=[forecasts_best;[median(curvesforecasts2(end-forecastingperiod+1:end,:),2) LB1(end-forecastingperiod+1:end)' UB1(end-forecastingperiod+1:end)']];
        end
        
        if printscreen1
            
            figure(11)
            %subplot(1,length(topmodels1),cc1)
            subplot(rows,cols,cc1)
            
            
            hold on
            
            h=area(timevect2',[LB1' UB1'-LB1'])
            hold on
            
            h(1).FaceColor = [1 1 1];
            h(2).FaceColor = [0.8 0.8 0.8];
            
            %line1=plot(timevect2,quantile(curvesforecasts2',0.5),'r-')
            
            line1=plot(timevect2,median(curvesforecasts2,2),'r-')
            
            set(line1,'LineWidth',2)
            
            if  1
                line1=plot(timevect2,quantile(curvesforecasts2',0.025),'k--')
                set(line1,'LineWidth',2)
                
                line1=plot(timevect2,quantile(curvesforecasts2',0.975),'k--')
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
            
            line1=plot(line2(:,1),line2(:,2),'k--')
            set(line1,'LineWidth',2)
            
            ylabel(strcat(caddisease,{' '},datatype))
            
            title(strcat(num2ordinal(rank1),' Ranked Model'))
            
            %title(strcat('Sub-epidemic Model Forecast-',{' '},getUSstateName(outbreakx),{' '},'- Reported by',{' '},caddate1))
            
            set(gca,'FontSize',24)
            set(gcf,'color','white')
            
        end
        
        % <=========================================================================================>
        % <================================ Save short-term forecast results ==================================>
        % <=========================================================================================>
        
        save(strcat('./output/Forecast-modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-forecastingperiod-',num2str(forecastingperiod),'-rank-',num2str(rankx),'.mat'),'curvesforecasts1','curvesforecasts2','datevecfirst1','datevecend1','timevect2','timelags','cadtemporal')
        
        % <=============================================================================================>
        % <=================== Plot data for the forecast period (if getperformance=1) ===================================>
        % <=============================================================================================>
        
        datenum1=datenum([str2num(caddate1(7:8))+2000 str2num(caddate1(1:2)) str2num(caddate1(4:5))]);
        
        datenum1=datenum1+1;
        
        
        if getperformance & forecastingperiod>0
            
            data2=getData(cadtemporal,datevecfirst1,datevecend1,datevec(datenum1),outbreak1,forecastingperiod);
                        
            timevect2=(data1(end,1)+1:(data1(end,1)+1+forecastingperiod-1))*DT;
                        
            if printscreen1
                
                line2=plot(timevect2,data2,'ro')
                set(line2,'LineWidth',2)
                
            end
        else
            
            timevect2=[];
            
            data2=[];
            
        end
        

        datalatest2=[data1;[timevect2' data2]];
        
        
        % <==================================================================================================>
        % <========== Get forecast performance metrics for the model (if getperformance=1) =====================================>
        % <==================================================================================================>


        if getperformance

            [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1]=computeforecastperformance(data1,datalatest2,curvesforecasts1,curvesforecasts2,forecastingperiod);

            [WISC,WISFS]=computeWIS(data1,datalatest2,curvesforecasts2,forecastingperiod);

            WISC_hash(rankx,1)=WISC;  %save WISC for the model to be used to weight the models in the ensemble in function getensemblesubepidemics

            WISF_hash(rankx,run_id+1)=WISFS(end,end);  %save WISF for the model to be used to weight the models in the ensemble in function getensemblesubepidemics
            
            % store metrics for calibration
            RMSECSS=[RMSECSS;[rankx outbreakx datenum(caddate1) RMSECS_model1(end,end)]];
            MSECSS=[MSECSS;[rankx outbreakx datenum(caddate1) MSECS_model1(end,end)]];
            MAECSS=[MAECSS;[rankx outbreakx datenum(caddate1) MAECS_model1(end,end)]];
            PICSS=[PICSS;[rankx outbreakx datenum(caddate1) PICS_model1(end,end)]];
            MISCSS=[MISCSS;[rankx outbreakx datenum(caddate1) MISCS_model1(end,end)]];
            
            WISCSS=[WISCSS;[rankx outbreakx datenum(caddate1) WISC(end,end)]];
            
            
            % store metrics for short-term forecasts
            if forecastingperiod>0
                
                RMSEFSS=[RMSEFSS;[rankx outbreakx datenum(caddate1) RMSEFS_model1(end,end)]];
                MSEFSS=[MSEFSS;[rankx outbreakx datenum(caddate1) MSEFS_model1(end,end)]];
                MAEFSS=[MAEFSS;[rankx outbreakx datenum(caddate1) MAEFS_model1(end,end)]];
                PIFSS=[PIFSS;[rankx outbreakx datenum(caddate1) PIFS_model1(end,end)]];
                MISFSS=[MISFSS;[rankx outbreakx datenum(caddate1) MISFS_model1(end,end)]];
                
                WISFSS=[WISFSS;[rankx outbreakx datenum(caddate1) WISFS(end,end)]];
                
            end
            
        end
        
        cc1=cc1+1;
        
        %end
        
        % <=============================================================================================>
        % <=================== Get performance metrics for ensemble model, Ensemble(K) ===============================>
        % <=============================================================================================>
        
        if rankx>1
            
            if run_id==0
                [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1 quantilesc quantilesf]=getensemblesubepidemicmodels(cadfilename2,datevecfirst1,npatches_fixed,onset_fixed,smoothfactor1,outbreakx,cadregion,caddate1,caddisease,datatype,flag1,method1,dist1,calibrationperiod1,1:rankx,forecastingperiod,getperformance,weight_type1,WISC_hash,WISC_hash,printscreen1);
            else
                [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 WISC RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1 WISFS forecast1 quantilesc quantilesf]=getensemblesubepidemicmodels(cadfilename2,datevecfirst1,npatches_fixed,onset_fixed,smoothfactor1,outbreakx,cadregion,caddate1,caddisease,datatype,flag1,method1,dist1,calibrationperiod1,1:rankx,forecastingperiod,getperformance,weight_type1,WISC_hash,WISF_hash(:,run_id),printscreen1);
            end

            quantilescs=[quantilescs;quantilesc];

            quantilesfs=[quantilesfs;quantilesf];


            %store calibration performance metrics
            RMSECSS=[RMSECSS;[100+rankx outbreakx datenum(caddate1) RMSECS_model1(end,end)]];
            MSECSS=[MSECSS;[100+rankx outbreakx datenum(caddate1) MSECS_model1(end,end)]];
            MAECSS=[MAECSS;[100+rankx outbreakx datenum(caddate1) MAECS_model1(end,end)]];
            PICSS=[PICSS;[100+rankx outbreakx datenum(caddate1) PICS_model1(end,end)]];
            MISCSS=[MISCSS;[100+rankx outbreakx datenum(caddate1) MISCS_model1(end,end)]];
            
            WISCSS=[WISCSS;[100+rankx outbreakx datenum(caddate1) WISC(end,end)]];
            
 
            if forecastingperiod>0
                
                %store metrics for short-term forecasts
                RMSEFSS=[RMSEFSS;[100+rankx outbreakx datenum(caddate1) RMSEFS_model1(end,end)]];
                MSEFSS=[MSEFSS;[100+rankx outbreakx datenum(caddate1) MSEFS_model1(end,end)]];
                MAEFSS=[MAEFSS;[100+rankx outbreakx datenum(caddate1) MAEFS_model1(end,end)]];
                PIFSS=[PIFSS;[100+rankx outbreakx datenum(caddate1) PIFS_model1(end,end)]];
                MISFSS=[MISFSS;[100+rankx outbreakx datenum(caddate1) MISFS_model1(end,end)]];
                
                WISFSS=[WISFSS;[100+rankx outbreakx datenum(caddate1) WISFS(end,end)]];
            end
            
            
            % store the first 3 ensemble forecasts including the 95% Prediction interval.
            switch rankx
                case 2
                    forecasts_ENS2=[forecasts_ENS2;[forecast1(end-forecastingperiod+1:end,:)]];
                    
                case 3
                    forecasts_ENS3=[forecasts_ENS3;[forecast1(end-forecastingperiod+1:end,:)]];
                    
                case 4
                    forecasts_ENS4=[forecasts_ENS4;[forecast1(end-forecastingperiod+1:end,:)]];
                    
            end
            
        end
        
    end
    
    if deletetempfiles %flag or indicator variable (1/0) to delete Forecast..mat files after use
        for j=topmodels1
            
            delete(strcat('./output/Forecast-modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-forecastingperiod-',num2str(forecastingperiod),'-rank-',num2str(j),'.mat'))
        end
    end
    
    %pause
    
end

if run_id>=0
    % save(strcat('./output/performStats-run_id-0-97-weight_type-',num2str(weight_type1),'-modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-dateini-',datestr(datenum(caddatex),'mm-dd-yy'),'-ndays-',num2str(ndays),'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-forecastingperiod-',num2str(forecastingperiod),'-topmodels-',num2str(topmodels1(end)),'.mat'))

end
