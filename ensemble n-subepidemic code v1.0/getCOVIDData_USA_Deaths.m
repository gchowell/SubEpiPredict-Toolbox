function data2=getCOVIDData_USA(cadtemporal,date1,outbreak1,forecastingperiod)


datevecfirst1=[2020 02 27]; % date corresponding to the first data point in time series data

filename1=strcat('./input/cumulative-',cadtemporal,'-coronavirus-deaths-USA-05-09-22.txt');

data=load(filename1);

%data=data(:,1:2:end);
dataprov=data';

data1=dataprov(outbreak1,:)';

data1=[data1(1);diff(data1)];

datenum1=datenum(date1)-datenum(datevecfirst1)+1

data2=data1(datenum1:datenum1+forecastingperiod-1);

%length(data2)

