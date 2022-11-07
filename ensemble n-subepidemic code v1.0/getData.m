function data2=getData(cadtemporal,datevecfirst1,datevecend1,date1,outbreak1,forecastingperiod)


filename1=strcat('./input/cumulative-',cadtemporal,'-coronavirus-deaths-USA-',datestr(datenum(datevecend1),'mm-dd-yy'),'.txt');

data=load(filename1);

dataprov=data';

data1=dataprov(outbreak1,:)';

data1=[data1(1);diff(data1)];

datenum1=datenum(date1)-datenum(datevecfirst1)+1

data2=data1(datenum1:datenum1+forecastingperiod-1);


