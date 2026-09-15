
global invasions
global timeinvasions
global Cinvasions

invasions(1)=1;
timeinvasions(1)=0;
Cinvasions(1)=0;

close all

legends=['A';'B';'C';'D'];

onset_thr=[100 400 3000 5000]

duration1=220;

npatches=3;

flag1=1;

rs=[0.18 0.18 0.18];

ps=[0.98 0.98 0.98];

as=1+zeros(npatches,1);

as=[];

Ks=[10000 5000 1000];

I0=1;

for cc1=1:4

    figure(1)
    subplot(2,2,cc1)

    curve=plot_nsubepidemic(flag1,rs,ps,as, Ks,npatches,onset_thr(cc1),I0,duration1)

    title(strcat('C_{thr}=',num2str(onset_thr(cc1))))

    text(5,max(curve)-5,strcat('\fontsize{24}',legends(cc1)))

end
