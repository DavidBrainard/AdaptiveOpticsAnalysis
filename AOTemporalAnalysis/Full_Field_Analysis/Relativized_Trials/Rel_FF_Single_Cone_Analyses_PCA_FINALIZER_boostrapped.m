

clear;

load('450nW_bootstrapped.mat');

AvgResp_450 = Avg_Resp;
StdResp_450 = Std_Resp;
Avg_StddevResp_450 = Avg_StddevResp;
Std_StddevResp_450 = Std_StddevResp;
Avg_MedianResp_450 = Avg_MedianResp;
Std_MedianResp_450 = Std_MedianResp;

load('50nW_bootstrapped.mat');

AvgResp_50 = Avg_Resp;
StdResp_50 = Std_Resp;
Avg_StddevResp_50 = Avg_StddevResp;
Std_StddevResp_50 = Std_StddevResp;
Avg_MedianResp_50 = Avg_MedianResp;
Std_MedianResp_50 = Std_MedianResp;

load('0nW_bootstrapped.mat');

AvgResp_0 = Avg_Resp;
StdResp_0 = Std_Resp;
Avg_StddevResp_0 = Avg_StddevResp;
Std_StddevResp_0 = Std_StddevResp;
Avg_MedianResp_0 = Avg_MedianResp;
Std_MedianResp_0 = Std_MedianResp;

%% Average Response
AllAvgResp = [AvgResp_450; AvgResp_50; AvgResp_0];
AllStdResp = [StdResp_450; StdResp_50; StdResp_0];

RespReg = [AllAvgResp ones(size(AllAvgResp))]\AllStdResp;

figure(1); clf; hold on;
plot( AvgResp_450, StdResp_450,'.');
plot( AvgResp_50, StdResp_50,'.');
plot( AvgResp_0, StdResp_0,'.');

r=-2:.1:12;
plot(r, r*RespReg(1)+RespReg(2),'k');

axis equal; grid on; axis([-1 10 0 2.5]);xlabel('Mean Response'); ylabel('Standard Deviation');

%% Std dev Response separate
AllStdAvgResp = [Avg_StddevResp_450; Avg_StddevResp_50; Avg_StddevResp_0];
AllStdStdResp = [Std_StddevResp_450; Std_StddevResp_50; Std_StddevResp_0];

StdRespReg = [AllStdAvgResp ones(size(AllStdAvgResp))]\AllStdStdResp;

figure(2); clf; subplot(2,1,1); hold on;
plot( Avg_StddevResp_450, Std_StddevResp_450,'.');
plot( Avg_StddevResp_50, Std_StddevResp_50,'.');
plot( Avg_StddevResp_0, Std_StddevResp_0,'.');

r=-2:.1:12;
plot(r, r*StdRespReg(1)+StdRespReg(2),'k');

axis equal; grid on; axis([-1 10 0 2.5]);xlabel('Mean StdDev Response'); ylabel('Standard Deviation');

subplot(2,1,2);
AllMedAvgResp = [Avg_MedianResp_450; Avg_MedianResp_50; Avg_MedianResp_0];
AllMedStdResp = [Std_MedianResp_450; Std_MedianResp_50; Std_MedianResp_0];

MedRespReg = [AllMedAvgResp ones(size(AllMedAvgResp))]\AllMedStdResp;

subplot(2,1,2); hold on;
plot( Avg_MedianResp_450, Std_MedianResp_450,'.');
plot( Avg_MedianResp_50, Std_MedianResp_50,'.');
plot( Avg_MedianResp_0, Std_MedianResp_0,'.');

r=-2:.1:12;
plot(r, r*MedRespReg(1)+MedRespReg(2),'k');
hold off;
axis equal; grid on; axis([-1 10 0 2.5]);xlabel('Mean Median Response'); ylabel('Standard Deviation');
