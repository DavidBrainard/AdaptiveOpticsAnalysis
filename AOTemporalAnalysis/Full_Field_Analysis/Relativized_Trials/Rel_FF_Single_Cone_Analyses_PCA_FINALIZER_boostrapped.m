

clear all;

load('450nW_bootstrapped.mat');

AvgResp_450 = Avg_Resp;
StdResp_450 = Std_Resp;

load('50nW_bootstrapped.mat');

AvgResp_50 = Avg_Resp;
StdResp_50 = Std_Resp;

load('0nW_bootstrapped.mat');

AvgResp_0 = Avg_Resp;
StdResp_0 = Std_Resp;


AllAvgResp = [AvgResp_450; AvgResp_50; AvgResp_0];
AllStdResp = [StdResp_450; StdResp_50; StdResp_0];

RespReg = [AllAvgResp ones(size(AllAvgResp))]\AllStdResp;

figure(1); clf; hold on;
plot( AvgResp_450, StdResp_450,'.');
plot( AvgResp_50, StdResp_50,'.');
plot( AvgResp_0, StdResp_0,'.');

r=-2:.1:12;
plot(r, r*RespReg(1)+RespReg(2),'k');

axis([-1 10 0 2.5]);