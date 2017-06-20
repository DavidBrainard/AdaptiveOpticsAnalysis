% clear;
% close all;

%% Average data from 11-28-2016
% w510 = [390 1.444887;
%         140 1.16780;
%         17  0.32;
%         2.4 0.121];
% w550 = [390 1.5476;
%         140 1.28898;
%         17  0.46411;
%         2.4 0.079];
w510 = [390 1.7;
        140 1.304;
        17  0.243;
        2.4 0];   
w550 = [390 1.78;
        140 1.53;
        17  0.567;
        2.4 0];   
w680 = [18900 1.087;
        2260 0.847;
        266  0.271;
        37 0.031];
    
wavelengths = [510 550 680];

%% Load CIE data

ciefunc = dlmread('/local_data/Projects/AdaptiveOpticsAnalysis/AOTemporalAnalysis/linCIE2008v2e_5.csv');
    
%%
interptarget = [.6];

x510 = interp1(w510(2:3,2), w510(2:3,1), interptarget,'linear');
x550 = interp1(w550(2:3,2), w550(2:3,1), interptarget,'linear');
x680 = interp1(w680(2:3,2), w680(2:3,1), interptarget,'linear');


figure(10); 
%           semilogx(w510(:,1), w510(:,2),'b',...
%                     x510, interptarget,'b*',...
%                     w550(:,1), w550(:,2),'g',...
%                     x550, interptarget,'g*',...
%                     w680(:,1), w680(:,2),'r',...
%                     x680, interptarget,'r*');
errorbar(w550(:,1), w550(:,2),std(w550(:,1)),std(w550(:,2)));
set(gca,'yscale','log')
title('Amplitude vs Irradiance'); xlabel('Irradiance'); ylabel('Amplitude');

x510 = (1./x510);
x550 = (1./x550);
x680 = (1./x680);

% xe510 = (1./x510);
% xe550 = (1./x550);
% xe680 = (1./x680);

x510 = x510./x550;
x680 = x680./x550;
x550 = x550./x550;

% xe510 = xe510./x510;
% xe680 = xe680./x680;
% xe550 = xe550./x550;

for i=1:length(interptarget)
    figure(2); semilogy(wavelengths, [x510(i), x550(i), x680(i)] ); hold on;
    title('Wavelength vs Normalized Irradiance'); xlabel('Wavelength'); ylabel('Normalized Irradiance');
end
semilogy(ciefunc(:,1), ciefunc(:,2),'k')
hold off;

