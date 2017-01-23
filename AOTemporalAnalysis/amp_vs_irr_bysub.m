clear;
close all;

a480 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_480_20170120.csv');
a510 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_510_20161130.csv');
a550 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_550_20161130.csv');
a590 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_590_20161219.csv');
a675 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_675_20161130.csv');

wavelengths = [480 510 550 590 675];

%% Load CIE data

ciefunc = dlmread('/local_data/Projects/AdaptiveOpticsAnalysis/AOTemporalAnalysis/linCIE2008v2e_5.csv');
    
%%http://www.bbc.com/news/live/world-us-canada-38682565
interptarget = .25;
% i=6;
for i=2:2:size(a510,2)
    
    x480(i/2) = interp1(a480(:,i), a480(:,i-1), interptarget,'linear');    
    x510(i/2) = interp1(a510(:,i), a510(:,i-1), interptarget,'linear');
    x550(i/2) = interp1(a550(:,i), a550(:,i-1), interptarget,'linear');
    x590(i/2) = interp1(a590(:,i), a590(:,i-1), interptarget,'linear');
    x675(i/2) = interp1(a675(:,i), a675(:,i-1), interptarget,'linear');
    
    figure(1); 
    semilogx(a480(:,i-1), a480(:,i),'b',...
            x480(i/2), interptarget,'b*',...
            a510(:,i-1), a510(:,i),'c',...
            x510(i/2), interptarget,'c*',...
            a550(:,i-1), a550(:,i),'g',...
            x550(i/2), interptarget,'g*',...
            a590(:,i-1), a590(:,i),'y',...
            x590(i/2), interptarget,'y*',...
            a675(:,i-1), a675(:,i),'r',...
            x675(i/2), interptarget,'r*'); %hold on;

end
title('Amplitude vs Irradiance'); xlabel('Irradiance'); ylabel('Amplitude');
hold off;

% the550x = a550(:,1:2:end);
% the550y = a550(:,2:2:end);
% 
% figure(10);
% errorbar(mean(the550x,2), mean(the550y,2), std(the550y,[],2), std(the550y,[],2) );
% set(gca,'xscale','log')


x480 = (1./x480);
x510 = (1./x510);
x550 = (1./x550);
x590 = (1./x590);
x675 = (1./x675);

x480 = x480./x550;
x510 = x510./x550;
x675 = x675./x550;
x590 = x590./x550;
x550 = x550./x550;

x480(4) = x480(3);
% i=3;
for i=1:length(x510)
    figure(2); semilogy(wavelengths, [x480(i), x510(i), x550(i), x590(i), x675(i)] ); hold on;
    title('Wavelength vs Normalized Irradiance'); xlabel('Wavelength'); ylabel('Normalized Irradiance');
end
semilogy(ciefunc(:,1), ciefunc(:,2),'k')
hold off;


figure(3); errorbar(wavelengths, [mean(x480), mean(x510), mean(x550), mean(x590), mean(x675)], ...
                                 [std(x480),  std(x510),  std(x550),  std(x590),  std(x675)] ); hold on;
title('Mean Wavelength vs Normalized Irradiance'); xlabel('Wavelength'); ylabel('Normalized Irradiance');
plot(ciefunc(:,1), ciefunc(:,2),'k')
set(gca,'yscale','log')
hold off;

