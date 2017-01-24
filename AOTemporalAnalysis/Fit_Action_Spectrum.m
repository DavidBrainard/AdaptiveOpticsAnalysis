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

%% Create the fits.

for i=2:2:size(a510,2)
    
%     figure(1); 
%     semilogx(a480(:,i-1), a480(:,i),'b',...
%             x480(i/2), interptarget,'b*',...
%             a510(:,i-1), a510(:,i),'c',...
%             x510(i/2), interptarget,'c*',...
%             a550(:,i-1), a550(:,i),'g',...
%             x550(i/2), interptarget,'g*',...
%             a590(:,i-1), a590(:,i),'y',...
%             x590(i/2), interptarget,'y*',...
%             a675(:,i-1), a675(:,i),'r',...
%             x675(i/2), interptarget,'r*'); %hold on;

    irradiances  = [{a480(:,i-1)}, {a510(:,i-1)}, {a550(:,i-1)}, {a590(:,i-1)}, {a675(:,i-1)}];
    amplitudes = [{a480(:,i)}, {a510(:,i)}, {a550(:,i)}, {a590(:,i)}, {a675(:,i)}];    

    
    rel_shifts(i/2,:) = action_spect_fit(wavelengths, irradiances, amplitudes);
end



rel_shifts = 1./rel_shifts;

figure; hold on;
for i=1:length(rel_shifts)
    semilogy(wavelengths,rel_shifts(i,:));
end
semilogy(ciefunc(:,1), ciefunc(:,2),'k')
hold off;

figure; errorbar(wavelengths, mean(rel_shifts), std(rel_shifts)); hold on;
title('Mean Wavelength vs Normalized Irradiance'); xlabel('Wavelength'); ylabel('Normalized Irradiance');
plot(ciefunc(:,1), ciefunc(:,2),'k')
set(gca,'yscale','log')
axis([450 700 10^-3 10^1]);
hold off;