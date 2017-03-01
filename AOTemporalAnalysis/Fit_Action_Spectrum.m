clear;
close all;

% a480 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_480_20170120.csv');
% a510 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_510_20161130.csv');
% a550 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_550_20161130.csv');
% a590 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_590_20161219.csv');
% a675 = dlmread('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/allsubs_manual_675_20161130.csv');

% wavelengths = [480 510 550 590 675];

[fname, pathname] = uigetfile('*.csv');

respdata = dlmread(fullfile(pathname,fname));

% The 1st column is the IDs.
IDs = unique( respdata(2:end,1) );
% Top row are the wavelengths.
wavelengths = unique( respdata(1,2:end) );
% Everything else is a response.


fit_type = 'sigmoid';

%% Load CIE data

ciefunc = dlmread( fullfile(pathname,'linCIE2008v2e_5.csv') );

%% Create the fits.
legendIDs = cell(length(IDs),1);
for i=1:length(IDs)
    [idrows, ~]=find(respdata==IDs(i));
    
    subdata = respdata(idrows,1:end);
    
    irradiances = cell(1,length(wavelengths));
    datavals = cell(1,length(wavelengths));
    
    for w=1:length(wavelengths)
        [~, wavecols]=find(respdata==wavelengths(w));
        % First column is the irradiance, second is the data value.
        
        % **There will never be a 0 irradiance case.
        irrdata = subdata(:,wavecols(1));
        valdata = subdata(:,wavecols(2));
        
        valdata = valdata(irrdata~=0);
        irrdata = irrdata(irrdata~=0);
        
        irradiances(w) = {irrdata};
        datavals(w) = {valdata};
        
    end
    
    rel_shifts(i,:) = action_spect_fit(wavelengths, irradiances, datavals, ciefunc, fit_type);
    legendIDs{i} = num2str(IDs(i));
    
    title( legendIDs{i} );%For getting the irradiance response as a function of subject
    saveas(gcf, fullfile(pathname, [legendIDs{i} 'irr_response.svg']));
end

% For getting the irradiance response as a function of wavelength
% for w=1:length(wavelengths)
%     figure(wavelengths(w));
%     legend(legendIDs);
%     saveas(gcf, fullfile(pathname, [num2str(wavelengths(w)) 'irr_response.svg']));
% end

figure; hold on;
for i=1:length(IDs)
    semilogy(wavelengths,rel_shifts(i,:));
end
semilogy(ciefunc(:,1), ciefunc(:,2),'k')
legend(legendIDs);
hold off;

figure; errorbar(wavelengths, mean(rel_shifts), std(rel_shifts)); hold on;
title('Mean Wavelength vs Normalized Irradiance'); xlabel('Wavelength'); ylabel('Action (relative)');
plot(ciefunc(:,1), ciefunc(:,2),'k')
set(gca,'yscale','log')
axis([450 700 10^-3 10^1]);
hold off;
% 
% for i=2:2:size(a510,2)
%     
% %     figure(1); 
% %     semilogx(a480(:,i-1), a480(:,i),'b',...
% %             x480(i/2), interptarget,'b*',...
% %             a510(:,i-1), a510(:,i),'c',...
% %             x510(i/2), interptarget,'c*',...
% %             a550(:,i-1), a550(:,i),'g',...
% %             x550(i/2), interptarget,'g*',...
% %             a590(:,i-1), a590(:,i),'y',...
% %             x590(i/2), interptarget,'y*',...
% %             a675(:,i-1), a675(:,i),'r',...
% %             x675(i/2), interptarget,'r*'); %hold on;
% 
%     irradiances  = [{a480(:,i-1)}, {a510(:,i-1)}, {a550(:,i-1)}, {a590(:,i-1)}, {a675(:,i-1)}];
%     datavals = [{a480(:,i)}, {a510(:,i)}, {a550(:,i)}, {a590(:,i)}, {a675(:,i)}];    
% 
%     
%     rel_shifts(i/2,:) = action_spect_fit(wavelengths, irradiances, datavals);
% end