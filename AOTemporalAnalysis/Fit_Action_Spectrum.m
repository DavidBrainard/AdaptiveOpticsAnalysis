clear;
% close all;

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
% ciefunc2 = dlmread( fullfile(pathname,'linss2_10e_5_M.csv') );

% temp for plot of luminosity vs response
wavecolors=['b' 'c', 'g', 'y','r'];
% isover9k = figure(9001); %ITS OVER NINE THOUSAAAAAAAAND

%% Create the fits.
legendIDs = cell(length(IDs),1);
luminosity = cell(1,length(wavelengths));
for i=1:length(IDs)
    [idrows, ~]=find(respdata==IDs(i));
    
    subdata = respdata(idrows,1:end);
    
    irradiances = cell(1,length(wavelengths));
    
    datavals = cell(1,length(wavelengths));
    
%     figure(i+100); hold on;
    for w=1:length(wavelengths)
        [~, wavecols]=find(respdata==wavelengths(w));
        % First column is the irradiance, second is the data value.
        
        % **There will never be a 0 irradiance case.
        irrdata = subdata(:,wavecols(1));
        valdata = subdata(:,wavecols(2));
        
        valdata = valdata(irrdata~=0);
        irrdata = irrdata(irrdata~=0);
        
        irradiances(w) = {irrdata};
%         luminosity{w} = [luminosity{w}; AOLightLevelConversions_Func(1, repmat(wavelengths(w),length(irradiances{w}),1), irradiances{w}/1000, false)];
        datavals(w) = {valdata};
        
%         plot(log10(luminosity{w}),datavals{w}, [wavecolors(w) 'o'],'MarkerFaceColor',wavecolors(w), 'MarkerSize',6);
        
    end
%     hold off;
%     title( num2str(IDs(i)) ); xlabel('Log candelas'); ylabel('Amplitude response');
    
%     saveas(gcf, fullfile(pathname, [num2str(IDs(i)) 'lum_response.png']));
    
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

figure(1); hold on;
for i=1:length(IDs)
    plot(wavelengths, log10(rel_shifts(i,:)),'.','MarkerSize',15);
end
plot(ciefunc(:,1), log10(ciefunc(:,2)),'k')
% set(gca,'yscale','log')
legend(legendIDs);
axis([450 700 -3 1]);
hold off;

figure; errorbar(wavelengths, mean( log10(rel_shifts) ), std(log10(rel_shifts))); hold on;
title('Mean Wavelength vs Normalized Irradiance'); xlabel('Wavelength'); ylabel('Action (relative)');
plot(ciefunc(:,1), log10(ciefunc(:,2)),'k')
% set(gca,'yscale','log')
axis([450 700 -3 1]);
hold off;
