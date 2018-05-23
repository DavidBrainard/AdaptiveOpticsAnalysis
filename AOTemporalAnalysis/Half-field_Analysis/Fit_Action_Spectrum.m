% Fit_Action_Spectrum
%
% This script is responsible for taking the summarized reflectance
% response data from each subject and wavelength, fitting each
% wavelength's reflectance response vs irradiance with
% a sigmoid, then determining each subject's individual action spectrum.
% It then outputs the average of all subjects' action spectra as well.
%
% Created by Robert F Cooper 01-24-2017
%
% The analyses performed in this script are from:
% Cooper RF, Tuten WS, Dubra A, Brainard BH, Morgan JIW. 
% "Non-invasive assessment of human cone photoreceptor function." 
% Biomed Opt Express 8(11): 5098-5112 and are
% encompassed in Figures 6/7, Equation 4.

clear;
close all;


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
        datavals(w) = {valdata};        
        
    end
%     hold off;
%     title( num2str(IDs(i)) ); xlabel('Log candelas'); ylabel('Amplitude response');
    
%     saveas(gcf, fullfile(pathname, [num2str(IDs(i)) 'lum_response.png']));
    
    rel_shifts(i,:) = action_spect_fit(wavelengths, irradiances, datavals, ciefunc, fit_type);
    legendIDs{i} = num2str(IDs(i));
    
    title( legendIDs{i} );%For getting the irradiance response as a function of subject
    axis([0 5 0 2.5]);
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
