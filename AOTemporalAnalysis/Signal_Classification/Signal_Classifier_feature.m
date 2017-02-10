
clear;
close all force;


controlBaseDir = fullfile(pwd,'control'); %uigetdir(pwd);
stimBaseDir = fullfile(pwd,'stim');

controlDataNames = read_folder_contents(controlBaseDir ,'mat');
stimDataNames = read_folder_contents(stimBaseDir ,'mat');

wavelet = 'gaus3';

prefilt = designfilt('lowpassiir', 'FilterOrder', 5, 'PassbandFrequency', .06, 'PassbandRipple', .2);
% prefilt = designfilt('bandpassiir', 'FilterOrder', 6, 'PassbandFrequency1', .02, 'PassbandFrequency2', .1, 'PassbandRipple', .2);

stimd=0;
stimless=0;
x45=[];
x45ind=[];
x92=[];
x92ind =[];
x04=[];
x04ind =[];

controlcoeffs=[];
controllabels=[];
stimdcoeffs=[];
stimdlabels=[];

allsigs=[];

clipind = 11; % Need to clip ends due to artifacts after filtering
stimind = 66;
maxind = 240;


affinitycurve = normpdf(0:maxind-stimind-1,16,33);
affinitycurve = affinitycurve./max(affinitycurve);

cutoff=0;
for j=1:1%length(controlDataNames)

    load(fullfile(controlBaseDir, controlDataNames{j}));
    
    % Remove the empty cells
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
    control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );
    
    for i=1:length(norm_control_cell_reflectance)
        
        times = control_cell_times{i};
        signal = norm_control_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);        

        interptimes = 0:maxind;
        interpsignal = interp1(times,signal,interptimes,'pchip');
        
        interpsignal = interpsignal(clipind:maxind);

        D5 = filtfilt(prefilt, interpsignal);
        
        % Remove the remaining signal before the stimulus delivery
        D5 = D5((stimind-clipind)+1:end-1);
        
        allsigs=[allsigs; D5];
                
        during = D5( 1:33 );
        after = D5( 34:end-clipind );

        % Determine distance from max response to stim?
        [maxrespval,maxrespind]=max( abs( diff([during after])) );
        
        stim_affinity = affinitycurve(maxrespind)*maxrespval;
        
        derivduring = diff(during);
        derivafter = diff(after);
        
        % Put together the feature lists
        controlcoeffs = [controlcoeffs; stim_affinity max(derivduring)-min(derivduring) ];
        
        controllabels = [controllabels; {'control'}];
        
        
        figure(1); title('Control cones'); hold on;
%         plot(D5);
%         hold on; 
        plot( diff([during after]) );
        axis([0 249 -1 1]);
%         axis([0 250 -10 10]);
        hold off;
        
        
%         N = length(D5);
%         xdft = fft(D5);
%         xdft = xdft(1:N/2+1);
%         psdx = (1/(2*pi*N)) * abs(xdft).^2;
%         psdx(2:end-1) = 2*psdx(2:end-1);
%         freq = 0:(2*pi)/N:pi;

%         figure(2);
%         plot(freq/pi,10*log10(psdx))
%         grid on
%         title('Periodogram Using FFT')
%         xlabel('Normalized Frequency (\times\pi rad/sample)')
%         ylabel('Power/Frequency (dB/rad/sample)')

        controlreconst(i,:) = D5;


    end
end

stimless=[];
stimd=[];
for j=1:1%length(stimDataNames)
    
    load(fullfile(stimBaseDir, stimDataNames{j}));
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );
    
    
    for i=1:length(norm_stim_cell_reflectance)
        times  = stim_cell_times{i};
        signal = norm_stim_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);
        
        interptimes = 0:maxind;
        interpsignal = interp1(times,signal,interptimes,'pchip');

        interpsignal = interpsignal(clipind:maxind);

        D5 = filtfilt(prefilt, interpsignal);
        
        % Remove the remaining signal before the stimulus delivery
        D5 = D5((stimind-clipind)+1:end-1);
        
        allsigs=[allsigs; D5];
                
        during = D5( 1:33 );
        after = D5( 34:end-clipind );

        % Determine distance from max response to stim?
        [maxrespval,maxrespind]=max( abs( diff([during after])) );
        
        stim_affinity = affinitycurve(maxrespind)*maxrespval;
        
        derivduring = diff(during);
        derivafter = diff(after);
        
        stimd = [stimd i];

        if( stim_affinity > .05)
            stimdcoeffs = [stimdcoeffs; stim_affinity max(derivduring)-min(derivduring)  ];
            stimdlabels = [stimdlabels; {'stimulus'}];
        end

        figure(3); title(['Stim cones']); hold on; 
% %         plot(D5); 
% %         hold on; 
        plot( diff([during after]) );
%         axis([0 length([during after]) -1 1]);
%         %plot(interpsignal); 
% %         axis([0 250 -10 10]);
        axis([0 250 -1 1]);
        hold off;
        
        
%         N = length(D5);
%         xdft = fft(D5);
%         xdft = xdft(1:N/2+1);
%         psdx = (1/(2*pi*N)) * abs(xdft).^2;
%         psdx(2:end-1) = 2*psdx(2:end-1);
%         freq = 0:(2*pi)/N:pi;
% 
%         figure(4);
%         plot(freq/pi,10*log10(psdx))
%         grid on
%         title('Periodogram Using FFT')
%         xlabel('Normalized Frequency (\times\pi rad/sample)')
%         ylabel('Power/Frequency (dB/rad/sample)')
        
    end
end

alllabels = [controllabels; stimdlabels];
allcoeffs = [controlcoeffs; stimdcoeffs];

% allcoeffs = allcoeffs( [1:3769 3771:end],: );
% alllabels = alllabels( [1:3769 3771:end],: );

[pcacoeffs, pcascore, ~, ~, explained] = pca( allcoeffs );
explained


projected = allcoeffs*pcacoeffs;

figure;
gscatter(projected(:,1), projected(:,2), alllabels);


SVMModel = fitcsvm(allcoeffs,alllabels,'KernelFunction','linear',...
                                       'KernelScale','auto',...
                                       'Standardize',true,...
                                       'BoxConstraint',10,'OutlierFraction',0.1);
sv = SVMModel.SupportVectors;

d = 0.02; % Step size of the grid
[x1Grid,x2Grid] = meshgrid( min(allcoeffs(:,1)):d:max(allcoeffs(:,1)),...
                            min(allcoeffs(:,2)):d:max(allcoeffs(:,2)) );
xGrid = [x1Grid(:),x2Grid(:)];        % The grid
[~,scores1] = predict(SVMModel,xGrid); % The scores

figure;
h(1:2) = gscatter(allcoeffs(:,1), allcoeffs(:,2), alllabels);
hold on
h(3) = plot(allcoeffs(SVMModel.IsSupportVector,1),...
            allcoeffs(SVMModel.IsSupportVector,2),'ko','MarkerSize',10);
    % Support vectors
contour(x1Grid, x2Grid, reshape(scores1(:,2), size(x1Grid)),[0 0],'k');
    % Decision boundary
title('Scatter Diagram with the Decision Boundary')
legend({SVMModel.ClassNames{1},SVMModel.ClassNames{2},'Support Vectors'},'Location','Best');
hold off


% d = 0.02; % Step size of the grid
% [x1Grid,x2Grid,x3Grid] = meshgrid( min(allcoeffs(:,1)):d:max(allcoeffs(:,1)),...
%                                    min(allcoeffs(:,2)):d:max(allcoeffs(:,2)),...
%                                    min(allcoeffs(:,3)):d:max(allcoeffs(:,3)) );
% xGrid = [x1Grid(:),x2Grid(:),x3Grid(:)]; % The grid
% [~,scores1] = predict(SVMModel,xGrid); % The scores
% 
% figure;
% h(1:2) = scatter3(allcoeffs(:,1), allcoeffs(:,2), allcoeffs(:,3));
% hold on
% h(3) = plot3(allcoeffs(SVMModel.IsSupportVector,1),...
%              allcoeffs(SVMModel.IsSupportVector,2),...
%              allcoeffs(SVMModel.IsSupportVector,3),'ko','MarkerSize',10);
