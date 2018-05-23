
% clear;
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
numpcacoeff=15;

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

        derivduraft = diff([during after]);
        
        if exist('pcacoeffs','var')
            transderv = derivduraft*pcacoeffs;

            controlcoeffs = [controlcoeffs; transderv(1:numpcacoeff)];
            controllabels = [controllabels; {'control'}];
        end
        
%         figure(1); title('Control cones'); hold on;
% %         plot(D5);
% %         hold on; 
%         plot( diff([during after]) );
%         axis([0 249 -1 1]);
% %         axis([0 250 -10 10]);
%         hold off;
%         
        controlsig(i,:) = derivduraft;

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

        derivduraft = diff([during after]);
                
        if exist('pcacoeffs','var')
            transderv = derivduraft*pcacoeffs;        
        

            stimdcoeffs = [stimdcoeffs; transderv(1:numpcacoeff)];
            stimdlabels = [stimdlabels; {'stimulus'}];
        end
%         figure(3); title(['Stim cones']); hold on; 
% % %         plot(D5); 
% % %         hold on; 
%         plot( diff([during after]) );
% %         axis([0 length([during after]) -1 1]);
% %         %plot(interpsignal); 
% % %         axis([0 250 -10 10]);
%         axis([0 250 -1 1]);
%         hold off;
        
        stimsig(i,:) = derivduraft;
    end
end


[pcacoeffs, pcascore, latent, ~, explained] = pca( [controlsig; stimsig] );
explained

alllabels = [controllabels; stimdlabels];
allcoeffs = [controlcoeffs; stimdcoeffs];

% dataSetInds = randperm(length(allcoeffs));
% dataSetInds(1)
SVMModel = fitcsvm(allcoeffs,alllabels,'KernelFunction','linear',...
                                             'KernelScale','auto',...
                                             'Standardize',true,...
                                             'BoxConstraint',10,...%'CrossVal','on','KFold',10,
                                             'OutlierFraction',0.1);

% kfoldPercentModelLoss = 100*kfoldLoss(SVMModel)

beta = SVMModel.Beta;
hyperplane = null(beta');
transformedTestData = normc([beta hyperplane])' * allcoeffs';
transformedTestData = transformedTestData';

plot(transformedTestData');
