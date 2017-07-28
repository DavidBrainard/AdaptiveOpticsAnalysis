
clear;
close all force;


controlBaseDir = uigetdir(pwd); %fullfile(pwd,'control'); 
stimBaseDir = uigetdir(pwd); %fullfile(pwd,'stim');
%%
controlDataNames = read_folder_contents(controlBaseDir ,'mat');
stimDataNames = read_folder_contents(stimBaseDir ,'mat');

wavelet = 'gaus3';

prefilt = designfilt('lowpassiir', 'FilterOrder', 5, 'PassbandFrequency', .06, 'PassbandRipple', .2);
% prefilt = designfilt('bandpassiir', 'FilterOrder', 6, 'PassbandFrequency1', .02, 'PassbandFrequency2', .1, 'PassbandRipple', .2);


clipind = 11; % Need to clip ends due to artifacts after filtering
stimind = 66;
maxind = 241;

allsigs=[];
affinitycurve = normpdf(0:maxind-stimind-1,16,33);
affinitycurve = affinitycurve./max(affinitycurve);

cutoff=0;
for j=1:length(controlDataNames)

    load(fullfile(controlBaseDir, controlDataNames{j}));
    
    % Legacy support
    if ~exist('norm_cell_reflectance','var')
        norm_cell_reflectance = norm_control_cell_reflectance;
        cell_times = control_cell_times;
    end
    
    % Remove the empty cells
    norm_control_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance)  );
    control_cell_times            = cell_times( ~cellfun(@isempty,cell_times) );
    
    naners = ~cellfun(@any, cellfun(@isnan, norm_control_cell_reflectance, 'UniformOutput',false));
    norm_control_cell_reflectance = norm_control_cell_reflectance( naners );
    control_cell_times            = control_cell_times( naners );
    
    controlcoeffs=[];
    controllabels=[];
    
    allcontrollabels{j} = [];
    allcontrolcoeffs{j} = [];
    
    for i=1:length(norm_control_cell_reflectance)
        
        times = control_cell_times{i};
        signal = norm_control_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);        

        interptimes = 0:maxind;
        interpsignal = interp1(times,signal,interptimes,'pchip');
        
        interpsignal = interpsignal(clipind:maxind);

        filtinterpsignal = filtfilt(prefilt, interpsignal);
        
        % Remove the remaining signal before the stimulus delivery
        interpsignal = interpsignal((stimind-clipind)+1:end);
        filtinterpsignal = filtinterpsignal((stimind-clipind)+1:end-1);
        
        allsigs=[allsigs; filtinterpsignal];
                
        during = filtinterpsignal( 1:33 );
        after = filtinterpsignal( 34:end-clipind );

        % Determine distance from max response to stim?
        [maxrespval,maxrespind]=max( abs( diff([during after])) );
        
        stim_affinity = affinitycurve(maxrespind)*maxrespval;
        
        derivduring = diff(during);
        derivafter = diff(after);
        
        [peakvals, peaks] = findpeaks([derivduring derivafter]);
        [troughvals, troughs] = findpeaks(-[derivduring derivafter]);
        
        if peaks(1) < troughs(1)
           listlen = length(peaks);
        else
           listlen = length(troughs);
        end
        
        for p=1:listlen
           if p > length(peaks) || p > length(troughs)
              break; 
           end
               
            pk_pk(p) = peakvals(p)+troughvals(p);
        end
        
        SWC = swt(interpsignal,4,'db4');
        
        % Features from Subasi et al 2007
        meanabscoeff = mean(abs(SWC'));
        meanpowercoeff = sum(SWC'.^2);
        stddevcoeff = std(SWC');
        coeffratio = [meanabscoeff(1)/meanabscoeff(2) meanabscoeff(2)/meanabscoeff(3) ...
                      meanabscoeff(3)/meanabscoeff(4) meanabscoeff(4)/meanabscoeff(5)];
        
        % Put together the feature lists
        controlcoeffs = [controlcoeffs; stim_affinity std(pk_pk) max(derivduring)-min(derivduring) meanabscoeff(4:5) meanpowercoeff(4:5) stddevcoeff(4:5) coeffratio(2:4)];        
%         controlcoeffs = [controlcoeffs; meanabscoeff(4:5) meanpowercoeff(4) stddevcoeff(3) coeffratio(3)];
        controllabels = [controllabels; {'control'}];
        
        
%         figure(1); title('Control cones'); hold on;
% %         plot(D5);
% %         hold on; 
%         plot( diff([during after]) );
%         axis([0 249 -1 1]);
% %         axis([0 250 -10 10]);
%         hold off;
        
        controlreconst(i,:) = filtinterpsignal;


    end
    clear norm_cell_reflectance;
    allcontrollabels{j} = [allcontrollabels{j}; controllabels];
    allcontrolcoeffs{j} = [allcontrolcoeffs{j}; controlcoeffs];
end




stimless=[];
stimd=[];
for j=1:length(stimDataNames)
    
    load(fullfile(stimBaseDir, stimDataNames{j}));
    
    % Legacy support
    if ~exist('norm_cell_reflectance','var')
        norm_cell_reflectance = norm_stim_cell_reflectance;
        cell_times = stim_cell_times;
    end
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance) );
    stim_cell_times            = cell_times(  ~cellfun(@isempty,cell_times) );
    
    naners = ~cellfun(@any, cellfun(@isnan, norm_stim_cell_reflectance, 'UniformOutput',false));
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( naners );
    stim_cell_times            = stim_cell_times( naners );
    
    stimdlabels=[];
    stimdcoeffs=[];
    
    allstimlabels{j} = [];
    allstimcoeffs{j} = []; 
    
    for i=1:length(norm_stim_cell_reflectance)
        times  = stim_cell_times{i};
        signal = norm_stim_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);
        
        interptimes = 0:maxind;
        interpsignal = interp1(times,signal,interptimes,'pchip');

        interpsignal = interpsignal(clipind:maxind);

        filtinterpsignal = filtfilt(prefilt, interpsignal);
        
        % Remove the remaining signal before the stimulus delivery
        interpsignal = interpsignal((stimind-clipind)+1:end);
        filtinterpsignal = filtinterpsignal((stimind-clipind)+1:end-1);
        
        allsigs=[allsigs; filtinterpsignal];
                
        during = filtinterpsignal( 1:33 );
        after  = filtinterpsignal( 34:end-clipind );

        % Determine distance from max response to stim?
        [maxrespval,maxrespind]=max( abs( diff([during after])) );
        
        stim_affinity = affinitycurve(maxrespind)*maxrespval;
        
        derivduring = diff(during);
        derivafter = diff(after);
        
        [peakvals, peaks] = findpeaks([derivduring derivafter]);
        [troughvals, troughs] = findpeaks(-[derivduring derivafter]);
        
        if peaks(1) < troughs(1)
           listlen = length(peaks);
        else
           listlen = length(troughs);
        end
        
        for p=1:listlen
           if p > length(peaks) || p > length(troughs)
              break; 
           end
               
            pk_pk(p) = peakvals(p)+troughvals(p);
        end
        
        SWC = swt(interpsignal,4,'db4');
        
        % Features from Subasi et al 2007
        meanabscoeff = mean(abs(SWC'));
        meanpowercoeff = sum(SWC'.^2);
        stddevcoeff = std(SWC');
        coeffratio = [meanabscoeff(1)/meanabscoeff(2) meanabscoeff(2)/meanabscoeff(3) ...
                      meanabscoeff(3)/meanabscoeff(4) meanabscoeff(4)/meanabscoeff(5)];
        
        stimd = [stimd i];

        if( stim_affinity > .05)
            stimdcoeffs = [stimdcoeffs; stim_affinity std(pk_pk) max(derivduring)-min(derivduring) meanabscoeff(4:5) meanpowercoeff(4:5) stddevcoeff(4:5) coeffratio(2:4)];
%             stimdcoeffs = [stimdcoeffs; meanabscoeff(4:5) meanpowercoeff(4) stddevcoeff(3) coeffratio(3)];
            stimdlabels = [stimdlabels; {'stimulus'}];
        end

%         figure(3); title(['Stim cones']); 
%         plot(interpsignal);hold on;
%         plot(filtinterpsignal); 
% % %         hold on; 
%         plot( diff([during after]) );
% %         axis([0 length([during after]) -1 1]);
% %         %plot(interpsignal); 
% % %         axis([0 250 -10 10]);
% %         axis([0 250 -1 1]);
%         hold off;
        
        
    end

    clear norm_cell_reflectance;
    allstimlabels{j} = [allstimlabels{j}; stimdlabels];
    allstimcoeffs{j} = [allstimcoeffs{j}; stimdcoeffs];
end



% Pick a random set from each type to fit from
numtrainingsets = 2;
controlSetInds = randperm( length(allcontrollabels), min( numtrainingsets, length(allcontrollabels)) )
stimSetInds = randperm( length(allstimlabels), min( numtrainingsets, length(allstimlabels)) )



%% Train our models.
% Aggregate all of the other data
trainingcoeff=[];
traininglabels=[];
validationcoeff=[];
validationlabels=[];
for o=1: max(length(allcontrollabels), length(allstimlabels)) 
    if all(o ~= controlSetInds) && all(o ~= stimSetInds) % If it's part of our training set, don't add it to the validation set.
        if o<=length(allcontrolcoeffs)
            validationcoeff  = [validationcoeff; allcontrolcoeffs{o}];    
            validationlabels = [validationlabels; allcontrollabels{o}];
        end
        if o<=length(allstimcoeffs)
            validationcoeff  = [validationcoeff; allstimcoeffs{o}];    
            validationlabels = [validationlabels; allstimlabels{o}];
        end
    end
end

for o=1:numtrainingsets
    
    trainingcoeff  = [trainingcoeff; allcontrolcoeffs{controlSetInds(o)}; allstimcoeffs{stimSetInds(o)}];    
    traininglabels = [traininglabels; allcontrollabels{controlSetInds(o)}; allstimlabels{stimSetInds(o)}];

end

% SVM
[pcacoeffs, pcascore, latent, ~, explained] = pca( trainingcoeff, 'VariableWeights','variance', 'Centered', true );
orthocoeffs = diag(std( trainingcoeff )) \ pcacoeffs; % Make the coefficients orthonormal
mu = mean(trainingcoeff);
stddev = std(trainingcoeff);

SVMModel = fitcsvm(pcascore(:,1:10),traininglabels,'KernelFunction','linear',...%'PolynomialOrder',2
                                             'KernelScale','auto',...
                                             'Standardize',true,...
                                             'BoxConstraint',10,'CrossVal','on','KFold',10,'OutlierFraction',0.05);

kfoldPercentModelLoss = 100*kfoldLoss(SVMModel,'mode','individual')

[minLoss, minInd]= min(kfoldPercentModelLoss)

minLossSVM = SVMModel.Trained{minInd};

centeredvalidationcoeff = (validationcoeff-repmat(mu,size(validationcoeff,1),1))./ repmat(stddev,size(validationcoeff,1),1);
validationscore = centeredvalidationcoeff*orthocoeffs;

% Estimate the classification error.
classificationPercentloss = 100*loss(minLossSVM,validationscore(:,1:10),validationlabels)


[predictions, score] = minLossSVM.predict(validationscore(:,1:10));

[x,y,t,auc,optrocpt]=perfcurve(predictions,max(score,[],2),'control');
optrocpt

figure(10); plot(x,y); hold on; %title(['SVM AUC: ' num2str(auc)]); 

confusionmat(validationlabels,predictions)

% Random forest
randforest = TreeBagger(200, trainingcoeff, traininglabels, 'ClassNames',{'control','stimulus'},...
                        'OOBPrediction','on','OOBPredictorImportance','on'); %,'SampleWithReplacement','off','InBagFraction',0.1);

figure(2); plot(oobError(randforest))

100*(1-error(randforest,validationcoeff,validationlabels));

[predictions, score] = randforest.predict(validationcoeff);

[x,y,t,auc,optrocpt]=perfcurve(predictions,max(score,[],2),'control');
optrocpt

figure(10); plot(x,y); xlabel('False Positive rate'); ylabel('True Positive rate'); hold off; %title(['Random forest AUC: ' num2str(auc)]);

confusionmat(validationlabels,predictions)