
clear;
close all force;


controlBaseDir = fullfile(pwd,'control'); %uigetdir(pwd); %
stimBaseDir = fullfile(pwd,'2p0'); %uigetdir(pwd); %
%%
controlDataNames = read_folder_contents(controlBaseDir ,'mat');
stimDataNames = read_folder_contents(stimBaseDir ,'mat');

allcontrollabels = cell(length(controlDataNames),1);
allcontrolcoeffs = cell(length(controlDataNames),1);
allstimlabels = cell(length(stimDataNames),1);
allstimcoeffs = cell(length(stimDataNames),1);

for j=1:length(controlDataNames)

    load(fullfile(controlBaseDir, controlDataNames{j}));
    
    % Legacy support
    if ~exist('norm_cell_reflectance','var')
        norm_cell_reflectance = norm_control_cell_reflectance;
        cell_times = control_cell_times;
    end
    %3
    [controlcoeffs, controllabels]=extract_features(cell_times, norm_cell_reflectance,[66 99],'control');
    
    clear norm_cell_reflectance;
    allcontrollabels{j} = [allcontrollabels{j}; controllabels];
    allcontrolcoeffs{j} = [allcontrolcoeffs{j}; controlcoeffs];
end


for j=1:length(stimDataNames)
    
    load(fullfile(stimBaseDir, stimDataNames{j}));
    
    % Legacy support
    if ~exist('norm_cell_reflectance','var')
        
        norm_cell_reflectance = norm_stim_cell_reflectance;
        cell_times = stim_cell_times;
    end
    
     [stimdcoeffs, stimdlabels]=extract_features(cell_times, norm_cell_reflectance,[66 99],'stimulus');    

    clear norm_cell_reflectance;
    allstimlabels{j} = [allstimlabels{j}; stimdlabels];
    allstimcoeffs{j} = [allstimcoeffs{j}; stimdcoeffs];
end





%% Train our models.
% Pick a random set from each type to fit from
numtrainingsets = 1;
% controlSetInds = randperm( length(allcontrollabels), min( numtrainingsets, length(allcontrollabels)) )
% stimSetInds = randperm( length(allstimlabels), min( numtrainingsets, length(allstimlabels)) )

controlSetInds=3;
stimSetInds=3;

trainingControlFiles=controlDataNames(controlSetInds)
trainingStimFiles=stimDataNames(stimSetInds)



for j=1:size(allstimcoeffs{1},2)

    figure(j); hold off; histogram(allcontrolcoeffs{controlSetInds}(:,j),100); hold on; histogram(allstimcoeffs{stimSetInds}(:,j),100);

end


% Aggregate all of the other data
trainingcoeff=[];
traininglabels=[];
validationcoeff=[];
validationlabels=[];
for o=1: max(length(allcontrollabels), length(allstimlabels)) 
    if all(o ~= controlSetInds)
        if o<=length(allcontrolcoeffs)
            validationcoeff  = [validationcoeff; allcontrolcoeffs{o}];    
            validationlabels = [validationlabels; allcontrollabels{o}];
        end
    end
    
    if all(o ~= stimSetInds) % If it's part of our training set, don't add it to the validation set.
        if o<=length(allstimcoeffs)
            validationcoeff  = [validationcoeff; allstimcoeffs{o}];    
            validationlabels = [validationlabels; allstimlabels{o}];
        end
    end
end

stimtrainingcoeff=[];
stimtraininglabels=[];
conttrainingcoeff=[];
conttraininglabels=[];
for o=1:numtrainingsets
    
    stimtrainingcoeff  = [stimtrainingcoeff; allstimcoeffs{stimSetInds(o)}];    
    stimtraininglabels = [stimtraininglabels; allstimlabels{stimSetInds(o)}];
    conttrainingcoeff  = [conttrainingcoeff; allcontrolcoeffs{controlSetInds(o)}];    
    conttraininglabels = [conttraininglabels; allcontrollabels{controlSetInds(o)}];
end

% Trim down the training sets so that they're the same size.
if size(stimtrainingcoeff,1) > size(conttrainingcoeff,1)
    trainingcoeff = [stimtrainingcoeff(1:size(conttrainingcoeff,1),:); conttrainingcoeff];
    traininglabels = [stimtraininglabels(1:size(conttraininglabels,1),:); conttraininglabels];
else
    trainingcoeff = [conttrainingcoeff(1:size(stimtrainingcoeff,1),:); stimtrainingcoeff];
    traininglabels = [conttraininglabels(1:size(stimtraininglabels,1),:); stimtraininglabels]; 
end


varNames = {'Affinity','AUC','Max_Deriv_Ampl', 'Pk_pk_std_dev', ...
            'Mean_Abs_1D','Mean_Abs_2D','Mean_Abs_3D','Mean_Abs_4D','Mean_Abs_4A',...%'Mean_Abs_6D','Mean_Abs_7D','Mean_Abs_8D',...
            'Mean_Power_1D','Mean_Power_2D','Mean_Power_3D','Mean_Power_4D','Mean_Power_4A',...%'Mean_Power_5D','Mean_Power_6D','Mean_Power_7D','Mean_Power_8D'...
            'Std_Dev_1D','Std_Dev_2D','Std_Dev_3D','Std_Dev_4D','Std_Dev_4A',...%'Std_Dev_5D','Std_Dev_6D','Std_Dev_7D','Std_Dev_8D',...
            'Coeff_Ratio_1D2D','Coeff_Ratio_2D3D','Coeff_Ratio_3D4D','Coeff_Ratio_4D4A','Training_Labels'};

trainingTable = [array2table(trainingcoeff) cell2table(traininglabels)];
trainingTable.Properties.VariableNames = varNames;
validationTable = [array2table(validationcoeff), cell2table(validationlabels)];
validationTable.Properties.VariableNames = varNames;
% Train the classifers

%% SVM

SVMModel = fitcsvm(trainingcoeff,traininglabels,'KernelFunction','polynomial','PolynomialOrder',2,...
                                             'KernelScale','auto',...
                                             'Standardize',true,...
                                             'BoxConstraint',1,...
                                             'KFold',5,...
                                             'ClassNames', {'control'; 'stimulus'});

% Perform cross-validation
partmod = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partmod, 'LossFun', 'ClassifError');

% Compute validation predictions and scores
[validationPredictions, validationScores] = kfoldPredict(partmod)

[predictions, score] = SVMModel.predict(validationcoeff);

[x,y,t,auc,optrocpt]=perfcurve(predictions,max(score,[],2),'stimulus');
optrocpt

figure(10); plot(x,y); hold on; %title(['SVM AUC: ' num2str(auc)]); 

confusionmat(validationlabels,predictions)

svm_classifier = minLossSVM;
% svm_classifier = [];

% Random forest
randforest = TreeBagger(200, trainingcoeff, traininglabels, 'ClassNames',{'stimulus','control'},...
                        'OOBPrediction','on','OOBPredictorImportance','on'); %,'SampleWithReplacement','off','InBagFraction',0.1);

figure(2); plot(oobError(randforest))

100*(1-error(randforest,validationcoeff,validationlabels))

[predictions, score] = randforest.predict(validationcoeff);

[x,y,t,auc,optrocpt]=perfcurve(predictions,max(score,[],2),'stimulus');
optrocpt
% auc

figure(10); plot(x,y); xlabel('False Positive rate'); ylabel('True Positive rate'); hold off; %title(['Random forest AUC: ' num2str(auc)]);

confusionmat(validationlabels,predictions)


randforest_classifier = trainedClassifier;
save(['trained_classifier_' date '.mat'],'randforest_classifier','svm_classifier','trainingControlFiles','trainingStimFiles');