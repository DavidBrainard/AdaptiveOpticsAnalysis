
clear;
close all force;


controlBaseDir = fullfile(pwd,'control'); %uigetdir(pwd); %
stimBaseDir = fullfile(pwd,'2p0'); %uigetdir(pwd); %
%%
controlDataNames = read_folder_contents(controlBaseDir ,'mat');
stimDataNames = read_folder_contents(stimBaseDir ,'mat');
%%
load(fullfile(controlBaseDir, controlDataNames{1}));

allcontrollabels = repmat({'control'},[length(norm_cell_reflectance),1]);
allcontrolcoeffs = nan(length(norm_cell_reflectance), length(controlDataNames)*9);
allcontrolcoords = coords_used;

for j=1:length(controlDataNames)
    controlDataNames{j}
    load(fullfile(controlBaseDir, controlDataNames{j}));
    
    allcontrolcoords = union(allcontrolcoords, coords_used,'rows');
    
    % These all must be the same length! (Same coordinate set)
    if size(norm_cell_reflectance,1) == size(allcontrolcoeffs,1)

        [coeffs, labels, varNames]=extract_features_v2(cell_times, norm_cell_reflectance,[66 132],'control');

        indme = ((j-1)*9)+1;
        allcontrolcoeffs(:,indme:indme+8) = coeffs;
    else
        error('Dataset length doesn''t match!');
    end
end

tokeep  = all(~isnan(allcontrolcoeffs),2);

allcontrolcoeffs = allcontrolcoeffs(tokeep,:);
allcontrollabels = allcontrollabels(tokeep);
allcontrolcoords = allcontrolcoords(tokeep,:);

%%
clear coords_used

load(fullfile(stimBaseDir, stimDataNames{1}));
allstimlabels = repmat({'stimulus'},[length(norm_cell_reflectance),1]);
allstimcoeffs = nan(length(norm_cell_reflectance), length(stimDataNames)*9);
allstimcoords = coords_used;

for j=1:length(stimDataNames)
    stimDataNames{j}
    load(fullfile(stimBaseDir, stimDataNames{j}));
    
    allstimcoords = union(allstimcoords, coords_used,'rows');
    
    % These all must be the same length! (Same coordinate set)
    if size(norm_cell_reflectance,1) == size(allstimcoeffs,1) && size(allstimcoords,1) == size(allstimcoeffs,1)

        [coeffs, labels, varNames]=extract_features_v2(cell_times, norm_cell_reflectance,[66 132],'stimulus');
        
        indme = ((j-1)*9)+1;
        allstimcoeffs(:,indme:indme+8) = coeffs;
    else
        error('Dataset length doesn''t match!');
    end
end

tokeep  = all(~isnan(allstimcoeffs),2);

allstimcoeffs = allstimcoeffs(tokeep,:);
allstimlabels = allstimlabels(tokeep);
allstimcoords = allstimcoords(tokeep,:);

%% Train our models.
% Pick a random set from each type to fit from
numtrainingsets = 2;
controlSetInds = randperm( length(allcontrollabels), min( numtrainingsets, length(allcontrollabels)) )
stimSetInds = randperm( length(allstimlabels), min( numtrainingsets, length(allstimlabels)) )

% controlSetInds=[1 2 4];
% stimSetInds=[1 2 4];

trainingControlFiles=controlDataNames(controlSetInds)
trainingStimFiles=stimDataNames(stimSetInds)



% for j=1:size(allstimcoeffs{1},2)
% 
%     figure(j); hold off; histogram(allcontrolcoeffs{controlSetInds}(:,j),100); hold on; histogram(allstimcoeffs{stimSetInds}(:,j),100);
% 
% end


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

trainingTable = [array2table(trainingcoeff) cell2table(traininglabels)];
trainingTable.Properties.VariableNames = varNames;
validationTable = [array2table(validationcoeff), cell2table(validationlabels)];
validationTable.Properties.VariableNames = varNames;
% Train the classifers

svm_classifier= trainedClassifier;
save(['trained_classifier_' date '.mat'],'svm_classifier','trainingControlFiles','trainingStimFiles');
%% SVM

SVMModel = fitcsvm(trainingcoeff,traininglabels,'KernelFunction','polynomial','PolynomialOrder',2,...
                                             'KernelScale','auto',...
                                             'Standardize',true,...
                                             'BoxConstraint',1,...                                             
                                             'ClassNames', {'control'; 'stimulus'});

% Perform cross-validation
partmod = crossval(SVMModel, 'KFold', 5);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partmod, 'LossFun', 'ClassifError');

% Compute validation predictions and scores
[validationPredictions, validationScores] = kfoldPredict(partmod)

[predictions, score] = SVMModel.predict(validationcoeff);

[x,y,t,auc,optrocpt]=perfcurve(predictions,max(score,[],2),'stimulus');
optrocpt

figure(10); plot(x,y); hold on; %title(['SVM AUC: ' num2str(auc)]); 

confusionmat(validationlabels,predictions)

svm_classifier = SVMModel;
% svm_classifier = [];

%% Random forest
randforest = TreeBagger(200, trainingcoeff, traininglabels, 'ClassNames',{'stimulus','control'},...
                        'OOBPrediction','on','OOBPredictorImportance','on'); %,'SampleWithReplacement','off','InBagFraction',0.1);

figure(2); plot(oobError(randforest))

100*(1-error(randforest,validationcoeff,validationlabels));

[predictions, score] = randforest.predict(validationcoeff);

[x,y,t,auc,optrocpt]=perfcurve(predictions,max(score,[],2),'stimulus');
optrocpt;
% auc

figure(10); plot(x,y); xlabel('False Positive rate'); ylabel('True Positive rate'); hold off; %title(['Random forest AUC: ' num2str(auc)]);

confusionmat(validationlabels,predictions)


randforest_classifier = randforest;
save(['trained_classifier_' date '.mat'],'randforest_classifier','svm_classifier','trainingControlFiles','trainingStimFiles');

%
