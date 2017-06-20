% SizePsychometric
%
% Little script to predict how size psychometric functions would look,
% given detection psychometric functions.
%
% Assume:
%   1) All spots are within a summation area, so that detection
%   psychometric functions based on energy all overlay.
%   2) A signal-detection y/n psychometric model, say yes if response on
%   each trial exceeds a criterion.
%   3) Additive Gaussian noise with signal independent variance.  Fix
%   variance at 1.
%
% Requires
%   BrainardLabToolbox. Run:
%   tbUse('BrainardLabToolbox');

% 6/20/17  dhb  Wrote it.

%% Clear
clear; close all;

%% Parameters

% Noise sd. This, in essence, defines the units of
% intensity for the simulation and we might as well set
% this to 1.
noiseSd = 1;

% Stimulus diameters are in are arbitrary units
testStimulusDiameters = [1 2 4 8 16];

% For the size judgments, there is a reference stimulus.
% Give its intensity and diameter here.
referenceStimulusDiameter = 4;
referenceStimulusIntensity = 10.^-0.585;

% This determines the false alarm rate.  We express as a
% factor of the noiseSd.
observerCriterion = 2.5*noiseSd;

% Intensities for psychometric functions
lowLogIntensity = -3;
highLogIntensity = 1;
psychometricIntensities = logspace(lowLogIntensity,highLogIntensity,100);

%% For convenience
testStimulusAreas = pi*(testStimulusDiameters/2).^2;
referenceStimulusArea = pi*(referenceStimulusDiameter/2).^2;

% Reference intensity probability detect
pRefYes = 1-normcdf(observerCriterion,referenceStimulusIntensity*referenceStimulusArea,noiseSd);
fprintf('Probability of detecting reference stimulus: %0.2f\n',pRefYes);

%% Simulate out y/n detection psychometric functions
%
% pOneYes is for detecting one test stimulus presented on a trial
% pTwoYes is for detecting both reference and test
ynFigure = figure; clf;
set(ynFigure,'Position',[100 100 400 850]);
for ii = 1:length(testStimulusDiameters)
    pOneYes(ii,:) = 1-normcdf(observerCriterion,psychometricIntensities*testStimulusAreas(ii),noiseSd);
    pTwoYes(ii,:) = pOneYes(ii,:).*pRefYes;
    subplot(length(testStimulusDiameters),1,ii); hold on;
    plot(log10(psychometricIntensities),pOneYes(ii,:),'-','LineWidth',4,'Color',[0 0.5 0.75]);
    plot(log10(psychometricIntensities),pTwoYes(ii,:),'-','LineWidth',4,'Color','r');
    xlim([lowLogIntensity highLogIntensity]);
    ylim([0 1]);
    xlabel('Stimulus Intensity');
    ylabel('Probability See');
    title(sprintf('Detection, Ref Size %d, Test Size %d',referenceStimulusDiameter,testStimulusDiameters(ii)));
    if (testStimulusDiameters(ii) < referenceStimulusDiameter)
        legend({'Detect Test', 'Detect Both'},'Location','NorthWest');
    else
        legend({'Detect Test', 'Detect Both'},'Location','SouthEast');
    end
end

%% Simulate out size judgment psychometric functions
%
% Whichever stimulus produces more energy on a trial is judged larger.
% But, only trials where both stimuli are seen are included.
%
% I'm sure there is an analytic formula, but simulation is so easy.
nSimulate = 10000;
sizeFigure = figure; clf;
set(sizeFigure,'Position',[100 200 400 850]);
for ii = 1:length(testStimulusDiameters)
    oneYes = zeros(nSimulate,length(psychometricIntensities));
    refYes = zeros(nSimulate,length(psychometricIntensities));
    oneLarger = zeros(nSimulate,length(psychometricIntensities));
    for jj = 1:length(psychometricIntensities)
        oneVals = normrnd(psychometricIntensities(jj)*testStimulusAreas(ii),noiseSd,1,nSimulate);
        refVals = normrnd(referenceStimulusIntensity*referenceStimulusArea,noiseSd,1,nSimulate);
        pOneYesSimulate(ii,jj) = sum(oneVals > observerCriterion)/length(oneVals);
        pRefYesSimulate(ii,jj) = sum(refVals > observerCriterion)/length(refVals);
        pTwoYesSimulate(ii,jj) = sum(oneVals > observerCriterion & refVals > observerCriterion)/length(oneVals);
        
        % Find trials where both are seen, and where both are seen and test
        % produces larger response.
        nSeeBoth = sum((oneVals > observerCriterion) & (refVals > observerCriterion));
        nTestBigger = sum((oneVals > observerCriterion) & (refVals > observerCriterion) & (oneVals > refVals));
        
        % Only estimate probabiity of seing test as bigger if we have
        % enough trials where both are seen.
        if (nSeeBoth > nSimulate/10)
            pSeeBoth(ii,jj) = nTestBigger/nSeeBoth;
        else
            pSeeBoth(ii,jj) = NaN;
        end
    end
    
    % Plot simulated psychometric on top of analytic psychometric
    % This is mainly to check the simulation logic
    figure(ynFigure);
    subplot(length(testStimulusDiameters),1,ii); hold on;
    plot(log10(psychometricIntensities),pOneYesSimulate(ii,:),'k:','LineWidth',2);
    plot(log10(psychometricIntensities),pTwoYesSimulate(ii,:),'k:','LineWidth',2);
    
    % Plot the size judgment psychometric
    figure(sizeFigure);
    subplot(length(testStimulusDiameters),1,ii); hold on;
    plot(log10(psychometricIntensities),pSeeBoth(ii,:),'-','LineWidth',4,'Color','g');
    xlim([lowLogIntensity highLogIntensity]);
    ylim([0 1]);
    xlabel('Stimulus Intensity');
    ylabel('Probability Test Bigger');
    title(sprintf('Size Judgment, Ref Size %d, Test Size %d',referenceStimulusDiameter,testStimulusDiameters(ii)));
end

%% Save the figures
FigureSave('ynPsychometric',ynFigure,'pdf');
FigureSave('sizePsychometric',sizeFigure,'pdf');

