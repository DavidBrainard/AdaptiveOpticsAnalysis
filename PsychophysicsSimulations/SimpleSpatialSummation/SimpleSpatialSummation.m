function SimpleSpatialSummation


%% Clear
clear; close all;

%% Parameters
%
% These are areas in terms the area of a single cone, not linear dimension
minSpotSize = 1;
nSpotSizes = 20;
maxSpotSize = 90;
nPositions = 100;
if (rem(nPositions,2) ~= 0)
    error('nPositions must be even');
end
detectorSizes = [2 5 10 20];

% Background retinal irradiance (isomerizations per cone per stimulus duration)
backgroundIsomerizations = 100;

% Criterion fraction correct
criterionFractionCorrect = 0.75;

%% Poisson calculation
%
% Loop over spot sizes and find threshold for each
weights = [];
spotSizes = unique(round(logspace(log10(minSpotSize),log10(maxSpotSize),nSpotSizes)));
for dd = 1:length(detectorSizes)
    detectorSize = detectorSizes(dd);
    
    for ss = 1:length(spotSizes)
        
        % Set spot size
        spotSize = spotSizes(ss);
        
        % Threshold guesses and bounds in irradiance expressed as
        % isomerizations per cone per stimulus duration
        if (ss == 1)
            initialThreshold = 1;
        else
            initialThreshold = threshold(dd,ss-1);
        end
        maxThreshold = 100000;
        minThreshold = 0;
        
        % Use fmincon to find threshold irradiance
        options = optimset('fmincon');
        options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
        threshold(dd,ss) = fmincon(@(x)FindThresholdFunction(x,backgroundIsomerizations,detectorSize,spotSize,criterionFractionCorrect,weights),initialThreshold,[],[],[],[],minThreshold,maxThreshold,[],options);
        thresholdEnergy(dd,ss) = spotSize*threshold(dd,ss);
        [~,fractionCorrect(dd,ss)] = FindThresholdFunction(threshold(dd,ss),backgroundIsomerizations,detectorSize,spotSize,criterionFractionCorrect,weights);
        
        fprintf('Detector size %d, spot size %d, threshold irradiance %0.2f, percent correct: %0.2f\n',detectorSize,spotSize,threshold(dd,ss),round(100*fractionCorrect(dd,ss)));
        
    end
end

% Plot
theColors = ['c' 'r' 'g' 'b' 'k' 'y'];
theFig = figure; clf; hold on;
set(gcf,'Position',[100 100 1050 560]);
subplot(1,2,1); hold on;
for dd = 1:length(detectorSizes)
    theColor = theColors(rem(dd,length(theColors)) + 1);
    plot(log10(spotSizes),log10(threshold(dd,:)),[theColor 'o'],'MarkerSize',10);
    plot(log10(spotSizes),log10(threshold(dd,:)),[theColor ],'LineWidth',2);
    plot([log10(detectorSizes(dd)) log10(detectorSizes(dd))],[0 0.25],[theColor ],'LineWidth',3);
end
xlim([0 3]);
ylim([0 3]);
axis('square');
xlabel('Log Stimulus Area');
ylabel('Log Threshold Irradiance');
subplot(1,2,2); hold on;
for dd = 1:length(detectorSizes)
    theColor = theColors(rem(dd,length(theColors)) + 1);
    plot(log10(spotSizes),log10(thresholdEnergy(dd,:)),[theColor 'o'],'MarkerSize',10);
    plot(log10(spotSizes),log10(thresholdEnergy(dd,:)),[theColor ],'LineWidth',2);
    plot([log10(detectorSizes(dd)) log10(detectorSizes(dd))],[0.5 0.75],[theColor ],'LineWidth',3);
    
end
xlim([0 3]);
ylim([0.5 3.5]);
axis('square');
xlabel('Log Stimulus Area');
ylabel('Log Threshold Energy');
FigureSave('SimplePoissonSpatialSummation',gcf,'pdf');


%% Gaussian calculation
%
% Loop over spot sizes and find threshold for each
spotSizes = unique(round(logspace(log10(minSpotSize),log10(maxSpotSize),nSpotSizes)));
for dd = 1:length(detectorSizes)
    detectorSize = detectorSizes(dd);
    positions = 1:nPositions;
    
    % Set Gaussian weights
    weights = normpdf(positions-nPositions/2-0.5,0,detectorSize);
    weights = detectorSize*weights/sum(weights);
    
    % This redoes the hard aperture using the Gaussian flow.
    DOSPOT = false;
    if (DOSPOT)
        [~,maxIndex] = sort(abs(weights),'descend');
        detectorPositionIndex = maxIndex(1:detectorSize);
        weights = zeros(size(weights));
        weights(detectorPositionIndex) = 1;
    end 
    % weights = zeros(size(positions));
    % firstPosition = floor((nPositions/2-detectorSize/2));
    % weights(firstPosition:firstPosition+detectorSize-1) = 1;
    
    for ss = 1:length(spotSizes)
        
        % Set spot size
        spotSize = spotSizes(ss);
        
        % Threshold guesses and bounds in irradiance expressed as
        % isomerizations per cone per stimulus duration
        if (ss == 1)
            initialThreshold = 1;
        else
            initialThreshold = thresholdGauss(dd,ss-1);
        end
        maxThreshold = 100000;
        minThreshold = 0;
        
        % Use fmincon to find threshold irradiance
        options = optimset('fmincon');
        options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
        thresholdGauss(dd,ss) = fmincon(@(x)FindThresholdFunction(x,backgroundIsomerizations,detectorSize,spotSize,criterionFractionCorrect,weights),initialThreshold,[],[],[],[],minThreshold,maxThreshold,[],options);
        thresholdEnergyGauss(dd,ss) = spotSize*thresholdGauss(dd,ss);
        [~,fractionCorrect(dd,ss)] = FindThresholdFunction(thresholdGauss(dd,ss),backgroundIsomerizations,detectorSize,spotSize,criterionFractionCorrect,weights);
        
        fprintf('Detector size %d, spot size %d, threshold irradiance %0.2f, percent correct: %0.2f\n',detectorSize,spotSize,thresholdGauss(dd,ss),round(100*fractionCorrect(dd,ss)));
        
    end
end

% Plot
theColors = ['c' 'r' 'g' 'b' 'k' 'y'];
theFig = figure; clf; hold on;
set(gcf,'Position',[100 100 1050 560]);
subplot(1,2,1); hold on;
for dd = 1:length(detectorSizes)
    theColor = theColors(rem(dd,length(theColors)) + 1);
    plot(log10(spotSizes),log10(thresholdGauss(dd,:)),[theColor 'o'],'MarkerSize',10);
    plot(log10(spotSizes),log10(thresholdGauss(dd,:)),[theColor ],'LineWidth',2);
    plot([log10(detectorSizes(dd)) log10(detectorSizes(dd))],[0 0.25],[theColor ],'LineWidth',3);
end
xlim([0 3]);
ylim([0 3]);
axis('square');
xlabel('Log Stimulus Area');
ylabel('Log Threshold Irradiance');
subplot(1,2,2); hold on;
for dd = 1:length(detectorSizes)
    theColor = theColors(rem(dd,length(theColors)) + 1);
    plot(log10(spotSizes),log10(thresholdEnergyGauss(dd,:)),[theColor 'o'],'MarkerSize',10);
    plot(log10(spotSizes),log10(thresholdEnergyGauss(dd,:)),[theColor ],'LineWidth',2);
    plot([log10(detectorSizes(dd)) log10(detectorSizes(dd))],[0.5 0.75],[theColor ],'LineWidth',3);
    
end
xlim([0 3]);
ylim([0.5 3.5]);
axis('square');
xlabel('Log Stimulus Area');
ylabel('Log Threshold Energy');
FigureSave('SimpleGaussianSpatialSummation',gcf,'pdf');

end


%% Error function.
%
% Comp
function [f,fractionCorrect] = FindThresholdFunction(x,backgroundIsomerizations,detectorSize,spotSize,criterionFractionCorrect,weights)

% If weights is empty, assume a hard summation error and do the Poisson
% calculation.  Otherwise weight up the cones and do a Gaussian
% calculations.
if (isempty(weights))
    % Build up background and test stimuli as vectors.  We don't model the
    % spatial geometry, just count cones.  The vector just specifies the cones
    % within the detector, becaues those are the only ones we pay attention to.
    % The stimulated cones get set up at the top of the vector on down,
    % depending on the spot size.
    stimulatedInDetectorSize = min([spotSize detectorSize]);
    backgroundMeanResponses = backgroundIsomerizations*ones(detectorSize,1);
    testMeanResponses = backgroundMeanResponses;
    testMeanResponses(1:stimulatedInDetectorSize) = testMeanResponses(1:stimulatedInDetectorSize) + x*ones(stimulatedInDetectorSize,1);
    
    % Compute fraction correct based on summing responses in detector
    fractionCorrect = analyticPoissonIdealObserver(sum(backgroundMeanResponses),sum(testMeanResponses));
    f = (fractionCorrect-criterionFractionCorrect).^2;
    
else
    nPositions = length(weights);
    backgroundIsomerizationsOverPosition = backgroundIsomerizations*ones(size(weights));
    backgroundMean = sum(backgroundIsomerizationsOverPosition.*weights);
    backgroundVar = sum((backgroundIsomerizationsOverPosition).*weights);
    
    testMeanResponsesOverPosition = backgroundIsomerizationsOverPosition;
    positionVals = (1:nPositions)-nPositions/2-0.5;
    [~,minIndex] = sort(abs(positionVals),'ascend');
    testPositionIndex = minIndex(1:spotSize);
    testMeanResponsesOverPosition(testPositionIndex) = testMeanResponsesOverPosition(testPositionIndex) + x*ones(size(1:spotSize));
    testMean = sum(testMeanResponsesOverPosition.*weights);
    testVar = sum((testMeanResponsesOverPosition).*weights);
    
    dPrime = (testMean-backgroundMean)/sqrt((testVar+backgroundVar)/2);
    fractionCorrect = dPrimeToTAFCFractionCorrect(dPrime);
    f = 1000*(fractionCorrect-criterionFractionCorrect).^2;
end


end