function SimpleSpatialSummation2D

%% Clear
clear; close all; clc

%% Parameters
%
% These are dimensions now linear diameters; for now, consider 1 pixel = 1
% arcmin
minSpotSize = 1;
nSpotSizes = 30;
maxSpotSize = 50;
nPositions = 100;
if (rem(nPositions,2) ~= 0)
    error('nPositions must be even');
end
detectorSizes = [6 12 18]; % Human foveal parasol cell diamter ~30 microns (about 0.1 deg or 6 arcmin): Dacey and Petersen (1992)

% Detector array specifications
gaussianFlag = 1; % Set to 1 if want detecting units to have Gaussian profile
centerSpacings = 2*detectorSizes; % Make these equal to detector sizes for a tightly-packed array
imageSize = maxSpotSize*3; % Pad things out a bit

% Background retinal irradiance (isomerizations per cone per stimulus
% duration)
backgroundIsomerizations = 100;

% Criterion fraction correct
criterionFractionCorrect = 0.78;

%% Poisson calculation
%
% Loop over spot sizes and find threshold for each
weights = [];
spotSizes = unique(round(logspace(log10(minSpotSize),log10(maxSpotSize),nSpotSizes)));

% Compute the number of pixels in the spot; to be used for more precise energy plotting since small circles might be drawn as squares
spotPixels = zeros(size(spotSizes));

% Loop through the detector sizes
for dd = 1:length(detectorSizes)
    detectorSize = detectorSizes(dd);
    
    % Get hexagonal coordinates based on the detector size and spacing
    [hexCenters] = makeHexagonalArray(centerSpacings(dd), imageSize);
    
    % Make the detector array
    [detectorImageLayers] = makeDetectorImage(hexCenters, centerSpacings(dd), detectorSize, imageSize, gaussianFlag);
    
    numIterations = length(spotSizes);
    bump = 1;    
    h = waitbar(0, ['Running simulation for detector size ' num2str(detectorSize)]);
    
    for ss = 1:length(spotSizes)
        
        % Set spot size
        spotSize = spotSizes(ss);
        
        % Make an image of the 2D stimulus to be passed to the
        % threshold-finding algorithm; everything that follows is centered
        % on the center pixel
        centerLoc = imageSize/2;
        spotImage = zeros(imageSize, imageSize);
        backgroundImage = spotImage+backgroundIsomerizations;
        
        % Place the stimulus in the blank canvas
        startLoc = centerLoc-floor(spotSize/2);
        stopLoc = startLoc+spotSize-1;
        spotImage(startLoc:stopLoc, startLoc:stopLoc) = double(Circle(spotSize/2));      
        if (any(isnan(spotImage(:))))
            error('Got a NaN in spot image after Circle call');
        end
        
        if dd == 1 % Count up the number of pixels the first time through;
            spotPixels(ss) = nansum(spotImage(:));
        end
        % Sub-pixel shifting, if necessary
        if rem(spotSize,2)==0 % Spot size even
            % Shift spotImage back onto centerLoc pixel
            [xx, yy] = meshgrid(1:imageSize, 1:imageSize);
            spotImage = interp2(xx,yy,spotImage, xx-0.5, yy-0.5,'linear',0);
        end   
        if (any(isnan(spotImage(:))))
            error('Got a NaN in spot image after interp');
        end
        
        % Threshold guesses and bounds in irradiance expressed as
        % isomerizations per cone per stimulus duration
        
        % Setting the initial threshold guess seems to be a bit finnicky 
        if (ss == 1)
            initialThreshold = 1;
        else % Set the next threshold guess by assuming linear summation holds
            initialThreshold = (threshold(dd,ss-1).*(spotSizes(ss-1).^2))./(spotSize.^2);
        end
        maxThreshold = 100000;
        minThreshold = 0;
        
        % Use fmincon to find threshold irradiance
        options = optimset('fmincon');
        options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
        threshold(dd,ss) = fmincon(@(x)FindThresholdFunction(x,backgroundIsomerizations,detectorImageLayers,spotImage,criterionFractionCorrect,weights),initialThreshold,[],[],[],[],minThreshold,maxThreshold,[],options);
        thresholdEnergy(dd,ss) = pi*((spotSize/2).^2)*threshold(dd,ss);
        thresholdEnergyPixels(dd,ss) = spotPixels(ss).*threshold(dd,ss);
        if isfinite(threshold(dd,ss))
            [~,fractionCorrect(dd,ss)] = FindThresholdFunction(threshold(dd,ss),backgroundIsomerizations,detectorImageLayers,spotImage,criterionFractionCorrect,weights);
            % Flag the data that get close to the criterion fraction
            % correct; don't plot data where the algorithm got stuck
            if abs(fractionCorrect(dd,ss)-criterionFractionCorrect)<0.015
                plotFlag(dd,ss) = 1;
            else
                plotFlag(dd,ss) = 0;
            end
        else
            fractionCorrect(dd,ss) = NaN;
            plotFlag(dd,ss) = 0;
        end

        % Update the waitbar
        waitbar(bump/numIterations,h)
        bump = bump+1;
        
        % NOW FIGURE OUT HOW TO COMBINE MULTIPLE PFs TO DETERMINE OVERALL
        % DETECTION PROBABILITY
        
        fprintf('Detector size %d, spot size %d, threshold irradiance %0.2f, percent correct: %0.2f\n',detectorSize,spotSize,threshold(dd,ss),round(100*fractionCorrect(dd,ss)));
    end
    close(h);
end

%% Plot
theColors = ['c' 'r' 'g' 'b' 'k' 'y'];
theFig = figure; clf; hold on;
set(gcf,'Position',[100 100 1050 560]);
% subplot(1,2,1); hold on;
% for dd = 1:length(detectorSizes)
%     theColor = theColors(rem(dd,length(theColors)) + 1);
%     plot(log10(pi.*((spotSizes(plotFlag(dd,:)==1)./2).^2)),log10(threshold(dd,plotFlag(dd,:)==1)),[theColor 'o'],'MarkerSize',10);
%     plot(log10(pi.*((spotSizes(plotFlag(dd,:)==1)./2).^2)),log10(threshold(dd,plotFlag(dd,:)==1)),[theColor ],'LineWidth',2);
%     plot([log10(pi.*((detectorSizes(dd)./2).^2)) log10(pi.*((detectorSizes(dd)./2).^2))],[1 1.25],[theColor ],'LineWidth',3);
% end
xlim([0 3.5]);
% ylim([1 4.5]);
axis('square');
% xlabel('Log Stimulus Area');
% ylabel('Log Threshold Irradiance');

riccosAreaFits = zeros(size(detectorSizes));

% Fit this panel with the two-segment linear regression
subplot(1,2,1); hold on;
for dd = 1:length(detectorSizes)
    theColor = theColors(rem(dd,length(theColors)) + 1);
    % Plot data in number of pixels
    plot(log10(spotPixels(plotFlag(dd,:)==1)),log10(thresholdEnergyPixels(dd,plotFlag(dd,:)==1)),[theColor 'o'],'MarkerSize',10);
    plot([log10(pi.*((detectorSizes(dd)./2).^2)) log10(pi.*((detectorSizes(dd)./2).^2))],[2.5 2.75],[theColor ],'LineWidth',3);
    
    % Fit the data via two-segment linear regression
    xData = log10(spotPixels(plotFlag(dd,:)==1));
    yData = log10(thresholdEnergyPixels(dd,plotFlag(dd,:)==1));
    beta = [min(yData) mean(xData) 1]; % Starting fit parameters: [Branch1_height Riccos_intercept Branch2_slope]
    riccosFit = fitnlm(xData,yData,@riccosFittingFunction, beta);
    xEval = linspace(min(xData), max(xData), 1000)';
    [yEval] = riccosFittingFunction(riccosFit.Coefficients.Estimate,xEval);
    plot(xEval,yEval, '-', 'Color', theColor, 'LineWidth', 2);
    riccosAreaFits(dd) = riccosFit.Coefficients.Estimate(2);
    plot([riccosFit.Coefficients.Estimate(2) riccosFit.Coefficients.Estimate(2)], [2.5 riccosFit.Coefficients.Estimate(1)], ':', 'Color', theColor', 'LineWidth', 2);
    
    fprintf('Slope of rising limb is %0.2f\n',riccosFit.Coefficients.Estimate(3));

end
xlim([0 3.5]);
ylim([2.5 6]);
axis('square');
xlabel('Log Stimulus Area');
ylabel('Log Threshold Energy');
% FigureSave('SimplePoissonSpatialSummation',gcf,'pdf');

ylim([0 3.5])
axis('equal'); grid on


end

%% Error function.
function [f,fractionCorrect] = FindThresholdFunction(x,backgroundIsomerizations,detectorImage,spotImage,criterionFractionCorrect,weights)

% If weights is empty, assume a hard summation error and do the Poisson
% calculation.  Otherwise weight up the cones and do a Gaussian
% calculations.
if (isempty(weights))
    % Build up background and test stimuli as vectors.  We don't model the
    % spatial geometry, just count cones.  The vector just specifies the
    % cones within the detector, becaues those are the only ones we pay
    % attention to. The stimulated cones get set up at the top of the
    % vector on down, depending on the spot size.
    backgroundResponses = zeros(size(detectorImage,3),1);
    testResponses = backgroundResponses;
    stimImage = x.*spotImage; % We are solving for x here by minimizing "f" below
    
    % Each detector resides in its own layer; loop through them here to
    % compute the background and test isomerizations for each detector and
    % compile into vectors "backgroundResponses" and "testResponses"
    for layerNum = 1:size(detectorImage,3)
        detectorResponseBackground = detectorImage(:,:,layerNum).*backgroundIsomerizations;
        backgroundResponses(layerNum) = sum(detectorResponseBackground(:));
        detectorResponseTest = (stimImage+backgroundIsomerizations).*detectorImage(:,:,layerNum);
        testResponses(layerNum) = sum(detectorResponseTest(:));
    end
    
    % Check for horrible things
    if (any(isnan(testResponses)))
        error('Got a NaN in test responses');
    end
    if (any(isnan(backgroundResponses)))
        error('Got a NaN in background responses');
    end
    
    % Do the computation
    efficiency = 0.5;
    fractionCorrect = analyticPoissonIdealObserver(efficiency*backgroundResponses,efficiency*testResponses);
    f = 1000.*(fractionCorrect-criterionFractionCorrect).^2;
    
else % This section is carried over from the 1D simulation but currently not used
    nPositions = length(weights);
    backgroundIsomerizationsOverPosition = backgroundIsomerizations*ones(size(weights));
    backgroundMean = sum(backgroundIsomerizationsOverPosition.*weights);
    backgroundVar = sum((backgroundIsomerizationsOverPosition).*weights);
    
    testMeanResponsesOverPosition = backgroundIsomerizationsOverPosition;
    positionVals = (1:nPositions)-nPositions/2-0.5;
    [~,minIndex] = sort(abs(positionVals),'ascend');
    testPositionIndex = minIndex(1:spotImage);
    testMeanResponsesOverPosition(testPositionIndex) = testMeanResponsesOverPosition(testPositionIndex) + x*ones(size(1:spotImage));
    testMean = sum(testMeanResponsesOverPosition.*weights);
    testVar = sum((testMeanResponsesOverPosition).*weights);
    
    dPrime = (testMean-backgroundMean)/sqrt((testVar+backgroundVar)/2);
    fractionCorrect = dPrimeToTAFCFractionCorrect(dPrime);
    f = 1000*(fractionCorrect-criterionFractionCorrect).^2;
end
end

%% Detector coordinates function
function [hexCenters] = makeHexagonalArray(centerSpacing, imageSize)

% This function generates x,y coordinates arranged in a hexagonal grid;
% these locations can be passed to another function that builds up an array
% of detectors whose locations are defined by these coordinates
%
% Input variables:
% centerSpacing     center-to-center spacing of the detectors, in pixels
% imageSize         size of the image, in pixels

if (rem(imageSize,2) ~= 0)
    error('imageSize must be even');
end

% Make sure things are centered in the image
centerLoc = imageSize/2;
minLoc = 0;
maxLoc = imageSize;

% Generate hexagonal grid
[X, Y] = meshgrid(minLoc:centerSpacing:maxLoc*(2/sqrt(3)));
X = (sqrt(3)/2).*X;
n = size(X,1);
if rem(n,2) == 0
    Y = Y + repmat([0 centerSpacing/2],[n,floor(n/2)]);
else
    Y = Y + [repmat([centerSpacing/2 0],[n,floor(n/2)]) ones(n,1).*centerSpacing/2];
end

% Shift X and Y locations to be centered on centerLoc
[~,closestColumn] = find(abs(X-centerLoc) == min(abs(X(:)-centerLoc)));
colShift = centerLoc-X(1,closestColumn(1,1));
X = X+colShift; % Apply the shift
[closestRow,~] = find(abs(Y(:,closestColumn(1,1))-centerLoc) == min(abs(Y(:,closestColumn(1,1))-centerLoc)));
rowShift = centerLoc-Y(closestRow(1,1),closestColumn(1,1));
Y = Y+rowShift; % Apply the shift


% Trim locations close to the edges of the image so we pass something sensible back
X(X<centerSpacing/2) = NaN;
Y(Y<centerSpacing/2) = NaN;
X(X>imageSize-centerSpacing/2) = NaN;
Y(Y>imageSize-centerSpacing/2) = NaN;

hexCenters(:,:,1) = X;
hexCenters(:,:,2) = Y;
end

%% Detector array assembly function

function [detectorImageLayers] = makeDetectorImage(hexCenters, centerSpacing, detectorSize, imageSize, gaussianFlag)

% This function should build a representation of the detector array in the
% form of [imageSize x imageSize x nDetectors]. There are likely practical
% advantages to keeping each detector in a separate layer, but detector
% array can be visually represented as a single image as follows:
% detectorImageSingle = max(detectorImageLayers, [],3); or:
% detectorImageSingle = sum(detectorImageLayers, 3);
%
% Input variables:
% hexCenters        two-layer matrix containing the X and Y coordinates of the detectors
% centerSpacing     center-to-center spacing of the detectors
% detectorSize      Size of detector: diameter for a circle; FHWM for a Gaussian detector
% imageSize         Desired size of detector array
% gaussianFlag      Flag for Gaussian-shaped detector (set to 0 for circular detectors with hard edges)

if (rem(imageSize,2) ~= 0)
    error('imageSize must be even');
end

if centerSpacing<detectorSize
    warning('Detector size exceeds detector spacing');
end

centerLoc = imageSize/2;

% Center the detector in the center of the image canvas
if gaussianFlag == 0 % Hard-edged circle
    if floor(detectorSize)~=detectorSize
        error('Circle dimensions must be an integer')
    end
    detectorImageCropped = double(Circle(detectorSize/2)); % Requires Circle.m from Psychtoolbox
    detectorImage = zeros(imageSize, imageSize);
    startRow = imageSize/2-floor(detectorSize/2);
    stopRow = imageSize/2-floor(detectorSize/2)+detectorSize-1;
    detectorImage(startRow:stopRow,startRow:stopRow) = detectorImageCropped;
    if rem(detectorSize,2) == 0
        % Shift detectorImage back onto centerLoc pixel
        [xx, yy] = meshgrid(1:imageSize, 1:imageSize);
        detectorImage = interp2(xx,yy,detectorImage, xx-0.5, yy-0.5);
    end
    
elseif gaussianFlag == 1 % Gaussian detector
    % FWHM = 2*sqrt(2*log(2))*sigma
    sigma = detectorSize/(2*sqrt(2*log(2)));
    detectorImage = fspecial('gaussian', imageSize+1, sigma); % Making this odd-sized ensures the peak is confined to a single pixel
    detectorImage(1,:) = []; % Trim the edges to get back to imageSize
    detectorImage(:,1) = []; % Trim the edges to get back to imageSize
    detectorImage = detectorImage./max(detectorImage(:)); % Normalize peak to 1
end

% Now shift the template detector to each [x,y] location
bump = 1;
X = hexCenters(:,:,1);
Y = hexCenters(:,:,2);
detectorImageLayers = zeros(imageSize, imageSize, sum(sum(isfinite(X.*Y))));
[xx, yy] = meshgrid(1:imageSize, 1:imageSize);
for n = 1:size(X,1)
    for j = 1:size(X,2)
        xShift = centerLoc-(X(n,j));
        yShift = centerLoc-(Y(n,j));
        if isfinite(xShift) && isfinite(yShift)
            detectorImageTemp = interp2(xx,yy,detectorImage, xx+xShift, yy+yShift);
            detectorImageTemp(isnan(detectorImageTemp)) = 0;
            detectorImageLayers(:,:,bump) = detectorImageTemp;
            bump = bump+1;
        end
    end
end
end

%% Ricco's fitting function
function [yEval]  = riccosFittingFunction(beta,xData)
% Fit a function to (IxA) vs A (in log units) where two lines intersect. The first branch has a slope of
% zero, while the slope of the second in unconstrained
%
% Input variables:
% xData     log stimulus areas
% beta      Starting fit parameters [Branch1_Height Branch_Intercept Branch2_Slope]

b1 = beta(1); % First branch height
b2 = beta(2); % Segment intercept (i.e. Ricco's area)
b3 = beta(3); % Second branch slope

% The function takes this form:
yEval = [0.*xData(xData<=b2)+b1; b3.*xData(xData>b2)+(b1-(b2*b3))];
end
