% modelingPositionUncertainty
%
% This simulation is built using a simple one-dimensional detector array;
% consider each value as a single cone or midget RGC. In general, things
% are set up so that one always proceeds from the top of the vector down
% (we don't fuss about spatial arrangment here; we allow for some wrap
% around behavior in the case where the stimulus template is centered on
% each pixel in the uncertainty window). The template used to compute
% responses here is either (1) a single pixel, in which case the response
% on the test and background conditions is taken as the max value within
% the uncertainty window; or (2) a template that is matched to the stimulus
% size. These two approaches have different effects but both produce
% summation-like behavior. In both cases D-prime is computed from simulated
% response distributions across a range of stimulus intensities using
% Equation 6 from Geisler, Psychometric Functions of Uncertain Template
% Matching Observers, JOV (2018). Like Geisler, we inject some response
% noise into the simulation and use the resultant response distributions to
% determine d-Prime. We estimate threshold in each case by interpolating to
% find the stimulus intensity that produces a specified D-prime.

% 5-15-2018     wst wrote it

%% Housekeeping
close all
clear all
clc

%% Simulation parameters

numSimulations = 2000; % Number of simulations to run per parameter combination

stimulusMatchedTemplate = 1; % Set this to zero if you want to keep the stimulus template a single pixel; set to one if you want to have the stimulus template matched to the stimulus size;

numDetectors = 20; % Number of cones you are modeling; they will be assembled in an Nx1 vector
backgroundReceptorActivityLevel = 100; % Arbitrary; can think of this as mean isomerizations in each cone
noiseRMSpercent = 5; % Will add some noise to each receptor here, defined as RMS in % from background

positionUncertaintyWindowDiameters = [1 3 5 10]; % Uncertainty window diameters, in pixels; 1 = no uncertainty
if max(positionUncertaintyWindowDiameters) > numDetectors
    error('Max uncertainty window exceeds number of detectors');
end

stimSizes = 1:2:13; % Stimulus sizes, in pixels/number of cones; it is easier if these are kept as odd numbers
if max(stimSizes) > numDetectors
    error('Max stim size exceeds number of detectors');
end

targetAmplitudes = [0.1 0.5 1 1.5 2 2.5 3 4 6.5 8 10]; % Use these to interpolate to find a fixed level of dPrime
dPrimeTarget = 0.95; % Corresponds to 75% correct on 2AFC, per earlier simulations

%% Compute dPrime for each parameter (stim size, stim amplitude, and uncertainty window size) combination

% Pre-allocate
dPrimeMatrix = zeros(length(targetAmplitudes),length(stimSizes),length(positionUncertaintyWindowDiameters));
threshAmplitude = zeros(1,length(stimSizes),length(positionUncertaintyWindowDiameters));

% Loop through each parameter combination to determine d-Prime
for uw = 1:length(positionUncertaintyWindowDiameters)
    windowSize = positionUncertaintyWindowDiameters(uw);
    for ss = 1:length(stimSizes)
        stimSize = stimSizes(ss);
        for amp = 1:length(targetAmplitudes)
            targetAmplitude = targetAmplitudes(amp);
            % Pre-allocate
            responseVectorBackground = zeros(numSimulations,1);
            responseVectorTest = responseVectorBackground;
            % Repeat a bunch of times to produce distributions for
            % background and test cases
            for ns = 1:numSimulations
                % Make some noisy vectors
                receptorVectorBackground = backgroundReceptorActivityLevel+(noiseRMSpercent).*randn(numDetectors,1);
                receptorVectorTest = backgroundReceptorActivityLevel+(noiseRMSpercent).*randn(numDetectors,1);
                
                % Add on the stimulus amplitude to the test vector
                receptorVectorTest(1:stimSize) = receptorVectorTest(1:stimSize)+targetAmplitude;
                
                % Compute the max response to the stimulus within the
                % uncertainty window
                if stimulusMatchedTemplate == 1
                    stimMidPoint = ceil(stimSize/2); % This is the stimulus pixel we'll center at each point in the uncertainty window
                    % Build the stimulus template; zeros and ones
                    % should be fine, I think
                    stimTemplate = zeros(numDetectors,1);
                    stimTemplate(1:stimSize) = ones;
                    
                    % Compute dot product after centering stimulus
                    % template on every pixel in the uncerainty window,
                    % then sum the response
                    templateResponseVectorBackground = zeros(windowSize,1);
                    templateResponseVectorTest = zeros(windowSize,1);
                    for centerPixelLoc = 1:windowSize
                        stimTemplateShifted = circshift(stimTemplate,[centerPixelLoc-stimMidPoint 0]);
                        templateResponseVectorBackground(centerPixelLoc,1) = sum(stimTemplateShifted.*receptorVectorBackground);
                        templateResponseVectorTest(centerPixelLoc,1) = sum(stimTemplateShifted.*receptorVectorTest);
                    end
                    % Take the max of the template response vectors
                    responseVectorBackground(ns)=max(templateResponseVectorBackground);
                    responseVectorTest(ns) = max(templateResponseVectorTest);
                    
                elseif stimulusMatchedTemplate == 0 % Just look for the highest single-cone (i.e. pixel) response in the uncertainty window
                    responseVectorBackground(ns)=max(receptorVectorBackground(1:windowSize));
                    responseVectorTest(ns) = max(receptorVectorTest(1:windowSize));
                end
            end
            % Compute d-Prime (Eq. 6 from Geisler, JOV (2018))
            dPrimeNumerator = (mean(responseVectorTest)-mean(responseVectorBackground)).*sqrt(2);
            dPrimeDenominator = sqrt(var(responseVectorTest)+var(responseVectorBackground));
            dPrimeMatrix(amp,ss,uw) = dPrimeNumerator./dPrimeDenominator;
        end
        % For each stimulus size, interpolate to determine the test amplitude that produces desired dPrime
        threshAmplitude(1,ss,uw) = interp1(dPrimeMatrix(:,ss,uw),targetAmplitudes', dPrimeTarget);
    end
end

%% Plot the results
figure, hold on
set(gcf, 'Units', 'inches', 'Position', [0.5 0.5 12 6], 'Color', [1 1 1])
theColors = [1 0 0; 0 0.7 0; 0 0 1; 0.7 0.7 0];
colorScaling = linspace(0.25,1,length(stimSizes));

for uw = 1:length(positionUncertaintyWindowDiameters)
    for ss = 1:length(stimSizes)
        if uw < 3
            subplot(2,4,uw)
        else
            subplot(2,4,uw+2)
        end
        hold on, plot(targetAmplitudes,dPrimeMatrix(:,ss,uw), '-o','Color', colorScaling(ss).*theColors(uw,:), 'MarkerFaceColor', colorScaling(ss).*theColors(uw,:),'MarkerSize',4, 'LineWidth', 2)
    end
    xlabel('Target amplitude (au)');
    ylabel('D-prime');
    axis square
    ylim([0 ceil(max(dPrimeMatrix(:)))]);
    if positionUncertaintyWindowDiameters(uw) == 1
        text(0.25,ceil(max(dPrimeMatrix(:))),'No uncertainty', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    else
        text(0.25,ceil(max(dPrimeMatrix(:))),[num2str(positionUncertaintyWindowDiameters(uw)) ' px diameter'], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end
end

% Now plot stimulus area vs energy on log-log axes
for uw = 1:length(positionUncertaintyWindowDiameters)
    subplot(2,4,[3 4 7 8]), hold on
    plot(log10(stimSizes), log10(stimSizes)+log10(threshAmplitude(:,:,uw)), '-o','Color', theColors(uw,:), 'MarkerFaceColor', theColors(uw,:),'MarkerSize',10, 'LineWidth', 2)
    if positionUncertaintyWindowDiameters(uw) == 1
        legendText{uw,1} = 'No uncertainty'; %#ok<SAGROW>
    else
        legendText{uw,1} = ['Uncertainty window: ' num2str(positionUncertaintyWindowDiameters(uw))]; %#ok<SAGROW>
    end
end

% Set the legend
hLegend = legend(legendText,'Location', 'NorthWest', 'AutoUpdate', 'off');
box on
axis equal
grid on

% Plot indicators of the uncertainty zone along the x-axis
yLimMin = min(get(gca, 'YLim'));
yLimMax = max(get(gca, 'YLim'));
for uw = 1:length(positionUncertaintyWindowDiameters)
    subplot(2,4,[3 4 7 8]), hold on
    plot([log10(positionUncertaintyWindowDiameters(uw)) log10(positionUncertaintyWindowDiameters(uw))], [yLimMin yLimMin+0.1], '-', 'Color', theColors(uw,:), 'LineWidth', 4)
end
ylim([yLimMin yLimMax]);
xlabel('Log stimulus area (au)', 'FontSize', 12);
ylabel('Log threshold energy (au)', 'FontSize', 12);
title(['RMS noise: ' num2str(noiseRMSpercent) '%; numSimulations: ' num2str(numSimulations)], 'FontSize', 12)

if stimulusMatchedTemplate == 1
    fileName = ['PositionUncertaintySummation_StimMatchedTemplate_rms_' num2str(noiseRMSpercent) '_nSim_' num2str(numSimulations) '.pdf'];
else
    fileName = ['PositionUncertaintySummation_SinglePixelTemplate_rms_' num2str(noiseRMSpercent) '_nSim_' num2str(numSimulations) '.pdf'];
end

% Save the figure
FigureSave(fileName,gcf,'pdf');

