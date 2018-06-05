% simpleConeReflectanceModel  Compute reflectance versus time for a swelling cone
%
% Description:
%   Very simple calculation of how light reflected from a cone varies with
%   cone length, taking into account both incoherent reflection and
%   interference between coherent component of light.
%
%   Based on a very simple model where there are two reflecting surfaces.  Light
%   reflects from each surface and interferes at the source side of the first surface.
%
%   The fraction of coherent (interfering) versus incoherent light is simply specified,
%   not computed.  Adding this computation would be a good extension.
%
%   The model has a number of parameters that have physical analogs. These
%   are set and commented in the code below.
%
% Notes:
% * [NOTE: DHB - I tried to understand coherence length from the web. In some places it
%   is defined a length at which phase is totally random.  In others, a
%   more graded definition is used.  We need to understand which the one we
%   calculate represents, if we are to add coherence length into this
%   cacluation.]
%

% History:
%   10/28/17  dhb  Wrote it.

%% Clear
clear; close all;

%% Parameters
%
% Imaging light wavelength
firstLightWavelengthNm = 720;
secondLightWavelengthNm = 1.25*firstLightWavelengthNm; % 850;
refractiveIndex1 = 1.331;
refractiveIndex2 = 1.328;

% Number of cones to simulate
nCones = 400;

% Length between reflecting surfaces
reflectingSurfaceDifferenceNm1 = 10*1e3;
reflectingSurfaceDifferenceNm2 = 7*1e3;

% How much variation in the starting length between surfces is there
% across cones (or trials within measurements for one cone), experessed as a
% fraction of the length between the two reflecting surfaces.
reflectingSurfaceRandomizationFraction = 0.1;

% Maximum length change between reflecting surfaces,
% in response to stimulus.
maxDeltaLengthNm = 150;

% What fraction of a full half sinusoid does length vary over
% in response to a stimulus, and how much trial-to-trial jitter
% is there in that fraction.
%
% Set the fraction to 1 if you think what happens is that the cone
% swells and then contracts to about its original length.  Set the
% fraction to 0.5 if you think it swells and then sits at the swollen
% length over the relevant time frame of the experiment.
%
% The randomization is also expressed as a fraction of the full half
% sinusoid
halfSinusoidFraction = 0.8;
halfSinusoidFractionRandomizationFraction = 0.5;

% Time in seconds for stimulus to have its effect. This just puts
% a time axis on the plots.
stimEffectTimeSecs = 4;

% Number sampling timepoints for simulation of what happens on a trial
nDeltaLengths = 100;
theDeltaTimesSecs = linspace(0,stimEffectTimeSecs,nDeltaLengths);

% Fraction of light reflected from the first and second surfaces.  Model below
% assumes that light not reflected from first surface is transmitted, and
% that light reflected from the second surface passes through the first
% surface both on the way in and on the way out.
firstSurfaceReflectance = 0.2;
secondSurfaceReflectance = 0.4;
thirdSurfaceReflectance = 0.4;

% Coherent fraction
%
% This depends on coherence length relative relative to distance between
% reflecting surfaces Here just making it up, because I don't really
% understand coherence length.
coherentFraction = 0.5;

%% Get intensities of light coming back from each surface
%
% This is calculated at the source side of the first surface; the light
% from the second surface goes through the first surface twice.
firstSurfaceIntensity = firstSurfaceReflectance;
secondSurfaceIntensity = (1-firstSurfaceReflectance)^2*secondSurfaceReflectance;

%% Different starting lengths, drawn at random from a uniform distribution
%
% Distribution is centered on the assumed length
randValues = rand(1,nCones);
theStartingDeltas = reflectingSurfaceRandomizationFraction*reflectingSurfaceDifferenceNm1*(randValues-0.5);

%% Set first phase to 0 without loss of generality
phase0 = 0;

%% Loop over cones/trials
for jj = 1:nCones
    % For each, we choose a random starting length and compute the
    % interferece effect. We also choose a random excursion, jittered
    % around a half sinusoid
    %
    % Get delta lengths as a function of time.
    theDeltaLengthsNm{jj} = maxDeltaLengthNm*sin(halfSinusoidFraction*(1+halfSinusoidFractionRandomizationFraction*(rand(1,1)-0.5))*pi*(0:nDeltaLengths-1)/nDeltaLengths);
    for ii = 1:nDeltaLengths
        
        % Get starting distance between surfaces for this cone, and compute phase of
        % second light at the point of interference.
        distanceNm1(jj,ii) = reflectingSurfaceDifferenceNm1 + theStartingDeltas(jj) + theDeltaLengthsNm{jj}(ii);
        distanceNm2(jj,ii) = reflectingSurfaceDifferenceNm2 + theStartingDeltas(jj) + theDeltaLengthsNm{jj}(ii);

        phaseWl1Sur1(jj,ii) = 2*pi*(2*distanceNm1(jj,ii))*refractiveIndex1/firstLightWavelengthNm;
        phaseWl2Sur1(jj,ii) = 2*pi*(2*distanceNm1(jj,ii))*refractiveIndex2/secondLightWavelengthNm;
        
        phaseWl1Sur2(jj,ii) = 2*pi*(2*distanceNm2(jj,ii))*refractiveIndex1/firstLightWavelengthNm;
        phaseWl2Sur2(jj,ii) = 2*pi*(2*distanceNm2(jj,ii))*refractiveIndex2/secondLightWavelengthNm;
        
        % Get amplitude of coherent component of reflected light, taking interference
        % into account.
        %
        % Formula from Hecht, Optics, 3rd ed, p. 291
        coherentAmplitude1(jj,ii) = sqrt(...
            firstSurfaceIntensity^2 + secondSurfaceIntensity^2 + ...
            2*firstSurfaceIntensity*secondSurfaceIntensity*cos(phase0-phaseWl1Sur1(jj,ii)));
        coherentAmplitude2(jj,ii) = sqrt(...
            firstSurfaceIntensity^2 + secondSurfaceIntensity^2 + ...
            2*firstSurfaceIntensity*secondSurfaceIntensity*cos(phase0-phaseWl2Sur1(jj,ii)));
        
        % Put together coherent and incoherent reflected light to get
        % intensity of total reflected light.  This is a model of what
        % we measure.
        totalAmplitude1(jj,ii) = coherentFraction*coherentAmplitude1(jj,ii) + (1-coherentFraction)*(firstSurfaceIntensity + secondSurfaceIntensity);
        totalAmplitude2(jj,ii) = coherentFraction*coherentAmplitude2(jj,ii) + (1-coherentFraction)*(firstSurfaceIntensity + secondSurfaceIntensity);
        
        
    end
    
    % Get initial slope of responses at both wavelengths
    nTimes = 50;
    initialSlope1(jj) = (totalAmplitude1(jj,nTimes+1)-totalAmplitude1(jj,1))/(theDeltaTimesSecs(nTimes+1)-theDeltaTimesSecs(1));
    initialSlope2(jj) = (totalAmplitude2(jj,nTimes+1)-totalAmplitude2(jj,1))/(theDeltaTimesSecs(nTimes+1)-theDeltaTimesSecs(1));
    
end

%% Make a plot of responses
nConesToPlot = 9;
theColors = ['r' 'g' 'b' 'k' 'y' 'c'];
nColors = length(theColors);
responseFig1 = figure; clf; hold on
responseFig2 = figure; clf; hold on
correlationFig = figure; clf; hold on
whichColor = 1;
for jj = 1:nConesToPlot
    figure(responseFig1);
    plot(theDeltaTimesSecs,totalAmplitude1(jj,:),theColors(whichColor),'LineWidth',3);
    
    figure(responseFig2);
    plot(theDeltaTimesSecs,totalAmplitude2(jj,:),theColors(whichColor),'LineWidth',3);

    figure(correlationFig);
    plot(totalAmplitude1(jj,:),totalAmplitude2(jj,:),'o','Color',theColors(whichColor),'MarkerFaceColor',theColors(whichColor), ...
        'MarkerSize',8);
    
    whichColor = whichColor + 1;
    if (whichColor > nColors)
        whichColor = 1;
    end
end
figure(responseFig1);
set(gca,'YTick',[0.2 0.3 0.4 0.5]);
set(gca,'FontName','Helvetica','FontSize', 16);
xlabel('Time','FontName','Helvetica','FontSize',18);
ylabel('Reflectance','FontName','Helvetica','FontSize',18);
xlim([0 2]);
ylim([0.2 0.5]);
set(gca,'XTickLabel',''); set(gca,'YTickLabel','');
FigureSave('ExampleResponses',responseFig1,'tif');
figure(responseFig2);
xlabel('Time (secs)');
ylabel('Reflectance');
ylim([0.2 0.6]);
figure(correlationFig);
set(gca,'FontName','Helvetica','FontSize', 16);
xlabel(sprintf('Standardized Reflectance %d nm',firstLightWavelengthNm),'FontName','Helvetica','FontSize',18);
ylabel(sprintf('Standardized Reflectance %d nm',secondLightWavelengthNm),'FontName','Helvetica','FontSize',18);
axis('square');
xlim([0.2 0.5]);
set(gca,'XTick',[0.2 0.3 0.4 0.5]);
set(gca,'YTick',[0.2 0.3 0.4 0.5]);
ylim([0.2 0.5]);
FigureSave('TwoWavelengthCompare',correlationFig,'tif');

%% Make a plot of initial slopes at two wavelengths
slopeFig = figure; clf;
plot(initialSlope1,initialSlope2,'ro','MarkerFaceColor','r','MarkerSize',14);
xlabel(sprintf('Initial response slope (%d nm)',firstLightWavelengthNm));
ylabel(sprintf('Initial response slope (%d nm)',secondLightWavelengthNm));
xlim([-0.3 0.3]); ylim([-0.3 0.3]); axis('square');

%% Make a plot of reflectance versus reflectance
slopeFig1 = figure; clf; hold on;
plot(totalAmplitude1(:,1),initialSlope1,'ro','MarkerFaceColor','r','MarkerSize',8);
plot(totalAmplitude2(:,1),initialSlope2,'bo','MarkerFaceColor','b','MarkerSize',8);
xlabel(sprintf('Initial reflectance amplitude',firstLightWavelengthNm));
ylabel(sprintf('Initial response slope',secondLightWavelengthNm));
xlim([0.2 0.6]); ylim([-0.3 0.3]); 
legend({sprintf('%d nm',round(firstLightWavelengthNm)),sprintf('%d nm',round(secondLightWavelengthNm))},'Location','NorthEast')

