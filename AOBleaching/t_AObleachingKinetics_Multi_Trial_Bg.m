%% Bleaching kinetics in ISETBio
%
% Description:
%    An early version of this reproduced Figure S1 (panel D) of Cooper et al. (2020), showing how the
%    fraction of L/M cone pigment bleaches evolves over time with a
%    stimulus train of flashes.  That version is t_AOBleachingKinetics in
%    ISETBio.
%     
%    This version has evolved.
%
% See also: t_AODisplay, t_radianceToCornealIrradiance

% History:
%   01/09/21  dhb  Wrote it.
%   11/06/21  dhb  Allow superimposing stim+bg on stim plot.  Based on rw
%                  version.

% Run using the code example below. This is 
% a quick and dirty way to superimpose the three
% conditions on a single figure.
%{ 
    close all; 
    figure; clf; hold on
    
    clear; whichStimCondition = 'both';
    t_AObleachingKinetics_Multi_Trial_Bg;

    clear; whichStimCondition = 'stimulus';
    t_AObleachingKinetics_Multi_Trial_Bg;

    clear; whichStimCondition = 'background';
    t_AObleachingKinetics_Multi_Trial_Bg;
%}


%% Initialize
%close all; clear; ieInit;

%% Parameters


%% Wavelength support
deltaWl = 1;
wls = (500:deltaWl:900)';
nWls = length(wls);

%% Stimulus spectral properties
%
% Keep stimulus and background separate from each other
stimulusWl = 545;                           % Wavelength of narrowband background field. This can be the imaging light
stimulusCornealPowerNw = 153;               % Stimulus power passing through pupil in nanowatts
stimulusCornealPowerUw = stimulusCornealPowerNw*(1e-3);
stimulusFWHM = 5;                         % Full-width at half max of stimlus/background (in nm).
relSpd = normpdf(wls,stimulusWl,FWHMToStd(stimulusFWHM));
unitSpd = relSpd/trapz(wls,relSpd);

stimulusWlBg = 785;                        % Wavelength of narrowband background field. This can be the imaging light
stimulusCornealPowerUwBg = 90;             % Stimulus power passing through pupil in microwatts
stimulusFWHMBg = 5;                         % Full-width at half max of stimlus/background (in nm).
relSpdBg = normpdf(wls,stimulusWlBg,FWHMToStd(stimulusFWHMBg));
unitSpdBg = relSpdBg/trapz(wls,relSpdBg);

%% Stimulus spatial parameters
backgroundLinearSizeDegs = 1;              % Linear side of square field in degs

%% Make relative spectral power distributions.
%
% Approximated by a Gaussian with specified center wavlength and FWHM,
% and with total power given by the corneal power specified above.
% The call to trapz takes wavelength spacing into account when normalizing the power.

%% Pupil size
pupilDiameterMm = 7;                        % Pupil diameter.
pupilAreaMm = pi*(pupilDiameterMm/2)^2;     % For convenience below, compute pupil area.

%% Convert stimulus power to trolands
%
% Get equivalent spectral radiance as a function of the wavelength support.
%
% The routine here finds the radiance on an external conventional display that
% produces the same retinal illuminance as the corneal power specified
% above.  This is purely geometric calculation; attenuation of light by
% occular media is not taken into account at this stage.  Note that this
% conversion routine expects power per wavelength band, not power per nm,
% as its input, but returns power per nm as its output.  A little
% confusing, but that is why the spd being passed in is multiplied by the
% wavelength spacing.
stimulusSpdCornealPowerUw = stimulusCornealPowerUw*unitSpd;
stimulusSpdRadiance = AOMonochromaticCornealPowerToRadiance(wls,wls,stimulusSpdCornealPowerUw*deltaWl,pupilDiameterMm,backgroundLinearSizeDegs^2);

stimulusSpdCornealPowerUwBg = stimulusCornealPowerUwBg*unitSpdBg;
stimulusSpdRadianceBg = AOMonochromaticCornealPowerToRadiance(wls,wls,stimulusSpdCornealPowerUwBg*deltaWl,pupilDiameterMm,backgroundLinearSizeDegs^2);

% Make sure our computed radiance yields the desired corneal
% irradiance when we go in the other direction.  The magic
% numbers (1e6) and (1e-3) in the call just below do unit conversions
% from units we're using here to those expected by RadianceAndDegrees2ToCornIrradiance
stimulusSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(stimulusSpdRadiance,backgroundLinearSizeDegs^2)*(1e6)*((1e-3)^2);
stimulusCornealPowerUWCheck = trapz(wls,stimulusSpdCornealIrradianceUWMm2Check)*pupilAreaMm;
if (abs(stimulusCornealPowerUWCheck-stimulusCornealPowerUw)/stimulusCornealPowerUw > 1e-4)
    error('Do not get right cornal power back from computed radiance');
end

%% Convert from radiance to luminance.
%
% The magic number 683 makes the unis of luminance cd/m2, given the
% Watts/[sr-m2-nm] units we're using for radiance.
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,wls);
T_y = T_xyz(2,:);
stimulusLuminance = trapz(wls,T_y'.*stimulusSpdRadiance);
fprintf('Stimulus luminance is %0.1f cd/m2\n',stimulusLuminance);

stimulusLuminanceBg = trapz(wls,T_y'.*stimulusSpdRadianceBg);
fprintf('Background luminance is %0.1f cd/m2\n',stimulusLuminanceBg);

%% Convert luminance to trolands for bleaching calcs
stimulusTrolands = stimulusLuminance*pupilAreaMm;
fprintf('Stimulus is %0.1f trolands\n',stimulusTrolands);

stimulusTrolandsBg = stimulusLuminanceBg*pupilAreaMm;
fprintf('Background is %0.1f trolands\n',stimulusTrolandsBg);

% Compute a temporal stimulus timecourse
timePerTrialSec = 6;
stimStartTimeSec = 1;
timePerStimSec = 1;
nSets = 5; %Number of repeated trials
nTrials = 10;
Dark_Adaptation_Sec = 120;
totalTimeSec = (timePerTrialSec*nTrials+Dark_Adaptation_Sec);
totalTimeMsec = 1000*totalTimeSec;
timeMsec = (1:totalTimeMsec) - 1;
timeSec = timeMsec/1000;

switch (whichStimCondition)
    case 'stimulus'
        trolands = zeros(size(timeMsec));
        for ii = 1:nTrials
            startTime = (ii-1)*timePerTrialSec*1000 + stimStartTimeSec*1000;
            finishTime = startTime + timePerStimSec*1000;
            index = find(timeMsec > startTime & timeMsec < finishTime);
            trolands(index) = stimulusTrolands;
        end
        plotColor = 'g';

    case 'background'
        trolands = stimulusTrolandsBg*ones(size(timeMsec));
        plotColor = 'r';

    case 'both'
        trolands = stimulusTrolandsBg*ones(size(timeMsec));
        for ii = 1:nTrials
            startTime = (ii-1)*timePerTrialSec*1000 + stimStartTimeSec*1000;
            finishTime = startTime + timePerStimSec*1000;
            index = find(timeMsec > startTime & timeMsec < finishTime);
            trolands(index) = trolands(index) + stimulusTrolands;
        end
        plotColor = 'b';

    otherwise
        error('Unknown stimulus condition specified');
end

%% Compute bleaching over time
for jj = 1:nSets
    if jj == 1
        initialFractionBleached = 0;
        fractionBleached = ComputePhotopigmentBleaching(trolands,'cones','trolands','Boynton',initialFractionBleached,'msec');
        fractionUnbleached = 1 - fractionBleached;%Original
        Recovery = fractionBleached(1,totalTimeMsec);
        plot(timeSec,fractionUnbleached,plotColor,'LineWidth',2)
        fractionBleachedSets = zeros(nSets, length(fractionBleached));
        fractionUnbleachedSets = zeros(nSets, length(fractionUnbleached));
        fractionBleachedSets(jj,:) =  fractionBleached;
        fractionUnbleachedSets(jj,:) =  fractionUnbleached;
    else
        initialFractionBleached = Recovery;
        fractionBleached = ComputePhotopigmentBleaching(trolands,'cones','trolands','Boynton',initialFractionBleached,'msec');
        fractionUnbleached = 1 - fractionBleached;%Original
        Recovery = fractionBleached(1,totalTimeMsec);
        plot(timeSec+(totalTimeSec*(jj-1)),fractionUnbleached,plotColor,'LineWidth',2)
        fractionBleachedSets(jj,:) =  fractionBleached;
        fractionUnbleachedSets(jj,:) =  fractionUnbleached;
    end
end
   
% Finish up plot
xlim([0 totalTimeSec*jj])
ylim([0 1]);
xlabel('Time (sec)');
ylabel('Fraction L/M Cone Pigment Unbleached');






