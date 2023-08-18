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
% conditions on a single figure.  Matched up to
% typical conditions in the Penn PCAM AOSLO.
%{
    % Initialize
    close all; clear;
    trolandFigure = figure; clf;
    set(gcf,'Position',[200 420 1200 660]);
    subplot(1,2,1); hold on

    % Stimulus energy from lowest power condition in our repeatability
    % paper.
    stimulusWl = 545; 
    stimulusCornealPowerNw = 153;
    bgWl = 785;
    bgCornealPowerUw = 90; 
    stimStartTimeSec = 1;
    timePerStimSec = 1;
    title(sprintf('Stim: %d nm, %0.1f nW, %0.4f sec, %0.1f nJ',stimulusWl,stimulusCornealPowerNw,timePerStimSec,stimulusCornealPowerNw*timePerStimSec)); 

    % Run for all three stimulus conditions
    whichStimCondition = 'both';
    t_AObleachingKinetics_Multi_Trial_Bg;

    whichStimCondition = 'stimulus';
    t_AObleachingKinetics_Multi_Trial_Bg;

    whichStimCondition = 'background';
    t_AObleachingKinetics_Multi_Trial_Bg;

    % Same stimulus energy with different time course
    figure(trolandFigure);
    subplot(1,2,2); hold on
    stimulusCornealPowerNw = 2750;
    stimStartTimeSec = 1;
    timePerStimSec = 153/2750;
    title(sprintf('Stim: %d nm, %0.1f nW, %0.4f sec, %0.1f nJ',stimulusWl,stimulusCornealPowerNw,timePerStimSec,stimulusCornealPowerNw*timePerStimSec)); 

    % Run again
    whichStimCondition = 'both';
    t_AObleachingKinetics_Multi_Trial_Bg;

    whichStimCondition = 'stimulus';
    t_AObleachingKinetics_Multi_Trial_Bg;

    whichStimCondition = 'background';
    t_AObleachingKinetics_Multi_Trial_Bg;
%}

% Zhang et al., 2019
%
% Our calculations match up well with the values in the paper.
%
% Note that we verified with Don Miller that their stimulus was a circular spot,
% not a square field.
%{
    % Initialize
    close all; clear;
    trolandFigure = figure; clf;
    set(gcf,'Position',[200 420 1200 1200]);

    % Linear size to produce same area with a square as their 2 deg circular stimulus
    diameterDegs = 2;
    linearSizeDegs = sqrt(pi)*diameterDegs/2;

    % Timing
    nVideos = 1;
    nSets = 1;
    timePerVideoSec = 30;

    % Stimulus wavelength info
    whichStimCondition = 'stimulus';

    % 3% bleach
    stimulusWl = 637; 
    stimulusFWHM = 3; 
    figure(trolandFigure);
    subplot(2,2,1); hold on
    timePerStimSec = 5e-3;
    stimulusCornealPowerNw = 1e3*0.53/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-3/100*[1 1],'b:','LineWidth',2);
    fprintf('Them %0.1f%%, us %0.1f%%\n\n',3,100*max(fractionBleached));

    % 16.7% bleach
    stimulusWl = 637; 
    stimulusFWHM = 3;
    figure(trolandFigure);
    subplot(2,2,2); hold on
    timePerStimSec = 10e-3;
    stimulusCornealPowerNw = 1e3*3.2/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-16.7/100*[1 1],'b:','LineWidth',2);
    fprintf('Them %0.1f%%, us %0.1f%%\n\n',16.7,100*max(fractionBleached));

    % 9.9% bleach
    stimulusWl = 528; 
    stimulusFWHM = 3;
    figure(trolandFigure);
    subplot(2,2,3); hold on
    timePerStimSec = 5e-3;
    stimulusCornealPowerNw = 1e3*0.5/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-9.9/100*[1 1],'b:','LineWidth',2);
    fprintf('Them %0.1f%%, us %0.1f%%\n\n',9.9,100*max(fractionBleached));

    % 1.1% bleach
    stimulusWl = 450; 
    stimulusFWHM = 20;
    figure(trolandFigure);
    subplot(2,2,4); hold on
    timePerStimSec = 5e-3;
    stimulusCornealPowerNw = 1e3*0.5/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-1.1/100*[1 1],'b:','LineWidth',2);
    fprintf('Them %0.1f%%, us %0.1f%%\n\n',1.1,100*max(fractionBleached));
%}

% Pandiyan et al., 2020
%{
% 
% Their stimulus energy and bleaching fractions.  We assume this is light
% entering the eye.  
%     1.2% = 0.09x10^6 photons/um^2
%     3.9% = 0.30x10^6 photons/um^2
%     9.5% = 0.76x10^6 photons/um^2
%     23.9% = 2.14x10^6 photons/um^2
%
% Their stimulus exposures range from 3-70 ms across all their bleaching
% but given we have energy we don't need to worry about actual stimulus
% duration very much.
%
% We don't match their published values.  Correspondence with the authors%
% indicates they are doing a different bleaching calculation, which we will
% need to work through.
%
% 528 ± 20 nm

    % Initialize
    close all; clear;
    trolandFigure = figure; clf;
    set(gcf,'Position',[200 420 1200 1200]);

    % Size in degrees of one micron
    micronsPerDeg = 290;
    linearSizeDegs = 1/micronsPerDeg;

    % Timing
    nVideos = 2;
    nSets = 1;
    timePerVideoSec = 30;
    timePerStimSec = 50e-3;

    % Stimulus wavelength info
    stimulusWl = 528; 
    stimulusFWHM = 20; 
    whichStimCondition = 'stimulus';

    % 1.2% bleach
    figure(trolandFigure);
    subplot(2,2,1); hold on
    stimulusCornealPowerNw = 1e9*QuantaToEnergy(stimulusWl,0.09e6)/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.2g nW/deg2, %0.4f sec, %0.2g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-1.2/100*[1 1],'b:','LineWidth',2);

    % 3.9% bleach
    figure(trolandFigure);
    subplot(2,2,2); hold on
    stimulusCornealPowerNw = 1e9*QuantaToEnergy(stimulusWl,0.30e6)/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.2g nW/deg2, %0.4f sec, %0.2g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-3.9/100*[1 1],'b:','LineWidth',2);

    % 9.5% bleach
    figure(trolandFigure);
    subplot(2,2,3); hold on
    stimulusCornealPowerNw = 1e9*QuantaToEnergy(stimulusWl,0.76e6)/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.2g nW/deg2, %0.4f sec, %0.2g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-9.5/100*[1 1],'b:','LineWidth',2);

    % 23.9% bleach
    figure(trolandFigure);
    subplot(2,2,4); hold on
    stimulusCornealPowerNw = 1e9*QuantaToEnergy(stimulusWl,2.14e6)/timePerStimSec;
    bgWl = 785;
    bgCornealPowerUw = 0; 
    stimStartTimeSec = 1;
    title(sprintf('Stim: %d nm, %0.2g nW/deg2, %0.4f sec, %0.2g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    t_AObleachingKinetics_Multi_Trial_Bg;
    plot([0 totalTimeSec*nSets],1-23.9/100*[1 1],'b:','LineWidth',2);
%}

% Cone class specific bleaching calculations. 
% This one is basically a check.
%
% The plots generated here show that 545 light
% and a field location of 10 deg, cone based 
% bleaching calculations for L and M cones match up with
% troland based calculations. This is expected, because
% the half-bleach constant for isomerizations was set
% from estimates of isomerizations at 10 deg, for a 560
% nm light.  
%{
    % Initialize
    close all; clear;
    trolandFigure = figure; clf;
    set(gcf,'Position',[200 420 1200 1200]);
    isomerizationFigure = figure; clf;
    set(gcf,'Position',[200 420 1200 1200]);

    % Linear size of stimulus
    linearSizeDegs = 1;
    
    % Field size for computing things like receptor aperture.
    % This is more or less field location.  Under the hood,
    % these calculations may round to 2 or 10 deg numbers
    % with the rounding break at 4 deg;
    fieldLocationDegs = 10;

    % Timing
    nVideos = 3;
    nSets = 1;
    timePerVideoSec = 10;

    % Stimulus wavelength info
    whichStimCondition = 'stimulus';
    computeConeBleaching = true;

    % Test condition.
    stimulusWl = 545; 
    stimulusFWHM = 10; 
    stimulusCornealPowerNw = 153;
    bgWl = 785;
    bgCornealPowerUw = 90; 
    stimStartTimeSec = 1;
    timePerStimSec = 1;

    % Compute
    t_AObleachingKinetics_Multi_Trial_Bg;

    % Figure title
    figure(trolandFigure); hold on
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    figure(isomerizationFigure); hold on
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
%}

% S cone bleaching calculations. 
%
% Using eccentric location calculation for comparison.
% We are a bit into made-up land with the S cones.
%{
    % Initialize
    close all; clear;
    trolandFigure = figure; clf;
    set(gcf,'Position',[200 420 1200 1200]);
    isomerizationFigure = figure; clf;
    set(gcf,'Position',[200 420 1200 1200]);

    % Linear size of stimulus
    linearSizeDegs = 1;
    
    % Field size for computing things like receptor aperture.
    % This is more or less field location.  Under the hood,
    % these calculations may round to 2 or 10 deg numbers
    % with the rounding break at 4 deg;
    fieldLocationDegs = 10;

    % Timing
    nVideos = 10;
    nSets = 5;
    timePerVideoSec = 30;

    % Stimulus wavelength info
    whichStimCondition = 'stimulus';
    computeConeBleaching = true;

    % Test condition. This is about a 3 bleach for luminance
    stimulusWl = 480; 
    stimulusCornealPowerNw = 100;
    stimulusFWHM = 30; 
    timePerStimSec = 1;
    % stimulusWl = 460; 
    % stimulusCornealPowerNw = 8.3;

    bgWl = 785;
    bgCornealPowerUw = 90; 
    stimStartTimeSec = 1;

    % Label what we're doing
    fprintf('Wavelength: %d; Corneal Power nW: %0.1f; Duration %0.1f sec\n',stimulusWl,stimulusCornealPowerNw,timePerStimSec);

    % Compute
    t_AObleachingKinetics_Multi_Trial_Bg;

    % Figure title
    figure(trolandFigure); hold on
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 
    figure(isomerizationFigure); hold on
    title(sprintf('Stim: %d nm, %0.4g nW/deg2, %0.4f sec, %0.4g nJ/deg2',stimulusWl,stimulusCornealPowerNw/(linearSizeDegs^2),timePerStimSec,stimulusCornealPowerNw*timePerStimSec/(linearSizeDegs^2))); 

    % Report out fraction bleached
    fprintf('Max fraction bleached\n');
    fprintf('\tL: %0.2f%%\n',100*max(fractionBleachedCones(1,:)));
    fprintf('\tM: %0.2f%%\n',100*max(fractionBleachedCones(2,:)));
    fprintf('\tS: %0.2f%%\n',100*max(fractionBleachedCones(3,:)));        
%}

%% Do cone isomerizations too?
if (~exist("computeConeBleaching","var"))
    computeConeBleaching = false;
end

%% What to compute for
if (~exist("whichStimCondition","var"))
    whichStimCondition = 'both';
end

%% Create figures?
if (~exist("trolandFigure","var"))
    trolandFigure = figure;
end
if (computeConeBleaching & ~exist("isomerizationFigure","var"))
    isomerizationFigure = figure;
end

%% Stimulus spectral properties
%
% Keep stimulus and background separate from each other
% These are set to default values if not set up for
% a particular purpose before this part of the code is
% reached.
%
% Wavelength of narrowband stimulus field. This can be the imaging light
if (~exist("stimulusWl","var"))
    stimulusWl = 545;
end

% Stimulus power passing through pupil in nanowatts
if (~exist("stimulusCornealPowerNw","var"))
    stimulusCornealPowerNw = 153;
end

% Full-width at half max of stimlus (in nm).
if (~exist("stimulusFWHM","var"))
    stimulusFWHM = 5;
end

% Wavelength of narrowband background field. This can be the imaging light
if (~exist("stimulusWlBg ","var"))
    bgWl = 785;
end

% Background power passing through pupil in microwatts
if (~exist("bgCornealPowerUw","var"))
    bgCornealPowerUw = 90;
end

% Full-width at half max of background (in nm).
if (~exist("bgFWHM","var"))
    bgFWHM = 5;
end

%% Stimulus spatial parameters
%
% Linear side of square field in degs.  Assumed to be the
% same size for stimulus and background
if (~exist("linearSizeDegs","var"))
    linearSizeDegs = 1;
end

%% Stimulus spatial parameters
%
% Linear side of square field in degs.  Assumed to be the
% same size for stimulus and background
if (~exist("fieldLocationDegs","var"))
    fieldLocationDegs = 1;
end

%% Temporal parameters for stimulus
%
% Background is always on except dark adaptation
%
% How long each video lasts
if (~exist("timePerVideoSec","var"))
    timePerVideoSec = 6;
end

% When does stimulus turn on in each video
if (~exist("stimStartTimeSec","var"))
    stimStartTimeSec = 1;
end

% How long is stimlus turned on for?
if (~exist("timePerStimSec","var"))
    timePerStimSec = 1;
end

% How many videos in a row?
if (~exist("nVideos","var"))
    nVideos = 10;
end

% Dark adaptation time at the start of each video
if (~exist("darkAdaptationTimeSecs","var"))
    darkAdaptationTimeSecs = 120;
end

% How many sets of nVideos vidoes happen, with dark adaptation
% at the start of each one
if (~exist("nSets","var"))
    nSets = 5;
end

%% Make relative spectral power distributions.
%
% Approximated by a Gaussian with specified center wavlength and FWHM,
% and with total power given by the corneal power specified above.
% The call to trapz takes wavelength spacing into account when normalizing the power.

% Set up stimulus
%
% Wavelength support
deltaWl = 1;
wls = (400:deltaWl:900)';
nWls = length(wls);

% Stimulus spd on wavelength support
stimulusCornealPowerUw = stimulusCornealPowerNw*(1e-3);
stimulusRelSpd = normpdf(wls,stimulusWl,FWHMToStd(stimulusFWHM));
stimulusUnitSpd = stimulusRelSpd/trapz(wls,stimulusRelSpd);

% Set up background
bgRelSpd = normpdf(wls,bgWl,FWHMToStd(bgFWHM));
bgUnitSpd = bgRelSpd/trapz(wls,bgRelSpd);

%% Pupil size and eye length
pupilDiameterMm = 7;                        % Pupil diameter.
pupilAreaMm = pi*(pupilDiameterMm/2)^2;     % For convenience below, compute pupil area.
eyeLengthMm = 17;

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
stimulusSpdCornealPowerUw = stimulusCornealPowerUw*stimulusUnitSpd;
stimulusSpdRadiance = AOMonochromaticCornealPowerToRadiance(wls,wls,stimulusSpdCornealPowerUw*deltaWl,pupilDiameterMm,linearSizeDegs^2);

bgSpdCornealPowerUw = bgCornealPowerUw*bgUnitSpd;
bgStimulusSpdRadiance = AOMonochromaticCornealPowerToRadiance(wls,wls,bgSpdCornealPowerUw*deltaWl,pupilDiameterMm,linearSizeDegs^2);

% Make sure our computed radiance yields the desired corneal
% irradiance when we go in the other direction.  The magic
% numbers (1e6) and (1e-3) in the call just below do unit conversions
% from units we're using here to those expected by RadianceAndDegrees2ToCornIrradiance
stimulusSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(stimulusSpdRadiance,linearSizeDegs^2)*(1e6)*((1e-3)^2);
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

bgLuminance = trapz(wls,T_y'.*bgStimulusSpdRadiance);
fprintf('Background luminance is %0.1f cd/m2\n',bgLuminance);

%% Convert luminance to trolands for bleaching calcs
stimulusTrolands = stimulusLuminance*pupilAreaMm;
fprintf('Stimulus is %0.1f trolands\n',stimulusTrolands);

bgTrolands = bgLuminance*pupilAreaMm;
fprintf('Background is %0.1f trolands\n',bgTrolands);

%% Get cone isomerization rates using a function we already have
%
% These come back as a function of wavelength and the underlying routine
% works at 1 nm spacing so we can just add them up.
if (computeConeBleaching)
    conversionVerboseFlag = false;
    [tempAllLumSpd, tempAllPhotTdSpd, stimulusIsomerizationsSecSpd, SFunc] = ...
        AOLightLevelConversions_Func(linearSizeDegs, wls, stimulusSpdCornealPowerUw, conversionVerboseFlag,pupilDiameterMm,eyeLengthMm,fieldLocationDegs);
    [~, ~, bgIsomerizationsSecSpd, ~] = ...
        AOLightLevelConversions_Func(linearSizeDegs, wls, bgSpdCornealPowerUw, conversionVerboseFlag,pupilDiameterMm,eyeLengthMm,fieldLocationDegs);
    if (SFunc(2) ~= 1)
        error('Utility function not working at 1 nm spacing');
    end
    allLumFuncCheck = sum(tempAllLumSpd);
    if (abs(allLumFuncCheck - stimulusLuminance)/stimulusLuminance > 1e-4)
        error('Luminance calculation mismatch in two places');
    end
    allPhotTdFuncCheck = sum(tempAllPhotTdSpd);
    if (abs(allPhotTdFuncCheck - stimulusTrolands)/stimulusTrolands > 1e-4)
        error('Photopic td calculation mismatch in two places');
    end
    stimulusIsomerizationsSec = sum(stimulusIsomerizationsSecSpd,2);
    bgIsomerizationsSec = sum(bgIsomerizationsSecSpd,2);

    fprintf('Cone isomerizations/sec, field location degs: %0.1f\n',fieldLocationDegs);
    fprintf('\tL: %0.2g\n',stimulusIsomerizationsSec(1));
    fprintf('\tM: %0.2g\n',stimulusIsomerizationsSec(2));
    fprintf('\tS: %0.2g\n',stimulusIsomerizationsSec(3));
else
    stimulusIsomerizationsSec = NaN*ones(3,1);
    bgIsomerizationsSec = NaN*ones(3,1);
end

%% Compute a temporal stimulus timecourse
%
% Compute some temporal paramaters
totalTimeSec = timePerVideoSec*nVideos+(nVideos-1)*darkAdaptationTimeSecs;
totalTimeMsec = 1000*totalTimeSec;
timeMsec = (1:totalTimeMsec) - 1;
timeSec = timeMsec/1000;

switch (whichStimCondition)
    case 'stimulus'
        trolandsTime = zeros(size(timeMsec));
        isomerizationsSecTime = zeros(3,length(timeMsec));
        for ii = 1:nVideos
            startTime = (ii-1)*timePerVideoSec*1000 + stimStartTimeSec*1000;
            finishTime = startTime + timePerStimSec*1000;
            index = find(timeMsec > startTime & timeMsec < finishTime);
            trolandsTime(index) = stimulusTrolands;
            isomerizationsSecTime(:,index) = stimulusIsomerizationsSec*ones(1,length(index));
        end
        plotColor = 'g';

    case 'background'
        trolandsTime = bgTrolands*ones(size(timeMsec));
        isomerizationsSecTime = bgIsomerizationsSec*ones(1,length(timeMsec));
        plotColor = 'r';

    case 'both'
        trolandsTime = bgTrolands*ones(size(timeMsec));
        isomerizationsSecTime = bgIsomerizationsSec*ones(1,length(timeMsec));
        for ii = 1:nVideos
            startTime = (ii-1)*timePerVideoSec*1000 + stimStartTimeSec*1000;
            finishTime = startTime + timePerStimSec*1000;
            index = find(timeMsec > startTime & timeMsec < finishTime);
            trolandsTime(index) = trolandsTime(index) + stimulusTrolands;
            isomerizationsSecTime(:,index) = isomerizationsSecTime(:,index) + stimulusIsomerizationsSec*ones(1,length(index));
        end
        plotColor = 'b';

    otherwise
        error('Unknown stimulus condition specified');
end

%% Compute bleaching over time
for jj = 1:nSets
    if jj == 1
        % Troland calc
        initialFractionBleached = 0;
        fractionBleached = ComputePhotopigmentBleaching(trolandsTime,'cones','trolands','Boynton',initialFractionBleached,'msec');
        fractionUnbleached = 1 - fractionBleached;
        recovery = fractionBleached(1,totalTimeMsec);
        figure(trolandFigure);
        plot(timeSec,fractionUnbleached,plotColor,'LineWidth',2)
        fractionBleachedSets = zeros(nSets, length(fractionBleached));
        fractionUnbleachedSets = zeros(nSets, length(fractionUnbleached));
        fractionBleachedSets(jj,:) =  fractionBleached;
        fractionUnbleachedSets(jj,:) =  fractionUnbleached;

        % Isomerization calc
        if (computeConeBleaching)
            initialFractionBleachedCones = zeros(3,1);
            for rr = 1:3
                fractionBleachedCones(rr,:) = ComputePhotopigmentBleaching(isomerizationsSecTime(rr,:),'cones','isomerizations','Boynton',initialFractionBleachedCones(rr),'msec');
            end
            fractionUnbleachedCones = 1 - fractionBleachedCones;
            recoveryCones = fractionBleachedCones(:,totalTimeMsec);
            figure(isomerizationFigure); hold on
            plot(timeSec,fractionUnbleachedCones(1,:),'r','LineWidth',2)
            plot(timeSec,fractionUnbleachedCones(2,:),'g','LineWidth',2)
            plot(timeSec,fractionUnbleachedCones(3,:),'b','LineWidth',2)

            fractionBleachedSetsCones = zeros(nSets,3,length(fractionBleached));
            fractionUnbleachedSetsCones = zeros(nSets,3,length(fractionUnbleached));
            fractionBleachedSetsCones(jj,:,:) =  fractionBleachedCones;
            fractionUnbleachedSetsCones(jj,:,:) =  fractionUnbleachedCones;
        end

    else
        % Troland calc
        initialFractionBleached = recovery;
        fractionBleached = ComputePhotopigmentBleaching(trolandsTime,'cones','trolands','Boynton',initialFractionBleached,'msec');
        fractionUnbleached = 1 - fractionBleached;
        recovery = fractionBleached(1,totalTimeMsec);
        figure(trolandFigure);
        plot(timeSec+(totalTimeSec*(jj-1)),fractionUnbleached,plotColor,'LineWidth',2)
        fractionBleachedSets(jj,:) =  fractionBleached;
        fractionUnbleachedSets(jj,:) =  fractionUnbleached;

        % Isomerization calc
        initialFractionBleachedCones = recoveryCones;

        for rr = 1:3
            fractionBleachedCones(rr,:) = ComputePhotopigmentBleaching(isomerizationsSecTime(rr,:),'cones','isomerizations','Boynton',initialFractionBleachedCones(rr),'msec');
        end
        fractionUnbleachedCones = 1 - fractionBleachedCones;
        recoveryCones = fractionBleachedCones(:,totalTimeMsec);
        figure(isomerizationFigure); hold on
        plot(timeSec+(totalTimeSec*(jj-1)),fractionUnbleachedCones(1,:),'r','LineWidth',2)
        plot(timeSec+(totalTimeSec*(jj-1)),fractionUnbleachedCones(2,:),'g','LineWidth',2)
        plot(timeSec+(totalTimeSec*(jj-1)),fractionUnbleachedCones(3,:),'b','LineWidth',2)

        fractionBleachedSetsCones = zeros(nSets,3,length(fractionBleached));
        fractionUnbleachedSetsCones = zeros(nSets,3,length(fractionUnbleached));
        fractionBleachedSetsCones(jj,:,:) =  fractionBleachedCones;
        fractionUnbleachedSetsCones(jj,:,:) =  fractionUnbleachedCones;
    end
end

% Finish up plots
figure(trolandFigure);
xlim([0 totalTimeSec*nSets])
ylim([0 1]);
xlabel('Time (sec)');
ylabel('Fraction L/M Cone Pigment Unbleached (from trolands)');

figure(isomerizationFigure);
xlim([0 totalTimeSec*nSets])
ylim([0 1]);
xlabel('Time (sec)');
ylabel('Fraction L, M, and S Cone Pigment Unbleached');
legend({'L', 'M', 'S'},'Location','SouthEast');







