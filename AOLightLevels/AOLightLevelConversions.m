% AOLightLevelConversions
%
% Take monochromatic retinal irradiance and convert it to many other
% equivalent units.
%
% 10/29/15  dhb  Wrote (again, since I lost the first version) from OLLightLevelCheck

%% Clear and close
clear; close all

%% Define the wavelength spacing that we will work with
S = [380 1 600];
if (S(2) ~= 1)
    error('We are assuming 1 nm spacing');
end

%% Get inputs and put into common format
theWavelength = GetWithDefault('Enter wavelength in nanometers',670);
pupilDiamMm = GetWithDefault('Enter pupil diameter in mm',7);2.54*(3.6*10-2)
eyeLengthMm = GetWithDefault('Enter assumed eye length in mm',17);
LENGTHORAREA = GetWithDefault('Enter stimulus side (1), stimulus diameter (2), or stimulus area (3)?',2);
switch (LENGTHORAREA)
    case 1
        rawStimulusSideDegIn = GetWithDefault('Enter square stimulus side in deg',1);
        stimulusAreaDegrees2 = rawStimulusSideDegIn^2;
    case 2
        rawStimulusSideDegIn = GetWithDefault('Enter circular stimulus diameter in deg',1.8);
        stimulusAreaDegrees2 = pi*(rawStimulusSideDegIn/2)^2;
    case 3
        stimulusAreaDegrees2 = GetWithDefault('Enter stimulus area in square degrees',1);
end

% Utility calculations
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
eyeLengthCm = eyeLengthMm*(10^-1);
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);

% Accept either corneal or retinal illuminance
RETORCORN = GetWithDefault('Enter retinal irradiance (1), corneal irradiance (2), or power entering eye (3)?',3);
switch (RETORCORN)
    case 1
        rawRetIrradianceMicrowattsPerCm2In = GetWithDefault('Enter retinal illuminance in microwatts/cm2',25.5);
        rawCornIrradianceMicrowattsPerCm2In = rawRetIrradianceMicrowattsPerCm2In*stimulusAreaDegrees2/pupilAreaCm2;
    case 2
        rawCornIrradianceMicrowattsPerCm2In = GetWithDefault('Enter corneal illuminance in microwatts/cm2',25.5);
        rawRadianceMicrowattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(rawCornIrradianceMicrowattsPerCm2In,stimulusAreaDegrees2);
        rawRetIrradianceMicrowattsPerCm2In = RadianceAndPupilAreaEyeLengthToRetIrradiance(rawRadianceMicrowattsPerCm2Sr,S,pupilAreaCm2,eyeLengthCm);
    case 3
        rawPowerIntoEyeIn = GetWithDefault('Enter power entering eye in microwatts',6000);
        rawCornIrradianceMicrowattsPerCm2In = rawPowerIntoEyeIn/pupilAreaCm2;
        rawRadianceMicrowattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(rawCornIrradianceMicrowattsPerCm2In,stimulusAreaDegrees2);
        rawRetIrradianceMicrowattsPerCm2In = RadianceAndPupilAreaEyeLengthToRetIrradiance(rawRadianceMicrowattsPerCm2Sr,S,pupilAreaCm2,eyeLengthCm);  
end
 
%% Given that we have retinal irradiance, corneal irradiance, the pupil area, and the stimulus area ...
% it is easy to check the retinal irradiance.
retIrradianceMicrowattsPerDeg2Check = rawCornIrradianceMicrowattsPerCm2In*pupilAreaCm2/stimulusAreaDegrees2;
retIrradianceMicrowattsPerCm2Check = retIrradianceMicrowattsPerDeg2Check/(0.01*DegreesToRetinalMM(1,eyeLengthCm*10)^2);
if (abs(rawRetIrradianceMicrowattsPerCm2In-retIrradianceMicrowattsPerCm2Check)/retIrradianceMicrowattsPerCm2Check > 0.05)
    error('Do not get same retinal irradiance two different ways');
end

%% Turn it into a spectral function
wls = SToWls(S);
retIrradianceMicrowattsPerCm2In = zeros(size(wls));
index = find(wls == theWavelength);
if (length(index) ~= 1)
    error('Specified wavelength is not in sampled wavelengths');
end
retIrradianceMicrowattsPerCm2In(index) = rawRetIrradianceMicrowattsPerCm2In;

%% Convert to some other units for convenience
retIrradianceMilliwattsPerCm2In = retIrradianceMicrowattsPerCm2In*(10^-3);
retIrradianceWattsPerCm2In = retIrradianceMilliwattsPerCm2In*(10^-3);
retIrradianceWattsPerUm2In = retIrradianceWattsPerCm2In*(10^-8);

%% Load CIE functions.   
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);

%% Load cone spectral sensitivities
load T_cones_ss2
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);

%% Load in a directly measured sunlight through window
% and off piece of white paper towel on floor for comparison.
% Surely that is safe to look at.
load spd_phillybright
spd_phillybright = SplineSpd(S_phillybright,spd_phillybright,S);
photopicLuminancePhillyBrightCdM2 = T_xyz(2,:)*spd_phillybright;

%% Compute irradiance to equivalent radiance
%
% This conversion depends on pupil area and eye length

% Do the conversion
radianceMilliwattsPerCm2Sr = RetIrradianceAndPupilAreaEyeLengthToRadiance(retIrradianceMilliwattsPerCm2In,S,pupilAreaCm2,eyeLengthCm);
radianceMicrowattsPerCm2Sr = radianceMilliwattsPerCm2Sr*(10^3);
radianceWattsPerCm2Sr = radianceMilliwattsPerCm2Sr*(10^-3);
radianceWattsPerM2Sr = radianceWattsPerCm2Sr*(10^4);

% Check by converting back using a different routine
retIrradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,eyeLengthMm);
retIrradianceWattsPerM2 = retIrradianceWattsPerUm2*(10^12);
retIrradianceWattsPerCm2 = retIrradianceWattsPerUm2*(10.^8);
retIrradianceMilliwattsPerCm2 = retIrradianceWattsPerCm2*(10^3);
retIrradianceMicrowattsPerCm2 = retIrradianceMilliwattsPerCm2*(10^3); 
tolerance = 1e-8;
if (max(abs(retIrradianceMicrowattsPerCm2(:)-retIrradianceMicrowattsPerCm2In(:))) > tolerance)
    error('Failure to go there and back');
end

%% Conver retinal irradiance to trolands, etc.
irradianceScotTrolands = RetIrradianceToTrolands(retIrradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));
irradiancePhotTrolands = RetIrradianceToTrolands(retIrradianceWattsPerUm2, S, 'Photopic', [], num2str(eyeLengthMm));
retIrradianceQuantaPerUm2Sec = EnergyToQuanta(S,retIrradianceWattsPerUm2);
retIrradianceQuantaPerCm2Sec = (10.^8)*retIrradianceQuantaPerUm2Sec;
retIrradianceQuantaPerDeg2Sec = (degPerMm^2)*(10.^-2)*retIrradianceQuantaPerCm2Sec;

% Yet more units
photopicLuminanceCdM2 = T_xyz(2,:)*radianceWattsPerM2Sr;

%% Pupil adjustment factor for Ansi MPE 
mpePupilDiamMm = 3;
mpePupilDiamMm  = GetWithDefault('Enter ANSI 2007 MPE caclulations assumed pupil diameter in mm',mpePupilDiamMm );
pupilAdjustFactor = (pupilDiamMm/mpePupilDiamMm).^2;

%% Get trolands another way.  For scotopic trolands, this just uses scotopic vlambda (in PTB as T_rods)
% and the magic factor of 1700 scotopic lumens per Watt from Wyszecki & Stiles (2cd edition),
% p. 257.  (This is the analog of 683 photopic lumens per Watt.  Then apply formula from
% page 103 of same book.
%
% Same idea for photopic trolands, although we already have luminance in cd/m2 from above so
% we can short cut a little.
%
% The agreement is good to integer scotopic trolands and I'm will to write off the rest
% as round off error.
load T_rods
T_scotopicVlambda = SplineCmf(S_rods,T_rods,S);
irradianceScotTrolands_check = pupilAreaMm2*1700*(T_scotopicVlambda*radianceWattsPerM2Sr);
irradiancePhotTrolands_check = pupilAreaMm2*photopicLuminanceCdM2;

%% Compute corneal irradiance
cornealIrradianceMicrowattsPerCm2 = RadianceAndDegrees2ToCornIrradiance(radianceMicrowattsPerCm2Sr,stimulusAreaDegrees2);
switch (RETORCORN)
    case 2
        index = find(cornealIrradianceMicrowattsPerCm2 ~= 0);
        if (length(index) ~= 1)
            error('Oops.  Somebody isn''t monochromatic');
        end
        if (abs(cornealIrradianceMicrowattsPerCm2(index)-rawCornIrradianceMicrowattsPerCm2In) > tolerance)
            error('Corneal irradiance computation did not invert');
        end
end
cornealIrradianceWattsPerM2 = cornealIrradianceMicrowattsPerCm2*(10^-6)*(10^4);
cornealIrradianceWattsPerCm2 = (10.^-4)*cornealIrradianceWattsPerM2;
cornealIrradianceQuantaPerCm2Sec = EnergyToQuanta(S,cornealIrradianceWattsPerCm2);


%% Report on stimulus
fprintf('\n');
fprintf('  * Analyzing light levels as input\n');
fprintf('  * Stimulus radiance %0.1f log10 watts/[m2-sr], %0.1f log10 watts/[cm2-sr]\n',log10(sum(radianceWattsPerM2Sr)),log10(sum(radianceWattsPerCm2Sr)));
fprintf('  * Stimulus luminance %0.1f candelas/m2\n',photopicLuminanceCdM2);
fprintf('    * For comparison, sunlight in Philly: %0.1f cd/m2\n',photopicLuminancePhillyBrightCdM2);
fprintf('  * Stimulus %0.0f (check val %0.0f) scotopic trolands, %0.0f photopic trolands (check val %0.0f)\n',irradianceScotTrolands,irradianceScotTrolands_check,...
    irradiancePhotTrolands,irradiancePhotTrolands_check);
fprintf('  * Stimulus %0.1f log10 scotopic trolands, %0.1f log10 photopic trolands\n',log10(irradianceScotTrolands),log10(irradiancePhotTrolands));
fprintf('  * Stimulus retinal irradiance %0.1f log10 watts/cm2\n',log10(sum(retIrradianceWattsPerCm2)));
fprintf('  * Stimulus retinal irradiance %0.1f log10 quanta/[cm2-sec]\n',log10(sum(retIrradianceQuantaPerCm2Sec)));
fprintf('  * Stimulus retinal irradiance %0.1f log10 quanta/[deg2-sec]\n',log10(sum(retIrradianceQuantaPerDeg2Sec)));
fprintf('  * Stimulus corneal irradiance %0.1f log10 watts/cm2\n',log10(sum(cornealIrradianceWattsPerCm2)));
fprintf('  * Stimulus corneal irradiance %0.1f log10 quanta/[cm2-sec]\n',log10(sum(cornealIrradianceQuantaPerCm2Sec)));
   

%% Get MPE from as a function of wavelength.  For each wavelength,
% take minimum radiance over specified sizes and durations.

% Specify what parameters to test
minLogSize = -1; maxLogSize = 2;
minLogDuration = -1; maxLogDuration = 4;
minLogYRad = -3; maxLogYRad = 2;
minLogYIrrad = -5; maxLogYIrrad = 0;
minLogYIntRad = 0; maxLogYIntRad = 3;
minLogYRadExp = -4; maxLogYRadExp = -1;

% Adjust info to be monochromatic
S_MPE = [theWavelength 1 1];
measuredWls = SToWls(S_MPE);
index = find(measuredWls >= 400);
stimulusWavelengthsNm = measuredWls(index);
stimulusSizesDeg = logspace(minLogSize,maxLogSize,5);
stimulusDurationsSec = logspace(minLogDuration,maxLogDuration,5);
%fprintf('Computing MPE over wavelengths from %0.1f to %0.1f deg\n',min(stimulusWavelengthsNm),max(stimulusWavelengthsNm));
clear MPELimitIntegratedRadiance_JoulesPerCm2Sr MPELimitRadiance_WattsPerCm2Sr MPELimitCornealIrradiance_WattsPerCm2 MPELimitCornealRadiantExposure_JoulesPerCm2
for w = 1:length(stimulusWavelengthsNm)
    stimulusWavelengthNm = stimulusWavelengthsNm(w);
    if (rem(w,10) == 0)  
        %fprintf('\tComputing minimum MPE for wavelength %d nm\n',stimulusWavelengthNm);
    end
    MPELimitIntegratedRadiance_JoulesPerCm2Sr(w) = Inf;
    MPELimitRadiance_WattsPerCm2Sr(w) = Inf;
    MPELimitCornealIrradiance_WattsPerCm2(w) = Inf;
    MPELimitCornealRadiantExposure_JoulesPerCm2(w) = Inf;
    for s = 1:length(stimulusSizesDeg)
        stimulusSizeDeg = stimulusSizesDeg(s);
        stimulusSizeMrad = DegToMrad(stimulusSizeDeg);
        for t = 1:length(stimulusDurationsSec)
            stimulusDurationSec = stimulusDurationsSec(t);
            
            % Compute MPE.  We don't understand how the cone limit computations fit in with
            % the standard, or not.  So, we run it both ways and take the lower limit returned.
            [temp1, temp2, temp3, temp4] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,0);
            [temp5, temp6, temp7, temp8] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,1);
            if (temp5 < temp1)
                temp1 = temp5;
            end
            if (temp6 < temp2)
                temp2 = temp6;
            end
            if (temp7 < temp3);
                temp3 = temp7;
            end
            if (temp8 < temp4)
                temp4 = temp8;
            end
            clear temp5 temp6 temp7 temp8
            
            % Store minimum at each wavelength.
            if (temp1 < MPELimitIntegratedRadiance_JoulesPerCm2Sr(w))
                MPELimitIntegratedRadiance_JoulesPerCm2Sr(w) = temp1;
            end
            if (temp2 < MPELimitRadiance_WattsPerCm2Sr(w))
                MPELimitRadiance_WattsPerCm2Sr(w) = temp2;
            end
            if (temp3 < MPELimitCornealIrradiance_WattsPerCm2(w))
                MPELimitCornealIrradiance_WattsPerCm2(w) = temp3;
            end
            if (temp4 < MPELimitCornealRadiantExposure_JoulesPerCm2(w))
                MPELimitCornealRadiantExposure_JoulesPerCm2(w) = temp4;
            end
        end
    end
end

%% Find how much total radiance we could tolerate if all our power was at the 
% wavelength with minimum MPE.
minMPERadiance = min(MPELimitRadiance_WattsPerCm2Sr(:));
fprintf('\n');
fprintf('  * Compute ANSI 2007 MPE as a function of wavelength.  For each wavelength, took minimum over size and duration\n');
fprintf('    * Size range: %0.1f to %0.1f degrees\n',min(stimulusSizesDeg),max(stimulusSizesDeg));
fprintf('    * Duration range: %0.1f to %0.1f seconds\n',min(stimulusDurationsSec),max(stimulusDurationsSec));
fprintf('  * Minimum ANSI MPE value over wavelengths: radiance %0.1f log W/[cm2-sr]\n',log10(minMPERadiance));
fprintf('    * Compare with total stimulus radiance %0.1f log  W/[cm2-sr]\n',log10(sum(radianceWattsPerCm2Sr)));
fprintf('    * Compare with total pupil adjusted radiance %0.1f log  W/[cm2-sr]\n',log10(sum(radianceWattsPerCm2Sr))+log10(pupilAdjustFactor));
fprintf('    * Pupil adjustment assumes observer pupil diameter of %0.1f mm, MPE standard diameter of %0.1f mm\n',pupilDiamMm,mpePupilDiamMm);
fprintf('\n');

%% Sum over wavelength of power divided by MPE
% Could put this back in, but would have to think
% a bit harder about wavelength spacing adjustment.
%
% index = find(stimulusWavelengthsNm >= 400);
% deltaMeasuredWls = measuredWls(2)-measuredWls(1);
% deltaMPEWls = stimulusWavelengthsNm(2)-stimulusWavelengthsNm(1);
% MPERatioSum = 0;
% for i = 1:length(stimulusWavelengthsNm)
%     index = find(measuredWls == stimulusWavelengthsNm(i));
%     MPERatioSum = MPERatioSum + radianceWattsPerCm2Sr(index)/MPELimitRadiance_WattsPerCm2Sr(i);
% end
% fprintf('MPERatioSum = %0.4f\n',MPERatioSum*deltaMPEWls/deltaMeasuredWls);

%% Now compare to the ISO Standard
stimulusDurationForISOMPESecs = 60*60;
[IsOverLimit,ISO2007MPEStruct] = ISO2007MPECheckType1ContinuousRadiance(S,radianceWattsPerM2Sr,stimulusDurationForISOMPESecs,stimulusAreaDegrees2,eyeLengthMm);
fprintf('  * ISO MPE Analysis\n');
ISO2007MPEPrintAnalysis(IsOverLimit,ISO2007MPEStruct);
fprintf('  * Assumed duration seconds %0.1f, hours %0.1f\n',stimulusDurationForISOMPESecs,stimulusDurationForISOMPESecs/3600);

