function [allLum, allPhotTd, isomerizationsSec, S] = ...
    AOLightLevelConversions_Func(stimulusLinearSizeDegs, wavelengths, um_power_to_eye, verbose, pupilDiamMm, eyeLengthMm,fieldSizeDeg)
%
% Take monochromatic retinal irradiance and convert it to many other
% equivalent units.  Each return variable comes back as a funciton of
% wavelengths

% History:
%    10/29/15  dhb  Wrote (again, since I lost the first version) from OLLightLevelCheck
%    08/13/23  dhb  Tune up

%% Set defaults for backwards compatibility
if (nargin < 5 | isempty(pupilDiamMm))
    pupilDiamMm = 7;  % 7mm typical
end
if (nargin < 6 | isempty(eyeLengthMm))
    eyeLengthMm = 17; % 17mm - lens to retina typical
end
if (nargin < 7 | isempty(fieldSizeDeg))
    fieldSizeDeg = 2;
end


%% Define the wavelength spacing that we will work with
S = [min(wavelengths) 1 max(wavelengths)];
if (S(2) ~= 1)
    error('We are assuming 1 nm spacing');
end

if length(wavelengths) ~= length(um_power_to_eye)
    error('List of wavelengths and power must be the same!');
end

% Utility calculations
rawStimulusSideDegIn = stimulusLinearSizeDegs;
stimulusAreaDegrees2 = rawStimulusSideDegIn^2;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
eyeLengthCm = eyeLengthMm*(10^-1);
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);

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

%% Loop over wavelengths.
for w=1:length(wavelengths)

    theWavelength = wavelengths(w);

    if (um_power_to_eye(w) ~= 0)
        rawPowerIntoEyeIn = um_power_to_eye(w);
        rawCornIrradianceMicrowattsPerCm2In = rawPowerIntoEyeIn./pupilAreaCm2;
        rawRadianceMicrowattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(rawCornIrradianceMicrowattsPerCm2In,stimulusAreaDegrees2);
        rawRetIrradianceMicrowattsPerCm2In = RadianceAndPupilAreaEyeLengthToRetIrradiance(rawRadianceMicrowattsPerCm2Sr,S,pupilAreaCm2,eyeLengthCm);

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
        % mpePupilDiamMm  = GetWithDefault('Enter ANSI 2007 MPE caclulations assumed pupil diameter in mm',mpePupilDiamMm );
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

        cornealIrradianceWattsPerM2 = cornealIrradianceMicrowattsPerCm2*(10^-6)*(10^4);
        cornealIrradianceWattsPerCm2 = (10.^-4)*cornealIrradianceWattsPerM2;
        cornealIrradianceQuantaPerCm2Sec = EnergyToQuanta(S,cornealIrradianceWattsPerCm2);


        %% Report on stimulus
        if verbose
            fprintf('\n\n');
            fprintf('  *** Analyzing light levels for wavelength: %0.1fnm and power: %0.1f uW/degree2: \n',wavelengths(w),um_power_to_eye(w));
            fprintf('  * Stimulus radiance %0.1f log10 watts/[m2-sr], %0.1f log10 watts/[cm2-sr]\n',log10(sum(radianceWattsPerM2Sr)),log10(sum(radianceWattsPerCm2Sr)));
            fprintf('  * Stimulus luminance %0.1f candelas/m2\n',photopicLuminanceCdM2);
            fprintf('    * For comparison, sunlight in Philly: %0.1f cd/m2\n',photopicLuminancePhillyBrightCdM2);
            fprintf('  * Stimulus %0.0f (check val %0.0f) scotopic trolands, %0.0f photopic trolands (check val %0.0f)\n',irradianceScotTrolands,irradianceScotTrolands_check,...
                irradiancePhotTrolands,irradiancePhotTrolands_check);
            fprintf('  * Stimulus %0.1f log10 scotopic trolands, %0.1f log10 photopic trolands\n',log10(irradianceScotTrolands),log10(irradiancePhotTrolands));
            fprintf('  * Stimulus retinal irradiance %0.1f log10 watts/cm2\n',log10(sum(retIrradianceWattsPerCm2)));
            fprintf('  * Stimulus retinal irradiance %0.2e quanta/[um2-sec]\n',sum(retIrradianceQuantaPerCm2Sec)/(1000*1000));
            fprintf('  * Stimulus retinal irradiance %0.1f log10 quanta/[deg2-sec]\n',log10(sum(retIrradianceQuantaPerDeg2Sec)));
            fprintf('  * Stimulus corneal irradiance %0.1f log10 watts/cm2\n',log10(sum(cornealIrradianceWattsPerCm2)));
            fprintf('  * Stimulus corneal irradiance %0.1f log10 quanta/[cm2-sec]\n',log10(sum(cornealIrradianceQuantaPerCm2Sec)));
        end

        allLum(w) = photopicLuminanceCdM2;
        allPhotTd(w) = irradiancePhotTrolands;

        %% Get receptor sensitivities and compute isomerization rates, fraction bleached
        %
        % We'll be fancy and use the SilentSubstitionToolbox, since we
        % see that as the future of this sort of calculation.
        observerAgeYears = 32;
        receptorObj = SSTReceptorHuman('obsAgeInYrs', observerAgeYears, 'fieldSizeDeg', fieldSizeDeg, 'obsPupilDiameter', pupilDiamMm, 'doPenumbralConesTrueFalse', false, 'verbosity', 'none');
        nReceptors = length(receptorObj.labels);

        % LMS cones should be 1 to 3. Trust but verify
        %
        % A key variable is inner segment diameter, about which there is reasonable
        % uncertainty outside the fovea.
        if (verbose)
            fprintf('\nComputing cone isomerizations and bleaching\n');
        end
        for rr = 1:3
            % Check cone types and set numbers.  This is set up to allow
            % eccentricity specific cone IS diameters, but the specific
            % numbers are pretty rough and ready.  May want to update SST
            % to carry these numbers around.
            switch (rr)
                case 1
                    if (~strcmp(receptorObj.labels{rr},'LCone'))
                        error('L cone is not where we expect it');
                    end
                    if (fieldSizeDeg <= 4)
                        ISdiameterUm(rr) = PhotoreceptorDimensions('FovealLCone','ISdiam','Human','Rodieck');
                    else
                        ISdiameterUm(rr) = PhotoreceptorDimensions('LCone','ISdiam','Human','Webvision');
                    end
                case 2
                    if (~strcmp(receptorObj.labels{rr},'MCone'))
                        error('M cone is not where we expect it');
                    end
                    if (fieldSizeDeg <= 4)
                        ISdiameterUm(rr) = PhotoreceptorDimensions('FovealMCone','ISdiam','Human','Rodieck');
                    else
                        ISdiameterUm(rr) = PhotoreceptorDimensions('MCone','ISdiam','Human','Webvision');
                    end
                case 3
                    if (~strcmp(receptorObj.labels{rr},'SCone'))
                        error('S cone is not where we expect it');
                    end
                    if (fieldSizeDeg <= 4)
                        ISdiameterUm(rr) = PhotoreceptorDimensions('FovealSCone','ISdiam','Human','Rodieck');
                    else
                        ISdiameterUm(rr) = PhotoreceptorDimensions('SCone','ISdiam','Human','Webvision');
                    end
            end

            % Get fundamentals in form that expects radiance in quantal units and
            % is in terms of probability of isomerization per incidenct quantum.
            T_quantalIsomerizations = SplineCmf(receptorObj.S,receptorObj.T.T_quantalIsomerizations(rr,:),S);

            % Need a cone collecting area.
            ISareaUm2 = pi*(ISdiameterUm(rr)/2)^2;

            % Get isom/sec and bleaching
            % 
            % Note that fraction bleached isn't additive over wavelength so these numbers
            % are not all that useful generally.  That is why we don't return them
            isomerizationsSec(rr,w) = ISareaUm2*T_quantalIsomerizations*retIrradianceQuantaPerUm2Sec;
            fractionBleachedFromIsomWl(rr) = ComputePhotopigmentBleaching(isomerizationsSec(rr),'cones','isomerizations','Boynton');
            if (verbose)
                fprintf('    * %s, inner segment diam %0.2f um, isom/sec %0.2f, fraction bleached %0.2f\n',receptorObj.labels{rr},ISdiameterUm(rr),isomerizationsSec(rr),fractionBleachedFromIsomWl(rr));
            end
        end

   % No power at this wavelength, just set zero at this wavelength in the
   % return variables.
    else
        allLum(w) = 0;
        allPhotTd(w) = 0;
        isomerizationsSec(1:3,w) = 0;
    end
end

