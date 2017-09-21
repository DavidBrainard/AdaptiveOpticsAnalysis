%%CorrectIrradiacneForLumEquivalence
%
% Description:
%    This script corrects light power for the finite bandwidth of our
%    nominally monochromatic stimuli.  
%
%    The particular correction is to find the factor that converts
%    stimulus irradiance to an equivalent monochromatic irradiance,
%    (emi) such that a purely monochromatic light of the emi will have the
%    same luminance as the actual narrowband stimulus.
%
%    The reason this is useful is a) to evaluate whether we want to apply
%    a correction factor to our data for comparison with the photopic
%    luminosity function and b) if so to compute what that correction
%    factor should be.
%
%    One could also consider applying a similar correction for the
%    spectral sensitivity of the radiometer we use to measure irradiance,
%    if that spectral sensitivity is known.

% 09/21/17  dhb  Wrote it, following conversation and code snippets from Rob Cooper.

%% Clear
clear; close;

%% Parameters
