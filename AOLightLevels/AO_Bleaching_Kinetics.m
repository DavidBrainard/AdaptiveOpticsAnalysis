% Robert F Cooper 2018-11-09
%
% This script is designed to calculate the bleaching kinetics of a variable
% number of stimuli.


I = 1; % Stimulus Intensity in Td
I_0 = 20000; %in Td %Stimulus intensity that bleaches at the rate of 1/N
            % 73.7 Td for rods, 20000 Td for cones

N = 120; % Scaling factor, where 400=rhodopsin, 120=L/M cones.
p = 1; % Percentage of bleached pigment that we start with


dp_dt = @(p, I) ((1-p)./N) - (( I.*p)./ (N.*I_0));

