%   Aggregate_Recursive_Run
%
%   Recursively calculates the pooled variance across a set of pre-analyzed 
%   signals from stimulus and control trials, performs the subtraction
%   between their standard deviations, and performs a piecewise fit of the 
%   resulting subtraction.
%
%   This script assumes that one mat file has valid stimulus *and*
%   control data.
%
% Outputs:
%
%   A date-stamped file containing information about the fit, named:
%   "Aggregate_Summary_[date].csv"
%
% Created by Robert F Cooper 12-31-2015
%
% The analyses performed by this script are from Cooper et al. "Non-invasive 
% assessment of human cone photoreceptor function", and
% encompassed in Figures 4B/5A, Equation 3.


clear;
close all force;
clc;

rootDir = uigetdir(pwd);

fPaths = read_folder_contents_rec(rootDir,'mat');

for i=1:size(fPaths,1)        
    [dataPath{i}, ref_image_fname] = getparent(fPaths{i});
    
end

dataPath = unique(dataPath)';

wbh = waitbar(0,['Aggregating dataset 0 of ' num2str(length(dataPath)) '.']);

for i=1:size(dataPath,1)
       
    waitbar(i/length(dataPath), wbh, ['Aggregating dataset (' num2str(i) ' of ' num2str(length(dataPath)) ').']);

    try    
        [fitData(i) residuals{i}] = Aggregate_Multiple_Temporal_Analyses(dataPath{i});
    catch ex
       disp([dataPath{i} ' failed to analyze:']);
       disp([ex.message ' ' ex.stack(1).name ': line ' num2str(ex.stack(1).line)] );
    end
end

close(wbh);

fitData = struct2table(fitData);

writetable(fitData,['Aggregate_Summary_' date '.csv']);