%   FF_aggregate_Multi_Trial_Run_Twosource
%
%   Recursively calculates the pooled variance across a set of pre-analyzed 
%   signals from stimulus and control trials, performs the subtraction
%   between their standard deviations, and performs a piecewise fit of the 
%   resulting subtraction.
%
%   If it cannot find a control video, then it simply runs the default
%   "Aggregate_Multiple_Temporal_Analyses", which assumes that one mat file
%   has valid stimulus and control data.
%
% Outputs:
%
%   A date-stamped file containing information about the fit, named:
%   "Aggregate_Summary_[date].csv"
%
% Created by Robert F Cooper 06-20-2017
%
% The analyses performed in this script are from:
% Cooper RF, Tuten WS, Dubra A, Brainard BH, Morgan JIW. 
% "Non-invasive assessment of human cone photoreceptor function." 
% Biomed Opt Express 8(11): 5098-5112. and are
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

fitData=[];
for i=1:length(dataPath)
       
    waitbar(i/length(dataPath), wbh, ['Aggregating dataset (' num2str(i) ' of ' num2str(length(dataPath)) ').']);

    try 
        % Find the control match to the stim
        [~, tail] = getparent(dataPath{i},1,'full');
        [parent, ~] = getparent(dataPath{i},3,'full');
        
        controlpath = fullfile(parent,'control',tail);
        
        if any( cellfun(@(s)strcmp(controlpath,s),dataPath) ) && ...
           ~strcmp(controlpath,dataPath{i})
           
            fitData = [fitData; FF_Aggregate_Analyses_TwoSource(dataPath{i},controlpath)];
        elseif ~strcmp(controlpath,dataPath{i})
            warning(['Unable to find paired control video for,' parent]);            
        else
            warning(['Not processing control video:' parent]);
        end
    catch ex
       disp([dataPath{i} ' failed to analyze:']);
       disp([ex.message ' ' ex.stack(1).name ': line ' num2str(ex.stack(1).line)] );
    end
    
end

close(wbh);

fitData = struct2table(fitData);

writetable(fitData,['Aggregate_Summary_' date '.csv']);
