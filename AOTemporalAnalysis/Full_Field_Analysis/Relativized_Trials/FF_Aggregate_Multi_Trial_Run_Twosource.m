
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

respData=[];
for i=1:length(dataPath)
       
    waitbar(i/length(dataPath), wbh, ['Aggregating dataset (' num2str(i) ' of ' num2str(length(dataPath)) ').']);

    try 
        % Find the control match to the stim
        [~, tail] = getparent(dataPath{i},2,'full');
        [parent, ~] = getparent(dataPath{i},4,'full');
        
        controlpath = fullfile(parent,'control',tail);
        
        if any( cellfun(@(s)strcmp(controlpath,s),dataPath) ) && ...
           ~strcmp(controlpath,dataPath{i})
           
            respData = [respData; FF_Aggregate_Analyses_TwoSource(dataPath{i},controlpath)];
        elseif ~strcmp(controlpath,dataPath{i})
            warning(['Unable to find paired control video for,' parent]);            
        else
            warning(['Not processing control video:' parent]);
        end
    catch ex
       warning([dataPath{i} ' failed to analyze:']);
       disp([ex.message ' ' ex.stack(1).name ': line ' num2str(ex.stack(1).line)] );
    end
    
%     pause;
end

close(wbh);

respData = struct2table(respData);

writetable(respData,['Aggregate_Summary_' date '.csv']);
% save(['Residuals_' date '.mat'], 'residuals');