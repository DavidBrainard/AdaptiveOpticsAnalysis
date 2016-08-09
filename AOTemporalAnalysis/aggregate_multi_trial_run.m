
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
        fitData(i) = Aggregate_Multiple_Temporal_Analyses(dataPath{i});
    catch ex
       disp([ref_image_fname ' failed to analyze:']);
       disp([ex.message ' ' ex.stack(1).name ': line ' num2str(ex.stack(1).line)] );
    end

    
end

close(wbh);

fitData = struct2table(fitData);

writetable(fitData,['Aggregate_Summary_' date '.csv']);