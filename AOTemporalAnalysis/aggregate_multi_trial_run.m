
clear;
close all force;
clc;

rootDir = uigetdir(pwd);

fPaths = read_folder_contents_rec(rootDir,'mat');

for i=1:size(fPaths,1)        
    [dataPath{i}, ref_image_fname] = getparent(fPaths{i});
end

dataPath = unique(dataPath)';

wbh = waitbar(0,['Aggregating dataset 0 of ' num2str(length(fPaths)) '.']);

for i=1:size(dataPath,1)
       
    waitbar(i/length(fPaths), wbh, ['Aggregating dataset 0 of ' ref_image_fname ' (' num2str(i) ' of ' num2str(length(fPaths)) ').']);

    try    
        fitData(i) = Aggregate_Multiple_Temporal_Analyses(dataPath);
    catch ex
       disp([ref_image_fname ' failed to analyze:']);
       disp([ex.message ': line ' num2str(ex.stack(1).line)] );
    end

    close all;
end

close(wbh);