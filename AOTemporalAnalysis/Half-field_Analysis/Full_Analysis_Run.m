% Created by Robert F Cooper 12-19-2017
%


clear;
close all force;
clc;

rootDir = uigetdir(pwd);

%% Multi-trial run temporal analysis
fPaths = read_folder_contents_rec(rootDir,'tif');

wbh = waitbar(0,['Converting image 0 of ' num2str(length(fPaths)) '.']);

for i=1:size(fPaths,1)
        
    [mov_path, ref_image_fname] = getparent(fPaths{i});
    
    waitbar(i/length(fPaths), wbh, ['Analyzing image ' ref_image_fname ' (' num2str(i) ' of ' num2str(length(fPaths)) ').']);
        
    if reprocess || ~exist(fullfile(mov_path,'Mat_Profile_Data',[ref_image_fname(1:end - length('AVG.tif') ) 'box_cutoff_regional_norm_prestimminusdiv_sub_90_profiledata.mat'] ), 'file' );
        
        try    
            Temporal_Reflectivity_Analysis(mov_path,ref_image_fname);
        catch ex
           disp([ref_image_fname ' failed to process:']);
           disp([ex.message ': line ' num2str(ex.stack(1).line)] );
        end
    end
    close all;
end

% Output: mat files for each acquisition.

%% Aggregate multi-trial-run BOOTSTRAP

fPaths = read_folder_contents_rec(rootDir,'mat');

for i=1:size(fPaths,1)        
    [dataPath{i}, ref_image_fname] = getparent(fPaths{i});    
end

dataPath = unique(dataPath)';

waitbar(0,wbh,['Aggregating dataset 0 of ' num2str(length(dataPath)) '.']);

for i=1:length(dataPath)
       
    waitbar(i/length(dataPath), wbh, ['Aggregating dataset (' num2str(i) ' of ' num2str(length(dataPath)) ').']);

    try 
        % Find the control match to the stim
        [~, tail] = getparent(dataPath{i},2,'full');
        [parent, ~] = getparent(dataPath{i},3,'full');
        
        controlpath = fullfile(parent,'control',tail);
        
        if any( cellfun(@(s)strcmp(controlpath,s),dataPath) ) && ...
           ~strcmp(controlpath,dataPath{i})
           
            fitData(i) = Aggregate_Analyses_TwoSource_bootstrap(dataPath{i},controlpath);
            [remain kid] = getparent(dataPath{i});
            [remain region] = getparent(remain);
            [remain stim_time] = getparent(remain);
            [remain stim_intensity] = getparent(remain);
            [remain stimwave] = getparent(remain);
            [~, id] = getparent(remain);

            outFname = [id '_' stimwave '_' stim_intensity '_' stim_time '_all_amps'];

            all_amps = fitData(i).all_amps;
            figure(1); hist(all_amps,20); title(outFname,'Interpreter','none');
            saveas(gcf,[outFname '.png']);
            save([outFname '.mat'],'all_amps');
        elseif ~strcmp(controlpath,dataPath{i})
            warning(['Unable to find paired control video for,' parent]);
            fitData(i) = Aggregate_Multiple_Temporal_Analyses_bootstrap(dataPath{i});
            
            [remain kid] = getparent(dataPath{i});
            [remain region] = getparent(remain);
            [remain stim_time] = getparent(remain);
            [remain stim_intensity] = getparent(remain);
            [remain stimwave] = getparent(remain);
            [~, id] = getparent(remain);

            outFname = [id '_' stimwave '_' stim_intensity '_' stim_time '_all_amps'];

            all_amps = fitData(i).all_amps;
            figure(1); hist(all_amps,20); title(outFname,'Interpreter','none');
            saveas(gcf,[outFname '.png']);
            save([outFname '.mat'],'all_amps');
        else
            warning(['Not processing control video:' parent]);
        end   
    catch ex
       disp([dataPath{i} ' failed to analyze:']);
       disp([ex.message ' ' ex.stack(1).name ': line ' num2str(ex.stack(1).line)] );
    end
    
end

fitData = struct2table(fitData);

writetable(fitData,['Bootstrapped_Aggregate_Summary_' date '.csv']);

