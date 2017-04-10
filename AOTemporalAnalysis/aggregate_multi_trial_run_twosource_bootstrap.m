
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
        elseif ~strcmp(controlpath,dataPath{i})
            warning(['Unable to find paired control video for,' parent]);
            fitData(i) = Aggregate_Multiple_Temporal_Analyses_bootstrap(dataPath{i});
        else
            warning(['Not processing control video:' parent]);
        end
        

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
    catch ex
       disp([dataPath{i} ' failed to analyze:']);
       disp([ex.message ' ' ex.stack(1).name ': line ' num2str(ex.stack(1).line)] );
    end
    
%     pause;
end

close(wbh);

fitData = struct2table(fitData);

writetable(fitData,['Bootstrapped_Aggregate_Summary_' date '.csv']);
