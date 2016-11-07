
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

ids = [];
signals = [];

colors=['r','g','b','c','m','y','k'];
c=0;
hz=16.66666666;
timeBase = (1:249)/hz;

for i=1:size(dataPath,1)
       
    waitbar(i/length(dataPath), wbh, ['Aggregating dataset (' num2str(i) ' of ' num2str(length(dataPath)) ').']);

    try    
        [pooled_control{i}, pooled_stim{i}, allmax, subid] = pooled_stddev_from_dir(dataPath{i});
        
        if isempty(ids) || ~strcmp(subid, ids(end,:))
            ids = [ids; subid];
            c = c+1;
            
            signals = [signals; nan(1,249)];
            
            for j=1:length( pooled_control{i} )
                signals(c,j) = pooled_control{i}(j);
            end
                        
            pooledcolors = colors(c);
        else
            
            for j=1:length( pooled_control{i} )
                if isnan(signals(c,j))
                    signals(c,j) = pooled_control{i}(j);
                else
                    if ~isnan(pooled_control{i}(j))
                        signals(c,j) = (signals(c,j)+pooled_control{i}(j))/2;
                    end
                end
            end
        end
        
        figure(10); hold on;
        plot( (1:length(pooled_control{i}))/hz,pooled_control{i},pooledcolors);
        
        
        
%         figure(10); 
%         plot( timeBase,sqrt(signals(1,:)),colors(1));
%         for j=2:1:size(signals,1)
%             hold on;
%             plot( timeBase,sqrt(signals(j,:)),colors(j));
%         end
%         hold off;
%         legend(ids)
                

        % Stim train
%         stimlen = str2double( strrep(stim_time(1:3),'p','.') );

%         trainlocs = 68/hz:1/hz:(68/hz+stimlen);
%         plot(trainlocs, max(pooled_stim)*ones(size(trainlocs)),'r*'); hold off;
        
        ylabel('Pooled Standard deviation'); xlabel('Time (s)');
        hold off;
    catch ex
       disp([ref_image_fname ' failed to analyze:']);
       disp([ex.message ' ' ex.stack(1).name ': line ' num2str(ex.stack(1).line)] );
    end

%     pause;
end

close(wbh);
