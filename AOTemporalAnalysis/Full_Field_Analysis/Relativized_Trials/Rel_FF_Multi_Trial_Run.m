
clear;
close all force;
clc;


% poolobj = gcp();
reprocess=true;

rootDir = uigetdir(pwd);

fPaths = read_folder_contents_rec(rootDir,'tif');

wbh = waitbar(0,['Processing trial 0 of ' num2str(length(fPaths)) '.']);

parfor i=1:size(fPaths,1)
    i
%     tic;
    [mov_path, ref_image_fname] = getparent(fPaths{i});

    if reprocess || ~exist(fullfile(mov_path,'Profile_Data',[ref_image_fname(1:end - length('AVG.tif') ) 'box_cutoff_regional_norm_prestimminusdiv_sub_90_profiledata.mat'] ), 'file' );
        
        vid_type= 'stimulus';
        if strcmpi(getparent(mov_path,0,'short'), 'control')
            vid_type = 'control';
        end
        
        try
%             f(i) = parfeval(@FF_Temporal_Reflectivity_Analysis, 0, mov_path, ref_image_fname,[67 99],vid_type);
            Rel_FF_Temporal_Reflectivity_Analysis(mov_path,ref_image_fname,[72 90],vid_type);
        catch ex
           disp([ref_image_fname ' failed to process:']);
           disp([ex.message ': line ' num2str(ex.stack(1).line)] );
        end
    end
    
    % If we've maxed out our number of workers wait until they've all
    % returned.
%     if mod(i,poolobj.NumWorkers) == 0
%         for m=1:poolobj.NumWorkers
%             fetchNext(f);
%             toc;
%             waitbar(i/length(fPaths), wbh, ['Finished trial ' ref_image_fname ' (' num2str(i) ' of ' num2str(length(fPaths)) ').']);        
%         end
%         tic;
%     end
    
    close all;
end

close(wbh);

