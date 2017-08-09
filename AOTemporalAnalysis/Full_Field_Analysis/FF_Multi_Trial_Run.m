
clear;
close all force;
clc;

reprocess=true;

rootDir = uigetdir(pwd);

fPaths = read_folder_contents_rec(rootDir,'tif');
wbh = waitbar(0,['Converting image 0 of ' num2str(length(fPaths)) '.']);

for i=1:size(fPaths,1)
    
    [mov_path, ref_image_fname] = getparent(fPaths{i});
    waitbar(i/length(fPaths), wbh, ['Finished image ' ref_image_fname ' (' num2str(i) ' of ' num2str(length(fPaths)) ').']);
    

    if reprocess || ~exist(fullfile(mov_path,'Profile_Data',[ref_image_fname(1:end - length('AVG.tif') ) 'box_cutoff_regional_norm_prestimminusdiv_sub_90_profiledata.mat'] ), 'file' );
        
        vid_type= 'stimulus';
        if strcmpi(getparent(mov_path,0,'short'), 'control')
            vid_type = 'control';
        end
        
        try
%             f(i) = parfeval(@FF_Temporal_Reflectivity_Analysis, 0, mov_path, ref_image_fname,[67 99],vid_type);
            FF_Temporal_Reflectivity_Analysis(mov_path,ref_image_fname,[67 99],vid_type);
        catch ex
           disp([ref_image_fname ' failed to process:']);
           disp([ex.message ': line ' num2str(ex.stack(1).line)] );
        end
    end
    close all;
end

% for i=1:size(fPaths,1)
%     waitbar(i/length(fPaths), wbh, ['Finished image ' ref_image_fname ' (' num2str(i) ' of ' num2str(length(fPaths)) ').']);
%     
%     fetchNext(f);
%     
% end
close(wbh);