
clear;
close all force;
clc;



reprocess=true;

rootDir = uigetdir(pwd);

fPaths = read_folder_contents_rec(rootDir,'tif');

wbh = waitbar(0,['Processing trial 0 of ' num2str(length(fPaths)) '.']);

p = gcp();
for i=1:size(fPaths,1)
    
    [mov_path, ref_image_fname] = getparent(fPaths{i});

    waitbar(i/length(fPaths), wbh, ['Submitted trial (' num2str(i) ' of ' num2str(length(fPaths)) ').']); 
    vid_type= 'stimulus';
    if strcmpi(getparent(mov_path,2,'short'), 'control')
        vid_type = 'control';
    end

    try
        f(i) = parfeval(@Rel_FF_Temporal_Reflectivity_Analysis, 0, mov_path, ref_image_fname,[72 90],vid_type);
%             Rel_FF_Temporal_Reflectivity_Analysis(mov_path,ref_image_fname,[72 90],vid_type);
    catch ex
       disp([ref_image_fname ' failed to process:']);
       disp([ex.message ': line ' num2str(ex.stack(1).line)] );
    end

end

waitbar(0, wbh, 'Waiting for trials to finish...'); 

results = cell(1,size(fPaths,1));
for i = 1:size(fPaths,1)
    
    try
        % fetchNext blocks until next results are available.
        [completedIdx] = fetchNext(f);
        [mov_path, ref_image_fname] = getparent(fPaths{completedIdx});
        waitbar(i/length(fPaths), wbh, ['Finished trial (' num2str(i) ' of ' num2str(length(fPaths)) ').']);    
    catch ex
       disp([ref_image_fname ' failed to process:']);
       disp([ex.message ': line ' num2str(ex.stack(1).line)] );
    end
end
delete(p);
close(wbh);

