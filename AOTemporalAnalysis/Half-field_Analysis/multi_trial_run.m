% Multi_Trial_Run
%
%   Recursively performs a temporal analysis of a dataset. This dataset
%   must be formatted in a specific way, all relative to the "reference",
%   or averaged, image from a single confocal trial.
%
%   For a given grayscale temporal video file [TEMPORAL_FILE].avi, 
%       There must be:
%       -A reference image [TEMPORAL_FILE]_AVG.tif
%       -A cone coordinate list [TEMPORAL_FILE]_AVG_coords.csv
%       -A frame list [TEMPORAL_FILE]_acceptable_frames.csv
%       -A video file (the same size as [TEMPORAL_FILE] containing the 
%        stimulus extent and duration. Stimulus extent is delineated by
%        above-noise values (20 A.U.).
%
%
% Inputs:
%       mov_path: The folder path of the files that will be analyzed.
%
%       ref_image_fname: The name of the averaged image (ending in _AVG.tif)                        
%
%
% Created by Robert F Cooper 11-3-2015
% 
% 
% The analyses performed in this script are from:
% Cooper RF, Tuten WS, Dubra A, Brainard BH, Morgan JIW. 
% "Non-invasive assessment of human cone photoreceptor function." 
% Biomed Opt Express 8(11): 5098-5112. and are
% encompassed in Figures 1-4A, Equations 1-2.


clear;
close all force;
clc;

reprocess=true;

rootDir = uigetdir(pwd);

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

close(wbh);