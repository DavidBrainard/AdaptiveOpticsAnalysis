% FF_Multi_Trial_Run
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


% Constants:
%
% The frames over which the stimulus was delivered. It is extremely important
% that the first number be correct, as the software will use that as the last
% frame included in normalization.
STIM_RANGE = [67 99]; 

rootDir = uigetdir(pwd);

fPaths = read_folder_contents_rec(rootDir,'tif');

wbh = waitbar(0,['Processing trial 0 of ' num2str(length(fPaths)) '.']);
for i=1:size(fPaths,1)
    tic;
    [mov_path, ref_image_fname] = getparent(fPaths{i});


        vid_type= 'stimulus';
        if strcmpi(getparent(mov_path,0,'short'), 'control')
            vid_type = 'control';
        end
        
        try
            FF_Temporal_Reflectivity_Analysis(mov_path,ref_image_fname,STIM_RANGE,vid_type);
        catch ex
           disp([ref_image_fname ' failed to process:']);
           disp([ex.message ': line ' num2str(ex.stack(1).line)] );
        end

    

    toc;
    waitbar(i/length(fPaths), wbh, ['Finished trial ' ref_image_fname ' (' num2str(i) ' of ' num2str(length(fPaths)) ').']);        

    
    close all;
end

close(wbh);
