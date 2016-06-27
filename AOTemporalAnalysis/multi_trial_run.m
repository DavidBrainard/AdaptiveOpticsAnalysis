
clear;
close all force;
clc;

[mov_path] = uigetdir(pwd);

[fnamelist, isdir ] = read_folder_contents(mov_path,'tif');

for i=1:size(fnamelist,1)
    
    ref_image_fname = fnamelist{i};
    
    Temporal_Reflectivity_Analysis(mov_path,ref_image_fname);
    
    
end