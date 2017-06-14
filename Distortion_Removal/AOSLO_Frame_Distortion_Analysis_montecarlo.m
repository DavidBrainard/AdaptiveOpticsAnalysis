

clear;
clc

motion_path = pwd;

fNames = read_folder_contents(motion_path,'csv');

% Load in the motion files
for i=1:length(fNames)        
    AOSLO_Frame_Distortion_Analysis_function(motion_path, fNames{i});
end
