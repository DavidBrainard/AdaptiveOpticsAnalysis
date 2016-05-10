% Robert F Cooper 10-05-2016
%
% This script enables a user to remove poor frames from a temporal confocal
% and visible series, while simultaneously updating the acceptable_frames
% file.

clear
clc

%% Filename determination and handling
[confocal_fname, pathname] = uigetfile('*.avi', 'Select the cropped confocal temporal video');

confind = strfind(confocal_fname,'confocal');

if isempty(confind)
   error('Could not find confocal in the filename. Needed for proper function of this script!'); 
end

visible_fname = strrep(confocal_fname, 'confocal', 'visible');

% Find where the filename should be cut off in the confocal videos, and
% determine our acceptable frame filename.
i=1;
[comb_str remain] = strtok(confocal_fname(confind:end), '_');
while ~isempty(remain)
    [tok remain] = strtok( remain, '_');

    if i==4
        confocal_fname_out = comb_str;
    elseif i==8
        acceptable_frame_fname = comb_str;
        break;
    end
    
    comb_str = [comb_str '_' tok];
    
    i=i+1;
end
% Create our expected acceptable frame filenames
acceptable_frame_fname = [confocal_fname(1:confind-1) acceptable_frame_fname '_acceptable_frames.csv'];
acceptable_frame_fname_out = [confocal_fname(1:confind-1) confocal_fname_out '_crop_affine_acceptable_frames.csv'];

% Create our visible output filename
visible_fname_out = [confocal_fname(1:confind-1) confocal_fname_out '_crop_affine.avi'];
visible_fname_out = strrep(visible_fname_out, 'confocal', 'visible');

% Create our confocal output filename - affine will need to be done outside
% MATLAB.
confocal_fname_out = [confocal_fname(1:confind-1) confocal_fname_out '_crop.avi'];

confocal_vidobj = VideoReader( fullfile(pathname, confocal_fname) );
visible_vidobj = VideoReader( fullfile(pathname, visible_fname) );

%% File loading
confocal_vid = cell(1, round(confocal_vidobj.Duration*confocal_vidobj.FrameRate));
frame_nums = cell(1, round(confocal_vidobj.Duration*confocal_vidobj.FrameRate));
visible_vid = cell(1, round(confocal_vidobj.Duration*confocal_vidobj.FrameRate));

i=1;
while hasFrame(confocal_vidobj)
    
    confocal_vid{i} = readFrame(confocal_vidobj);
    visible_vid{i} = readFrame(visible_vidobj);
    frame_nums{i} = ['Frame ' num2str(i)];
    i=i+1;
end

acc_frame_list = dlmread( fullfile(pathname, acceptable_frame_fname) );

if length(acc_frame_list) ~= length(confocal_vid)
   error('Acceptable frames and confocal video list lengths do not match!');
end

[update, selected, multihand] = linkedImageFrameSingle(confocal_vid, frame_nums);


