% Robert F Cooper 10-05-2016
%
% This script enables a user to remove poor frames from a temporal confocal
% and visible series, while simultaneously updating the acceptable_frames
% file.

clear
close all force
clc

%% Filename determination and handling
[confocal_fname, pathname] = uigetfile('*.avi', 'Select the confocal temporal video');

confind = strfind(confocal_fname,'confocal');

if isempty(confind)
   error('Could not find confocal in the filename. Needed for proper function of this script!'); 
end

visible_fname = strrep(confocal_fname, 'confocal', 'visible');
split_fname = strrep(confocal_fname, 'confocal', 'split_det');

if exist(fullfile(pathname,split_fname),'file')
   loadsplit = 1;
else
   loadsplit = 0;
end

% Find where the filename should be cut off in the confocal videos, and
% determine our acceptable frame filename.
i=1;
[comb_str remain] = strtok(confocal_fname(confind:end), '_');
acceptable_frame_fname = [];
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

if ~exist(fullfile(pathname, acceptable_frame_fname),'file')
    reply = input('Unable to find acceptable frames csv! Search for it? Y/N [Y]:','s');
    if isempty(reply)
       reply = 'Y';
    end
    
    if strcmpi(reply,'Y')
        [acceptable_frame_fname, af_pathname] = uigetfile(fullfile(pathname, '*.csv'), 'Select the acceptable frames csv.');
    else
        error('Unable to find acceptable frames csv!');
    end
end

acceptable_frame_fname_out = [confocal_fname(1:confind-1) confocal_fname_out '_crop_affine_acceptable_frames.csv'];

% Create our visible output filename
visible_fname_out = [confocal_fname(1:confind-1) confocal_fname_out '_crop_affine.avi'];
visible_fname_out = strrep(visible_fname_out, 'confocal', 'visible');

if loadsplit
    split_fname_out = [confocal_fname(1:confind-1) confocal_fname_out '_crop.avi'];
    split_fname_out = strrep(split_fname_out, 'confocal', 'split_det');
end

% Create our confocal output filename - affine will need to be done outside
% MATLAB.
confocal_fname_out = [confocal_fname(1:confind-1) confocal_fname_out '_crop.avi'];

confocal_vidobj = VideoReader( fullfile(pathname, confocal_fname) );
visible_vidobj = VideoReader( fullfile(pathname, visible_fname) );

if loadsplit
    split_vidobj = VideoReader( fullfile(pathname, split_fname) );
end

%% File loading
confocal_vid = cell(1, round(confocal_vidobj.Duration*confocal_vidobj.FrameRate));
frame_nums = cell(1, round(confocal_vidobj.Duration*confocal_vidobj.FrameRate));
visible_vid = cell(1, round(confocal_vidobj.Duration*confocal_vidobj.FrameRate));

if loadsplit
    split_vid = cell(1, round(confocal_vidobj.Duration*confocal_vidobj.FrameRate));
end

i=1;
while hasFrame(confocal_vidobj)
    
    confocal_vid{i} = readFrame(confocal_vidobj);
    visible_vid{i} = readFrame(visible_vidobj);
    if loadsplit
        split_vid{i} = readFrame(split_vidobj);
    end
    frame_nums{i} = ['Frame ' num2str(i) ' of: ' num2str(size(confocal_vid,2))];
    i=i+1;
end

acc_frame_list = dlmread( fullfile(pathname, acceptable_frame_fname) );

acc_frame_list = sort(acc_frame_list);

if length(acc_frame_list) ~= length(confocal_vid)
   error(['Acceptable frames and confocal video list lengths do not match! (' num2str(length(acc_frame_list)) ' vs ' num2str(length(confocal_vid)) ')']);
end

[selected, cropregion] = linkedImageFrameSingle(confocal_vid, frame_nums);

close all;



% Remove the frames that weren't in the list
confocal_vid = confocal_vid(selected);
visible_vid = visible_vid(selected);
if loadsplit
    split_vid   = split_vid(selected);
end
acc_frame_list = acc_frame_list(selected);

confocal_vid_out = uint8( zeros( cropregion(3,2)-cropregion(2,2)+1, cropregion(2,1)-cropregion(1,1)+1, length(confocal_vid) ));
visible_vid_out = uint8( zeros( cropregion(3,2)-cropregion(2,2)+1, cropregion(2,1)-cropregion(1,1)+1, length(confocal_vid) ));
if loadsplit
    split_vid_out = uint8( zeros( cropregion(3,2)-cropregion(2,2)+1, cropregion(2,1)-cropregion(1,1)+1, length(confocal_vid) ));
end

for i=1:length(confocal_vid)

    confocal_vid_out(:,:,i) = uint8(confocal_vid{i}( cropregion(2,2):cropregion(3,2), cropregion(1,1):cropregion(2,1)));
    visible_vid_out(:,:,i) = uint8(visible_vid{i}( cropregion(2,2):cropregion(3,2), cropregion(1,1):cropregion(2,1)));
    if loadsplit
        split_vid_out(:,:,i) = uint8(split_vid{i}( cropregion(2,2):cropregion(3,2), cropregion(1,1):cropregion(2,1)));
    end

end

outfolder = 'region_cropped';
mkdir(pathname, outfolder);

dlmwrite( fullfile(pathname, outfolder,acceptable_frame_fname_out),acc_frame_list);
confocal_vidobj = VideoWriter( fullfile(pathname, outfolder, confocal_fname_out), 'Grayscale AVI' );
visible_vidobj = VideoWriter( fullfile(pathname, outfolder, visible_fname_out), 'Grayscale AVI' );
if loadsplit
    split_vidobj = VideoWriter( fullfile(pathname, outfolder, split_fname_out), 'Grayscale AVI' );
end


open(confocal_vidobj);
writeVideo(confocal_vidobj,confocal_vid_out);
close(confocal_vidobj);

open(visible_vidobj);
writeVideo(visible_vidobj,visible_vid_out);
close(visible_vidobj);

if loadsplit
    open(split_vidobj);
    writeVideo(split_vidobj,split_vid_out);
    close(split_vidobj);
end

