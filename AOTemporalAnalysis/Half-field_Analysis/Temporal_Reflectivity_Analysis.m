function []=Temporal_Reflectivity_Analysis(mov_path, ref_image_fname)
% []=TEMPORAL_REFLECTIVITY_ANALYSIS(mov_path, ref_image_fname)
%
%   Performs a temporal analysis of a single trial's dataset. This dataset
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

% Created by Robert F Cooper 11-3-2015 10:40AM

profile_method = 'box';
norm_type = 'regional_norm_prestimminusdiv'; % The normalization approach.
cutoff = 0.9; % The percentage of time a cone must be stimulated relative to all stimulus in order to be included for analysis

if ~exist('mov_path','var') || ~exist('ref_image_fname','var')
    close all force;
    [ref_image_fname, mov_path]  = uigetfile(fullfile(pwd,'*.tif'));
end
ref_coords_fname = [ref_image_fname(1:end-4) '_coords.csv'];
stack_fnames = read_folder_contents( mov_path,'avi' );

for i=1:length(stack_fnames)
    if ~isempty( strfind( stack_fnames{i}, ref_image_fname(1:end - length('_AVG.tif') ) ) )
        temporal_stack_fname = stack_fnames{i};
        acceptable_frames_fname = [stack_fnames{i}(1:end-4) '_acceptable_frames.csv'];
        visible_stack_fname = strrep(temporal_stack_fname,'confocal','visible');
        visible_stack_fname = strrep(visible_stack_fname,'split_det','visible');
        break;
    end
end

%% Load the dataset(s)
ref_image  = double(imread(  fullfile(mov_path, ref_image_fname) ));
accept_images = dlmread( fullfile(mov_path,acceptable_frames_fname) );
accept_images = sort(accept_images)' +1; % For some dumb reason it doesn't store them in the order they're put in the avi.

ref_coords = dlmread( fullfile(mov_path, ref_coords_fname));

ref_coords = [ref_coords(:,1) ref_coords(:,2)];

temporal_stack_reader = VideoReader( fullfile(mov_path,temporal_stack_fname) );
if exist(fullfile(mov_path,visible_stack_fname),'file')
    visible_stack_reader = VideoReader( fullfile(mov_path,visible_stack_fname) );
end

i=1;
while(hasFrame(temporal_stack_reader))
    temporal_stack(:,:,i) = double(readFrame(temporal_stack_reader));
    if exist(fullfile(mov_path,visible_stack_fname),'file')
        visible_stack(:,:,i)  = readFrame(visible_stack_reader);
    else
        visible_stack(:,:,i)  = zeros(size(temporal_stack(:,:,i)));
    end
    i = i+1;
end


%% Find the frames where the stimulus was on, and create the stimulus
% "masks"
visible_signal = zeros(size(visible_stack,3),1);

noise_floor = 20; %mean(max(visible_stack(:,:,1)));
for v=1:size(visible_stack,3)
    
    vis_frm = visible_stack(:,:,v);
    visible_signal(v) = mean(vis_frm(:));
end

stim_times = find( visible_signal > noise_floor );
% Make them relative to their actual temporal locations.
stim_locs = accept_images(stim_times);

% Make the masks where visible light fell (and was detected)
vis_masks = zeros(size(visible_stack,1), size(visible_stack,2), size(stim_locs,1)+1);

% The first frame recieves no stimulus
for i=1:length(stim_locs)
    
    vis_frm = visible_stack(:,:,stim_times(i));        

    vis_frm( vis_frm<=2*noise_floor ) = 0;
    vis_frm( vis_frm>2*noise_floor )  = 1;

    vis_masks(:,:,i+1) = imclose(vis_frm, strel('disk',13) );
    vis_masks(:,:,i+1) = imopen(vis_masks(:,:,i+1), strel('disk',9) );

end

if ~isempty(stim_locs) % If there were stimulus frames, find them and set up the masks to use, as well as the normalization frames

    stim_mask = max(vis_masks,[],3);
    
    control_mask = ~stim_mask; % The control region is where no stimulus fell.
    
    vis_masks = sum(vis_masks,3);
    
    
    figure(2); imagesc(vis_masks); axis image;
    
    % Only take regions with more than 90% of the stimulus falling on it.
    vis_masks( vis_masks < cutoff*max(vis_masks(:)) ) = 0;
    
    stim_mask( vis_masks < cutoff*max(vis_masks(:)) ) = 0;
    
else % If there were no detected stimuli frames.
    
    stim_mask  = zeros( size(temporal_stack,1), size(temporal_stack,2) );
    stim_locs = 67; %size(temporal_stack,3)+1;
    vis_masks = [];
    control_mask = ~stim_mask;
    
end

%% Create the capillary mask- only use the data before the stimulus fires to do so.
if ~isempty(stim_times)
    capillary_mask = double(~tam_etal_capillary_func( temporal_stack(:,:,1:stim_times(1)-1) ));
else
    capillary_mask = double(~tam_etal_capillary_func( temporal_stack(:,:,1:66) ));
end

if ~exist( fullfile(mov_path, 'Capillary_Maps'), 'dir' )
    mkdir(fullfile(mov_path, 'Capillary_Maps'))
end
imwrite(uint8(capillary_mask*255), fullfile(mov_path, 'Capillary_Maps' ,[ref_image_fname(1:end - length('_AVG.tif') ) '_cap_map.png' ] ) );

% Mask the images to exclude zones with capillaries on top of them.
capillary_masks = repmat(capillary_mask, [1 1 size(temporal_stack,3)]);

var_ref_image = stdfilt(ref_image,ones(3));

ref_image = ref_image.*capillary_mask;
temporal_stack = temporal_stack.*capillary_masks;

%% Isolate individual profiles
ref_coords = round(ref_coords);

cellseg = cell(size(ref_coords,1),1);
cellseg_inds = cell(size(ref_coords,1),1);

wbh = waitbar(0,'Segmenting coordinate 0');
for i=1:size(ref_coords,1)

    waitbar(i/size(ref_coords,1),wbh, ['Segmenting coordinate ' num2str(i)]);
    
    switch( profile_method )
        case 'box'
            roiradius = 1;
            
            if (ref_coords(i,1) - roiradius) > 1 && (ref_coords(i,1) + roiradius) < size(ref_image,2) &&...
               (ref_coords(i,2) - roiradius) > 1 && (ref_coords(i,2) + roiradius) < size(ref_image,1)
           
                [R, C ] = meshgrid((ref_coords(i,2) - roiradius) : (ref_coords(i,2) + roiradius), ...
                                   (ref_coords(i,1) - roiradius) : (ref_coords(i,1) + roiradius));
           
                cellseg_inds{i} = sub2ind( size(ref_image), R, C );

                cellseg_inds{i} = cellseg_inds{i}(:);
                
                ref_image(cellseg_inds{i})= -1;
            end
        case 'cross'
            roiradius = 2;
            
            if (ref_coords(i,1) - roiradius) > 1 && (ref_coords(i,1) + roiradius) < size(ref_image,2) &&...
               (ref_coords(i,2) - roiradius) > 1 && (ref_coords(i,2) + roiradius) < size(ref_image,1)
           
                cellseg_inds{i} = sub2ind( size(ref_image), (ref_coords(i,2) - roiradius) : (ref_coords(i,2) + roiradius), repmat(ref_coords(i,1), [1 2*roiradius+1]) );
                
                cellseg_inds{i} = [cellseg_inds{i} sub2ind( size(ref_image), repmat(ref_coords(i,2), [1 2*roiradius+1]) , ref_coords(i,1) - roiradius : ref_coords(i,1) + roiradius) ];

                cellseg_inds{i} = cellseg_inds{i}(:);
                
                ref_image(cellseg_inds{i})= -1;
            end
    end

end

%% Code for viewing the segmented/masked cones
colorcoded_im = repmat(ref_image,[1 1 3]);

max_overlap = max(vis_masks(:));
max_red_mult = max(ref_image(:))/max_overlap;
            
seg_mask = (ref_image == -1);

if ~isempty( vis_masks )
                        
    colorcoded_im(:,:,1) = colorcoded_im(:,:,1) + (capillary_mask.*seg_mask.* (vis_masks.*max_red_mult));
    colorcoded_im(:,:,3) = colorcoded_im(:,:,3) + (capillary_mask.*seg_mask.* (control_mask.*max(ref_image(:))) );

else
    colorcoded_im(:,:,3) = colorcoded_im(:,:,3) + (capillary_mask.*seg_mask.* (control_mask.*max(ref_image(:))) );
end

cropfig = figure(1); 
imagesc( uint8(colorcoded_im) ); axis image; %hold on; plot(ref_coords(:,1),ref_coords(:,2),'.'); hold off;

if ~exist( fullfile(mov_path, 'Stim_Maps'), 'dir' )
    mkdir(fullfile(mov_path, 'Stim_Maps'))
end
imwrite(uint8(colorcoded_im), fullfile(mov_path, 'Stim_Maps' ,[ref_image_fname(1:end - length('_AVG.svg') ) '_stim_map.png' ] ) );

%% Crop each area to a certain size
% cropsize = 7.5;
% croprect = [0 0 cropsize/.4678 cropsize/.4678];
% 
% figure(cropfig); 
% title('Select the crop region for the STIMULUS cones');
% h=imrect(gca, croprect);
% stim_mask_rect = wait(h);
% close(cropfig)
% stim_mask = poly2mask([stim_mask_rect(1) stim_mask_rect(1)                   stim_mask_rect(1)+stim_mask_rect(3) stim_mask_rect(1)+stim_mask_rect(3) stim_mask_rect(1)],...
%                       [stim_mask_rect(2) stim_mask_rect(2)+stim_mask_rect(4) stim_mask_rect(2)+stim_mask_rect(4) stim_mask_rect(2)                   stim_mask_rect(2)],...
%                       size(colorcoded_im,1), size(colorcoded_im,2));
% 
% cropfig = figure(1); 
% imagesc( uint8(colorcoded_im) ); axis image; title('Select the crop region for the CONTROL cones');
% h=imrect(gca, croprect);
% control_mask_rect = wait(h);
% close(cropfig)
% control_mask = poly2mask([control_mask_rect(1) control_mask_rect(1)                      control_mask_rect(1)+control_mask_rect(3) control_mask_rect(1)+control_mask_rect(3) control_mask_rect(1)],...
%                          [control_mask_rect(2) control_mask_rect(2)+control_mask_rect(4) control_mask_rect(2)+control_mask_rect(4) control_mask_rect(2)                      control_mask_rect(2)],...
%                          size(colorcoded_im,1), size(colorcoded_im,2));

%% Extract the raw reflectance of each cell.
coords_used = ref_coords(~cellfun(@isempty,cellseg_inds), :);
cellseg_inds = cellseg_inds(~cellfun(@isempty,cellseg_inds));

stim_cell_reflectance = cell( length(cellseg_inds),1  );
stim_cell_times = cell( length(cellseg_inds),1  );

control_cell_reflectance = cell( length(cellseg_inds),1  );
control_cell_times = cell( length(cellseg_inds),1  );

j=1;

if ~ishandle(wbh)
    wbh = waitbar(0, 'Creating reflectance profile for cell: 0');
end

usefulmask = stim_mask+control_mask;

for i=1:length(cellseg_inds)
    waitbar(i/length(cellseg_inds),wbh, ['Creating reflectance profile for cell: ' num2str(i)]);

    stim_cell_times{i} = accept_images;
    control_cell_times{i} = accept_images;
    
    stim_cell_reflectance{i}    = zeros(1, size(temporal_stack,3));
    control_cell_reflectance{i} = zeros(1, size(temporal_stack,3));
    
    if any( capillary_mask(cellseg_inds{i}) ~= 1 ) || any( usefulmask(cellseg_inds{i}) ~= 1 )
        stim_cell_reflectance{i} = nan(1,size(temporal_stack,3));
        control_cell_reflectance{i} = nan(1,size(temporal_stack,3));
    else
        
        if all( stim_mask(cellseg_inds{i}) == 1 )
            control_cell_reflectance{i} = nan(1,size(temporal_stack,3));
            
            for t=1:size(temporal_stack,3)
                stim_masked_timepoint = temporal_stack(:,:,t);

                if all( stim_masked_timepoint(cellseg_inds{i}) ~= 0 )
                    % Store if a cell is a stimulus-region cell, or if it is a
                    % control region cell?
                    stim_cell_reflectance{i}(t) = mean( stim_masked_timepoint(cellseg_inds{i}));% ./  mean(stim_norm_timepoint(cellseg_inds{i}) );            
                else
                    stim_cell_reflectance{i}(t) = NaN;
                end
            end
        end
        
        if all( control_mask(cellseg_inds{i}) == 1 )
            stim_cell_reflectance{i} = nan(1,size(temporal_stack,3));
            for t=1:size(temporal_stack,3)
                control_masked_timepoint = temporal_stack(:,:,t);
                if all( control_masked_timepoint(cellseg_inds{i}) ~= 0 )
                    control_cell_reflectance{i}(t) = mean( control_masked_timepoint(cellseg_inds{i}));
                else
                    control_cell_reflectance{i}(t) =  NaN;
                end
            end            
        end
    end
end
close(wbh);


%% Normalize the intensities of each cone to the average value of the control cones
c_cell_ref = cell2mat(control_cell_reflectance);
c_cell_ref = c_cell_ref( ~all(isnan(c_cell_ref),2), :);
s_cell_ref = cell2mat(stim_cell_reflectance);

% Get the indexes of the stimulated cells.
stimcellinds = find( ~all(isnan(s_cell_ref),2) );
s_cell_ref = s_cell_ref( ~all(isnan(s_cell_ref),2), :);


%% Find the means
for t=1:size(c_cell_ref,2)
    c_ref_mean(t) = mean(c_cell_ref( ~isnan(c_cell_ref(:,t)) ,t));    
    s_ref_mean(t) = mean(s_cell_ref( ~isnan(s_cell_ref(:,t)) ,t));
    a_ref_mean(t) = mean([c_cell_ref(~isnan(c_cell_ref(:,t)), t); s_cell_ref( ~isnan(s_cell_ref(:,t)), t) ]);
    
    c_ref_hist(t,:) = hist(c_cell_ref( ~isnan(c_cell_ref(:,t)) ,t),255);
    s_ref_hist(t,:) = hist(s_cell_ref( ~isnan(s_cell_ref(:,t)) ,t),255);
    
    c_ref_std(t) = std(c_cell_ref( ~isnan(c_cell_ref(:,t)) ,t));
    s_ref_std(t) = std(s_cell_ref( ~isnan(s_cell_ref(:,t)) ,t));
end

if ~isempty( strfind( norm_type, 'scaled_control_norm') )
    mean_ratio = s_ref_mean./c_ref_mean;
    s_ref_mean = s_ref_mean./mean(mean_ratio);
end

figure(90); plot(c_ref_std,'b'); hold on; plot(s_ref_std,'r'); hold off; title('Stimulus std dev vs control std dev');
legend('Control std dev','Stimulus std dev')
if ~exist( fullfile(mov_path, 'Frame_Stddev_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Frame_Stddev_Plots'))
end
saveas(gcf, fullfile(mov_path, 'Frame_Stddev_Plots' , [ref_image_fname(1:end - length('_STD_DEV.tif') ) '_' profile_method '_cutoff_' norm_type '_' num2str(cutoff*100) '_mean_plot.png' ] ) );

%% Normalization
norm_stim_cell_reflectance = cell( size(stim_cell_reflectance) );
norm_control_cell_reflectance = cell( size(control_cell_reflectance) );

for i=1:length( stim_cell_reflectance )
    
    if ~isempty( strfind( norm_type, 'control_norm') )        
        norm_stim_cell_reflectance{i} = stim_cell_reflectance{i} ./ c_ref_mean;
    elseif ~isempty( strfind( norm_type, 'regional_norm') )
        norm_stim_cell_reflectance{i} = stim_cell_reflectance{i} ./ s_ref_mean;
    elseif ~isempty( strfind( norm_type, 'no_norm') )
        norm_stim_cell_reflectance{i} = stim_cell_reflectance{i};
    elseif ~isempty( strfind( norm_type, 'global_norm') )                
        norm_stim_cell_reflectance{i} = stim_cell_reflectance{i} ./ a_ref_mean;
    elseif ~isempty( strfind( norm_type, 'stimulus_norm') )                
        norm_stim_cell_reflectance{i} = stim_cell_reflectance{i} ./ s_ref_mean;
    else
        error('No normalization selected!')
    end
    
    no_ref = ~isnan(norm_stim_cell_reflectance{i});
    
    norm_stim_cell_reflectance{i} = norm_stim_cell_reflectance{i}(no_ref); 
    stim_cell_times{i}            = stim_cell_times{i}(no_ref);

end


for i=1:length( control_cell_reflectance )
    
    if ~isempty( strfind( norm_type, 'control_norm') )
        norm_control_cell_reflectance{i} = control_cell_reflectance{i} ./ c_ref_mean;
    elseif ~isempty( strfind( norm_type, 'regional_norm') )
        norm_control_cell_reflectance{i} = control_cell_reflectance{i} ./ c_ref_mean;
    elseif ~isempty( strfind( norm_type, 'no_norm') )
        norm_control_cell_reflectance{i} = control_cell_reflectance{i};
    elseif ~isempty( strfind( norm_type, 'global_norm') )
        norm_control_cell_reflectance{i} = control_cell_reflectance{i} ./ a_ref_mean;
    elseif ~isempty( strfind( norm_type, 'stimulus_norm') )                
        norm_control_cell_reflectance{i} = control_cell_reflectance{i} ./ s_ref_mean;
    else
        error('No normalization selected!')
    end

    no_ref = ~isnan(norm_control_cell_reflectance{i});
    
    norm_control_cell_reflectance{i} = norm_control_cell_reflectance{i}(no_ref);
    control_cell_times{i}       = control_cell_times{i}(no_ref);

end


%% Standardization
if ~isempty( strfind(norm_type, 'prestimminusdiv'))
    % Then normalize to the average intensity of each cone BEFORE stimulus.
    for i=1:length( norm_stim_cell_reflectance ) % STIM
        
        prestim_mean(i) = mean( norm_stim_cell_reflectance{i}(stim_cell_times{i}<stim_locs(1) & ~isnan( norm_stim_cell_reflectance{i} )) );
        
        norm_stim_cell_reflectance{i} = norm_stim_cell_reflectance{i}-prestim_mean(i);
        
        prestim_std(i) = std( norm_stim_cell_reflectance{i}(stim_cell_times{i}<stim_locs(1) & ~isnan( norm_stim_cell_reflectance{i} )) );
        
        norm_stim_cell_reflectance{i} = norm_stim_cell_reflectance{i}/( prestim_std(i) ); % /sqrt(length(norm_stim_cell_reflectance{i})) );
        plot(norm_stim_cell_reflectance{i});
    end
    prestim_mean_stim = prestim_mean;
    prestim_std_stim = prestim_std;
    prestim_std=[];
    prestim_mean=[];
    
    for i=1:length( norm_control_cell_reflectance ) % CONTROL

        prestim_mean(i) = mean( norm_control_cell_reflectance{i}( control_cell_times{i}<stim_locs(1) & ~isnan( norm_control_cell_reflectance{i} ) ) );
        
        norm_control_cell_reflectance{i} = norm_control_cell_reflectance{i}-prestim_mean(i);
        
        prestim_std(i) = std( norm_control_cell_reflectance{i}( control_cell_times{i}<stim_locs(1) & ~isnan( norm_control_cell_reflectance{i} ) ) );
        
        norm_control_cell_reflectance{i} = norm_control_cell_reflectance{i}/( prestim_std(i) ); % /sqrt(length(norm_control_cell_reflectance{i})) );
    end
    prestim_mean_control = prestim_mean;
    prestim_std_control = prestim_std;

elseif ~isempty( strfind(norm_type, 'prestimminus'))
    % Then normalize to the average intensity of each cone BEFORE stimulus.
    for i=1:length( norm_stim_cell_reflectance ) % STIM

        prestim_mean(i) = mean( norm_stim_cell_reflectance{i}(stim_cell_times{i}<stim_locs(1) & ~isnan( norm_stim_cell_reflectance{i} )) );

        norm_stim_cell_reflectance{i} = (norm_stim_cell_reflectance{i}-prestim_mean(i));
    end
    prestim_std=[];
    prestim_mean=[];
    
    for i=1:length( norm_control_cell_reflectance ) % CONTROL

        prestim_mean(i) = mean( norm_control_cell_reflectance{i}( control_cell_times{i}<stim_locs(1) & ~isnan( norm_control_cell_reflectance{i} ) ) );        

        norm_control_cell_reflectance{i} = (norm_control_cell_reflectance{i}-prestim_mean(i));
    end
elseif ~isempty( strfind(norm_type, 'minus'))
    % Then normalize to the  intensity of each cone at it's starting point.
    for i=1:length( norm_stim_cell_reflectance ) % STIM
        if ~isempty( norm_stim_cell_reflectance{i} )
            norm_stim_cell_reflectance{i} = (norm_stim_cell_reflectance{i}-norm_stim_cell_reflectance{i}(1));
        end
    end
    
    for i=1:length( norm_control_cell_reflectance ) % CONTROL
        if ~isempty( norm_control_cell_reflectance{i})
            norm_control_cell_reflectance{i} = (norm_control_cell_reflectance{i}-norm_control_cell_reflectance{i}(1));
        end
    end

            
elseif ~isempty( strfind(norm_type, 'prestim'))
    % Then normalize to the average intensity of each cone BEFORE stimulus.
    for i=1:length( norm_stim_cell_reflectance ) % STIM

        prestim_mean = mean( norm_stim_cell_reflectance{i}(stim_cell_times{i}<stim_locs(1) & ~isnan( norm_stim_cell_reflectance{i} )) );

        norm_stim_cell_reflectance{i} = norm_stim_cell_reflectance{i}./prestim_mean;
    end
    for i=1:length( norm_control_cell_reflectance ) % CONTROL

        prestim_mean = mean( norm_control_cell_reflectance{i}( control_cell_times{i}<stim_locs(1) & ~isnan( norm_control_cell_reflectance{i} ) ) );

        norm_control_cell_reflectance{i} = norm_control_cell_reflectance{i}./prestim_mean;

    end    
elseif ~isempty( strfind(norm_type, 'poststim'))
    % Then normalize to the average intensity of each cone AFTER stimulus.
    for i=1:length( norm_stim_cell_reflectance ) % STIM

        prestim_mean = mean( norm_stim_cell_reflectance{i}(stim_cell_times{i}>stim_locs(end) & ~isnan( norm_stim_cell_reflectance{i} )) );

        norm_stim_cell_reflectance{i} = norm_stim_cell_reflectance{i}./prestim_mean;
    end
    for i=1:length( norm_control_cell_reflectance ) % CONTROL

        prestim_mean = mean( norm_control_cell_reflectance{i}( control_cell_times{i}>stim_locs(end) & ~isnan( norm_control_cell_reflectance{i} ) ) );

        norm_control_cell_reflectance{i} = norm_control_cell_reflectance{i}./prestim_mean;

    end    
elseif ~isempty( strfind(norm_type, 'meanall'))
    % Then normalize to the average intensity of each cone's average value.
    for i=1:length( norm_stim_cell_reflectance ) % STIM

        prestim_mean = mean( norm_stim_cell_reflectance{i}(~isnan( norm_stim_cell_reflectance{i} )) );

        norm_stim_cell_reflectance{i} = norm_stim_cell_reflectance{i}./prestim_mean;
    end
    for i=1:length( norm_control_cell_reflectance ) % CONTROL

        prestim_mean = mean( norm_control_cell_reflectance{i}( ~isnan( norm_control_cell_reflectance{i} ) ) );

        norm_control_cell_reflectance{i} = norm_control_cell_reflectance{i}./prestim_mean;
    end
else
    
end

%% Standard deviation of all cells before first stimulus

% Find their max sizes    
thatstimmax = max( cellfun(@max,stim_cell_times( ~cellfun(@isempty,stim_cell_times) ) ) );
thatcontrolmax = max( cellfun(@max, control_cell_times( ~cellfun(@isempty,control_cell_times)) ) );   

[ ref_stddev_stim, ref_stim_times ] = reflectance_std_dev( stim_cell_times( ~cellfun(@isempty,stim_cell_times) ), ...
                                                              norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) ), thatstimmax );
[ ref_stddev_control, ref_control_times ] = reflectance_std_dev( control_cell_times( ~cellfun(@isempty,control_cell_times) ), ...
                                                                norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance) ), thatcontrolmax );

ref_times = [];

if ~isempty(ref_stim_times) && ~isempty(ref_control_times)
    i=1;
    while i<= length( ref_control_times )

        % Remove times from both stim and control that are NaN
        if isnan(ref_stim_times(i)) || isnan(ref_control_times(i))

            ref_stim_times(i) = [];
            ref_control_times(i) = [];

            ref_stddev_stim(i) = [];
            ref_stddev_control(i) = [];        
        else
            ref_times = [ref_times; ref_stim_times(i)];
            i = i+1;
        end
    end
    
elseif ~isempty(ref_stim_times)
    i=1;
    while i<= length( ref_stim_times )

        % Remove times from both stim and control that are NaN
        if isnan(ref_stim_times(i))
            ref_stim_times(i) = [];
            ref_stddev_stim(i) = [];
        else
            ref_times = [ref_times; ref_stim_times(i)];
            i = i+1;
        end
    end
    
elseif ~isempty(ref_control_times)
    i=1;
    while i<= length( ref_control_times )

        % Remove times from both stim and control that are NaN
        if isnan(ref_control_times(i))
            ref_control_times(i) = [];
            ref_stddev_control(i) = [];
        else
            ref_times = [ref_times; ref_control_times(i)];
            i = i+1;
        end
    end
end

%%
% If its in the normalization, subtract the control value from the stimulus
% value
if ~isempty( strfind(norm_type, 'sub') )
    if ~isempty(ref_stddev_stim)
        ref_stddev_stim    = ref_stddev_stim-ref_stddev_control;
    end
    if ~isempty(ref_stddev_control)
        ref_stddev_control = ref_stddev_control-ref_stddev_control;
    end
end

hz=16.6;

figure(9); plot(ref_times/hz,100*((s_ref_mean-s_ref_mean(1))./s_ref_mean(1)),'r'); axis([0 15 -100 100]);
if ~exist( fullfile(mov_path, 'Frame_Mean_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Frame_Mean_Plots'))
end
saveas(gcf, fullfile(mov_path, 'Frame_Mean_Plots' , [ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_cutoff_' norm_type '_' num2str(cutoff*100) '_mean_plot.svg' ] ) );

figure(10); hold off;
if ~isempty(ref_stddev_stim)
    plot( ref_times/hz,ref_stddev_stim,'r'); hold on;
end

if ~isempty(ref_stddev_control)
    plot( ref_times/hz,ref_stddev_control,'b'); hold on;
end
legend('Stimulus cones','Control cones');
plot(stim_locs/hz, max([ref_stddev_stim; ref_stddev_control])*ones(size(stim_locs)),'r*'); hold off;
ylabel('Standard deviation'); xlabel('Time (s)'); title( strrep( [ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_stddev_ref_plot' ], '_',' ' ) );
axis([0 15 -1 3])

%% Mean to starting value correlation
mean_ratio = s_ref_mean./c_ref_mean;


if ~exist( fullfile(mov_path, 'Std_Dev_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Std_Dev_Plots'))
end
saveas(gcf, fullfile(mov_path, 'Std_Dev_Plots' , [ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_cutoff_' norm_type '_stddev.png' ] ) );

if ~exist( fullfile(mov_path, 'Mat_Profile_Data'), 'dir' )
    mkdir(fullfile(mov_path, 'Mat_Profile_Data'))
end

if ~isempty(stim_times)
    prestim_ref = ref_stddev_stim(ref_times<stim_locs(1) & ~isnan( ref_stddev_stim ));
    prestim_ratio = mean_ratio( ref_times<stim_locs(1) & ~isnan(ref_stddev_stim) )-1;
    save(fullfile(mov_path, 'Mat_Profile_Data' ,[ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_cutoff_' norm_type '_' num2str(cutoff*100) '_profiledata.mat']), 'stim_cell_times', 'norm_stim_cell_reflectance', ...
      'control_cell_times', 'norm_control_cell_reflectance','stimcellinds','mean_ratio','prestim_ref','prestim_ratio','c_ref_mean','s_ref_mean' );
else
    save(fullfile(mov_path, 'Mat_Profile_Data' ,[ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_cutoff_' norm_type '_' num2str(cutoff*100) '_profiledata.mat']), 'stim_cell_times', 'norm_stim_cell_reflectance', ...
      'control_cell_times', 'norm_control_cell_reflectance','stimcellinds','mean_ratio','c_ref_mean','s_ref_mean' );
end


  
%% Remove the empty control cells
norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );

figure(11);
if ~isempty(ref_stddev_control)
    for i=1:length(norm_control_cell_reflectance) % Plot raw
        plot(control_cell_times{i}/hz, norm_control_cell_reflectance{i}-norm_control_cell_reflectance{i}(1),'b' ); hold on;
    end
end

%% Remove the empty stim cells
norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );

figure(11); hold on;
if ~isempty(ref_stddev_stim)
    for i=1:length(norm_stim_cell_reflectance) % Plot raw
        plot(stim_cell_times{i}/hz, norm_stim_cell_reflectance{i}-norm_stim_cell_reflectance{i}(1),'r' ); 
    end
end

%% Save the plots
if ~exist( fullfile(mov_path, 'Raw_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Raw_Plots'))
end
saveas(gcf, fullfile(mov_path,  'Raw_Plots' ,[ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method  '_cutoff_' norm_type '_' num2str(cutoff*100) '_raw_plot.svg' ] ) );
