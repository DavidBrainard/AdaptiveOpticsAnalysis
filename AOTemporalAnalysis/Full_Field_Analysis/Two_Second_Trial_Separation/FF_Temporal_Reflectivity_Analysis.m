function []=FF_Temporal_Reflectivity_Analysis(mov_path, ref_image_fname, stimulus_frames, vid_type)
% []=FF_Temporal_Reflectivity_Analysis_twosource(mov_path, ref_image_fname)
% Robert F Cooper 06-20-2017
%
% @param mov_path - The path to the video which stores the reflectance
% information. Must have an associated file that ends with
% '_acceptable_frames.csv'. That file contains a list of the frame numbers 
% (time points) which were extracted from the original, raw video. If, 
% for example, frame 4 from a 10 frame video did not register, 
% the file would contain the list: [0 1 2 4 5 6 7 8 9].
%
% @param ref_image_fname - The path to the reference image. Must have an
% associated coordinate file that ends in '_coords.csv'.
%
% @param stimulus_frames - A 2 element vector containing the temporal bound
% when the stimulus was delivered, inclusive. (ex: A [67 99] means that a
% stimulus was delivered from the 67th to the 99th frame.
%
% @param vid_type - The type of the video. This can either be the string
% 'control' or 'stimulus'.
%
% This software is responsible for the data processing of temporal
% datasets, by taking the data from separate stim folders and a control
% folders.
% This is designed to work with FULL FIELD datasets, and performs a
% a normalization across the ENTIRE frame, with a standardization relative
% to the prestimulus behavior of EACH cell.
%
% 
% The processing performed in this script is from:
% Cooper RF, Tuten WS, Dubra A, Brainard BH, Morgan JIW. 
% "Non-invasive assessment of human cone photoreceptor function." 
% Biomed Opt Express 8(11): 5098-5112 and are
% encompassed in Figures 1-4A, Equations 1-2.

% *** Constants ***
%
% The shape used to isolate the reflectance at each time point. Box is
% best choice, as it is both the fastest and most reliable.
profile_method = 'box';

% For release, this version only contains the normalizations used in the
% paper. However, the code is structured such that you can add more if desired.
norm_type = 'global_norm_prestim_stdiz'; 

if ~exist('mov_path','var') || ~exist('ref_image_fname','var')
    close all force;
    [ref_image_fname, mov_path]  = uigetfile(fullfile(pwd,'*.tif'));
    stimulus_frames=[67 99];
    
    rply = input('Stimulus (s) or Control (c)? [s]','s');
    if isempty(rply) || strcmpi(rply,'s')
        vid_type='stimulus';
    elseif strcmpi(rply,'c')
        vid_type='control';
    end
end
ref_coords_fname = [ref_image_fname(1:end-4) '_coords.csv'];
stack_fnames = read_folder_contents( mov_path,'avi' );

for i=1:length(stack_fnames)
    if ~isempty( strfind( stack_fnames{i}, ref_image_fname(1:end - length('_AVG.tif') ) ) )
        temporal_stack_fname = stack_fnames{i};
        acceptable_frames_fname = [stack_fnames{i}(1:end-4) '_acceptable_frames.csv'];
        break;
    end
end

%% Load the dataset(s)
ref_image  = double(imread(  fullfile(mov_path, ref_image_fname) ));
stimd_images = dlmread( fullfile(mov_path,acceptable_frames_fname) );
stimd_images = sort(stimd_images)' +1; % For some dumb reason it doesn't store them in the order they're put in the avi.

ref_coords = dlmread( fullfile(mov_path, ref_coords_fname));
ref_coords = [ref_coords(:,1) ref_coords(:,2)]; % If an offset is needed (if, for example you're coming from split coordinates)

temporal_stack_reader = VideoReader( fullfile(mov_path,temporal_stack_fname) );

i=1;
while(hasFrame(temporal_stack_reader))
    temporal_stack(:,:,i) = double(readFrame(temporal_stack_reader));
    i = i+1;
end


%% Find the frames where the stimulus was on, based on the acceptable frames.
stim_times = find(stimd_images>=stimulus_frames(1) & stimd_images<=stimulus_frames(2));

% Make them relative to their actual temporal locations.
stim_locs = stimd_images(stim_times);


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
capillary_masks = repmat(capillary_mask,[1 1 size(temporal_stack,3)]);

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
        case 'segmentation'
            roiradius = 8;
            
            if (ref_coords(i,1) - roiradius) > 1 && (ref_coords(i,1) + roiradius) < size(ref_image,2) &&...
               (ref_coords(i,2) - roiradius) > 1 && (ref_coords(i,2) + roiradius) < size(ref_image,1)

                roi = var_ref_image(ref_coords(i,2) - roiradius : ref_coords(i,2) + roiradius, ref_coords(i,1) - roiradius : ref_coords(i,1) + roiradius);
                
                polarroi = imcart2pseudopolar(roi,.25,2);

                
                diffpim = diff(polarroi,1,2);

                diffpim = diffpim - min(diffpim(:));
                
                [pad_roi, adj, rowcol]=segment_cell(diffpim);

                dg = digraph(adj);        

                shortpath = shortestpath(dg,sub2ind(size(pad_roi), 1, ceil(size(pad_roi,2)/3)), ...
                                            sub2ind(size(pad_roi), size(pad_roi,1), ceil(size(pad_roi,2)/3)) );

                cone_edge_pol = rowcol(shortpath,:);

                cone_edge_pol = cone_edge_pol((cone_edge_pol(:,1) > 3 & cone_edge_pol(:,1) < 360),:);
                
                [x,y] = pol2cart( cone_edge_pol(1:end,1)*pi/180 , ...
                                  cone_edge_pol(1:end,2)/4 );

                conv_inds = convhull(ceil(x),ceil(y));

                cellseg{i} = [x(conv_inds)+ref_coords(i,1), y(conv_inds)+ref_coords(i,2)];

                cellseg_mask = roipoly(ref_image, (cellseg{i}(:,1))+1, (cellseg{i}(:,2)));
                cellseg_inds{i} = find(cellseg_mask~=0);
                
                ref_image(cellseg_inds{i})= -1;
                
                %figure(10); imagesc(ref_image); axis image;
            end
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

%% Extract the raw reflectance of each cell.

coords_used = ref_coords(~cellfun(@isempty,cellseg_inds), :);
cellseg_inds = cellseg_inds(~cellfun(@isempty,cellseg_inds));


cell_reflectance = cell( length(cellseg_inds),1  );
cell_times = cell( length(cellseg_inds),1  );

j=1;

if ~ishandle(wbh)
    wbh = waitbar(0, 'Creating reflectance profile for cell: 0');
end

for i=1:length(cellseg_inds)
    waitbar(i/length(cellseg_inds),wbh, ['Creating reflectance profile for cell: ' num2str(i)]);

    cell_times{i} = stimd_images;
    cell_reflectance{i} = zeros(1, size(temporal_stack,3));
    
    if any( capillary_mask(cellseg_inds{i}) ~= 1 )
        
        cell_reflectance{i} = nan(1,size(temporal_stack,3));
    elseif strcmp(profile_method, 'box') % Shorthand, but way faster than the naive implementation, for box extraction only. (>2x speedup)

        [m,n] = ind2sub(vid_size, cellseg_inds{i});
        
        thisstack = temporal_stack(min(m):max(m), min(n):max(n),:);        
        thisstack(thisstack == 0) = NaN;
        
        thisstack = sum(thisstack,1);
        thisstack = squeeze(sum(thisstack,2));
        thisstack = thisstack./ (length(cellseg_inds{i})*ones(size(thisstack)));
        
        cell_reflectance{i} = thisstack';

    else

         for t=1:size(temporal_stack,3)
             masked_timepoint = temporal_stack(:,:,t); 
             if all( masked_timepoint(cellseg_inds{i}) ~= 0 )
                 cell_reflectance{i}(t) = mean( masked_timepoint(cellseg_inds{i}));
             else            
                 cell_reflectance{i}(t) = NaN;
             end
         end

    end
end
close(wbh);


cell_ref = cell2mat(cell_reflectance);

cellinds = find( ~all(isnan(cell_ref),2) );
coords_used = coords_used(cellinds,:);
cell_ref = cell_ref( cellinds, :); % Remove any cells that have NaNs.

%% Find the means / std devs
for t=1:size(cell_ref,2)
    ref_mean(t) = mean(cell_ref( ~isnan(cell_ref(:,t)) ,t));        
    ref_hist(t,:) = hist(cell_ref( ~isnan(cell_ref(:,t)) ,t),255);    
    ref_std(t) = std(cell_ref( ~isnan(cell_ref(:,t)) ,t));
end


figure(9); plot(ref_mean,'k'); title('Cell Mean Reflectance');
if ~exist( fullfile(mov_path, 'Frame_Mean_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Frame_Mean_Plots'))
end
saveas(gcf, fullfile(mov_path, 'Frame_Mean_Plots' , [ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_cutoff_' norm_type '_' vid_type '_mean_plot.svg' ] ) );

figure(90); plot(ref_std,'k'); title('Cell Reflectance Std Dev');
if ~exist( fullfile(mov_path, 'Frame_Stddev_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Frame_Stddev_Plots'))
end
saveas(gcf, fullfile(mov_path, 'Frame_Stddev_Plots' , [ref_image_fname(1:end - length('_STD_DEV.tif') ) '_' profile_method '_cutoff_' norm_type '_' vid_type '_mean_plot.png' ] ) );


%% Normalization to the mean
norm_cell_reflectance = cell( size(cell_reflectance) );

for i=1:length( cell_reflectance )
    
    if ~isempty( strfind( norm_type, 'global_norm') )
        norm_cell_reflectance{i} = cell_reflectance{i} ./ ref_mean;
    elseif ~isempty( strfind( norm_type, 'no_norm') )
        norm_cell_reflectance{i} = cell_reflectance{i};
        error('No normalization selected!')
    end

    no_ref = ~isnan(norm_cell_reflectance{i});
    
    norm_cell_reflectance{i} = norm_cell_reflectance{i}(no_ref);
    cell_times{i}       = cell_times{i}(no_ref);
    
end


%% Standardization
if ~isempty( strfind(norm_type, 'prestim_stdiz'))
    % Then normalize to the average intensity of each cone BEFORE stimulus.
    prestim_std=nan(1,length( norm_cell_reflectance ));
    prestim_mean=nan(1,length( norm_cell_reflectance ));
    
    for i=1:length( norm_cell_reflectance )

        prestim_mean(i) = mean( norm_cell_reflectance{i}( cell_times{i}<stim_locs(1) & ~isnan( norm_cell_reflectance{i} ) ) );
        
        norm_cell_reflectance{i} = norm_cell_reflectance{i}-prestim_mean(i);
        
        prestim_std(i) = std( norm_cell_reflectance{i}( cell_times{i}<stim_locs(1) & ~isnan( norm_cell_reflectance{i} ) ) );
        
        norm_cell_reflectance{i} = norm_cell_reflectance{i}/( prestim_std(i) ); % /sqrt(length(norm_control_cell_reflectance{i})) );
    end
else
    % NOP - leave stub in case another standardization is created.
end

%% Standard deviation of all cells before first stimulus

thatmax = max( cellfun(@max, cell_times( ~cellfun(@isempty,cell_times)) ) );   

[ ref_stddev,ref_times ] = reflectance_std_dev( cell_times( ~cellfun(@isempty,cell_times) ), ...
                                                norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance) ), thatmax );

clipped_ref_times = [];

if ~isempty(ref_times)
    i=1;
    while i<= length( ref_times )

        % Remove timepoints from cells that are NaN
        if isnan(ref_times(i))
            ref_times(i) = [];
            ref_stddev(i) = [];        
        else
            clipped_ref_times = [clipped_ref_times; ref_times(i)];
            i = i+1;
        end

    end    
end

hz=16.6;

figure(10); 
if ~isempty(ref_stddev)
    plot( clipped_ref_times/hz,ref_stddev,'k'); hold on;
end
if strcmp(vid_type,'stimulus')
    plot(stim_locs/hz, max(ref_stddev)*ones(size(stim_locs)),'r*'); 
end

hold off;
ylabel('Standard deviation'); xlabel('Time (s)'); title( strrep( [ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_stddev_ref_plot' ], '_',' ' ) );
axis([0 15 -1 4])

if ~exist( fullfile(mov_path, 'Std_Dev_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Std_Dev_Plots'))
end
saveas(gcf, fullfile(mov_path, 'Std_Dev_Plots' , [ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_' norm_type '_' vid_type '_stddev.png' ] ) );

if ~exist( fullfile(mov_path, 'Profile_Data'), 'dir' )
    mkdir(fullfile(mov_path, 'Profile_Data'))
end
% Dump all the analyzed data to disk
save(fullfile(mov_path, 'Profile_Data' ,[ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method '_' norm_type '_' vid_type '_profiledata.mat']), ...
     'cell_times', 'norm_cell_reflectance','coords_used','ref_image','ref_mean','ref_stddev','vid_type' );

  
%% Remove the empty cells
norm_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance)  );
cell_times            = cell_times( ~cellfun(@isempty,cell_times) );


figure(11);
if ~isempty(ref_stddev)
    for i=1:length(norm_cell_reflectance)
        plot(cell_times{i}, norm_cell_reflectance{i},'k' ); hold on;
    end
end
hold off;


%% Save the plots
if ~exist( fullfile(mov_path, 'Raw_Plots'), 'dir' )
    mkdir(fullfile(mov_path, 'Raw_Plots'))
end
saveas(gcf, fullfile(mov_path,  'Raw_Plots' ,[ref_image_fname(1:end - length('_AVG.tif') ) '_' profile_method  '_' norm_type '_' vid_type '_raw_plot.png' ] ) );
