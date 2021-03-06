% Robert F Cooper 16-19-2017
%
% This script enables a user to remove poor frames from a temporal confocal
% and
% series, while simultaneously updating the acceptable_frames
% file. With the inclusion of the shutter, there is no need to track
% visible frames.

clear
close all force
clc

%% Filename determination and handling
[confocal_fname, pathname] = uigetfile('*.avi', 'Select the confocal temporal video', 'MultiSelect','on');

if ~iscell(confocal_fname)
    confocal_fname={confocal_fname};
end

for k=1:length(confocal_fname)

    confind = strfind(confocal_fname{k},'confocal');

    if isempty(confind)
       error('Could not find confocal in the filename. Needed for proper function of this script!'); 
    end

    split_fname = strrep(confocal_fname{k}, 'confocal', 'split_det');

    if exist(fullfile(pathname,split_fname),'file')
       loadsplit = 1;
    else
       loadsplit = 0;
    end

    % Find where the filename should be cut off in the confocal videos, and
    % determine our acceptable frame filename.
    i=1;
    [comb_str remain] = strtok(confocal_fname{k}(confind:end), '_');
    acceptable_frame_fname = [];
    while ~isempty(remain)
        [tok remain] = strtok( remain, '_');

        if i==5
            confocal_fname_out = comb_str;
        elseif i==8
            acceptable_frame_fname = comb_str;
            if exist([confocal_fname{k}(1:confind-1) acceptable_frame_fname '_acceptable_frames.csv'],'file')
                break;
            end
        elseif i==9
            acceptable_frame_fname = comb_str;
            if exist([confocal_fname{k}(1:confind-1) acceptable_frame_fname '_acceptable_frames.csv'],'file')
                break;
            end
        end

        comb_str = [comb_str '_' tok];

        i=i+1;
    end
    % Create our expected acceptable frame filenames
    acceptable_frame_fname = [confocal_fname{k}(1:confind-1) acceptable_frame_fname '_acceptable_frames.csv'];

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

    acceptable_frame_fname_out = [confocal_fname{k}(1:confind-1) confocal_fname_out '_crop_affine'];


    if loadsplit
        split_fname_out = [confocal_fname{k}(1:confind-1) confocal_fname_out '_crop_affine'];
        split_fname_out = strrep(split_fname_out, 'confocal', 'split_det');
    end

    % Create our confocal output filename - affine will need to be done outside
    % MATLAB.
    confocal_fname_out = [confocal_fname{k}(1:confind-1) confocal_fname_out '_crop_affine'];

    confocal_vidobj = VideoReader( fullfile(pathname, confocal_fname{k}) );

    if loadsplit
        split_vidobj = VideoReader( fullfile(pathname, split_fname) );
    end

    %% File loading
    vid_length = round(confocal_vidobj.Duration*confocal_vidobj.FrameRate);

    confocal_vid = cell(1, vid_length);
    confocal_f_mean = zeros(vid_length,1);

    frame_nums = cell(1, vid_length);

    if loadsplit
        split_vid = cell(1, vid_length);
        split_f_mean = zeros(vid_length,1);
    end

    i=1;
    while hasFrame(confocal_vidobj)

        confocal_vid{i} = readFrame(confocal_vidobj);
        confocal_f_mean(i) = mean(double(confocal_vid{i}(confocal_vid{i}~=0)));
        if loadsplit
            split_vid{i} = readFrame(split_vidobj);
            split_f_mean(i) = mean(double(split_vid{i}(split_vid{i}~=0)));
        end
        frame_nums{i} = ['Frame ' num2str(i) ' of: ' num2str(size(confocal_vid,2))];
        i=i+1;
    end


    acc_frame_list = dlmread( fullfile(pathname, acceptable_frame_fname) );

    acc_frame_list = sort(acc_frame_list);

    if length(acc_frame_list) < length(confocal_vid)
       error(['Acceptable frames and confocal video list lengths do not match! (' num2str(length(acc_frame_list)) ' vs ' num2str(length(confocal_vid)) ')']);
    elseif length(acc_frame_list) > length(confocal_vid)
        acc_frame_list = acc_frame_list(1:length(confocal_vid));
    end

    %% Filter by image mean
    confocal_mean = mean(confocal_f_mean);
    confocal_dev = std(confocal_f_mean);
    if loadsplit
        split_mean = mean(split_f_mean);
        split_dev = std(split_f_mean);
    end

    contenders = false(1,length(frame_nums));
    for n=1:length(frame_nums)        
        contenders(n) =  (confocal_f_mean(n) > confocal_mean-2*confocal_dev);

        if loadsplit
           contenders(n) = contenders(n) & (split_f_mean(n) > split_mean-2*split_dev); 
        end
    end

    % Remove frames from contention.
    confocal_vid = confocal_vid(contenders);
    acc_frame_list = acc_frame_list(contenders);
    if loadsplit
        split_vid = split_vid(contenders);
    end


    %% Determine which frames have divisions
    frm_rp = cell(length(confocal_vid),1);
    cc_areas = cell(length(confocal_vid),1);
    contenders = false(length(confocal_vid),1);
    div_frms = false(length(confocal_vid),1);

    for n=1:length(confocal_vid)
        frm_nonzeros = (confocal_vid{n}>0);

        frm_nonzeros = imclose(frm_nonzeros, ones(5)); % There will always be noise in a simple threshold- get rid of it.

        frm_cc = bwconncomp(frm_nonzeros);    
        frm_rp{n} = regionprops(frm_cc,'Area','BoundingBox');

        cc_areas{n} = [frm_rp{n}.Area];

        % Assuming components that aren't more than an area of a few lines are
        % noise.
        big_comps = cc_areas{n} > size(confocal_vid{n},2)*6;

        small_comps = frm_rp{n}(~big_comps);

        % Mask out any small noisy areas
        for c=1:length(small_comps)
            maskbox = small_comps(c).BoundingBox;
            maskbox = round(maskbox);
            confocal_vid{n}(maskbox(2):(maskbox(2)+maskbox(4)), maskbox(1):(maskbox(1)+maskbox(3)) ) = 0;
        end

        % Remove the small areas from consideration.
        frm_rp{n} = frm_rp{n}(big_comps);
        cc_areas{n} = cc_areas{n}(big_comps);
    end

    mean_cc_area = mean([cc_areas{:}]);
    std_cc_area = std([cc_areas{:}]);

    for n=1:length(cc_areas)

        imagesc( imclose(confocal_vid{n}>0, ones(5)) ); colormap gray;
        % Remove components that aren't more than mu-2std dev pixels in area    
        big_enough_comps = cc_areas{n} > mean_cc_area-std_cc_area;

        % Put these in a list to track if we have to we can use them, but
        % dropping the smaller of the two components
        if sum(big_enough_comps) == 1 && length(big_enough_comps) == 1
            contenders(n) = 1;

        % If it has breaks in it (and has at least one piece that's big enough)
        % flag it to determine if it is worth keeping.
        elseif sum(big_enough_comps) > 0 && sum(big_enough_comps) <= length(big_enough_comps)
            div_frms(n) = 1;        
        end
    end

    % Pull out the divided frames and see how workable they are.
    div_confocal_vid = confocal_vid(div_frms);
    div_acc_frame_list = acc_frame_list(div_frms);
    if loadsplit
        div_split_vid = split_vid(div_frms);
    end

    % Remove the divided frames from contention (for now).
    confocal_vid = confocal_vid(contenders);
    acc_frame_list = acc_frame_list(contenders);
    if loadsplit
        split_vid = split_vid(contenders);
    end

    %% Make the ideal area mask from the unregistered videos.

    contenders = false(length(confocal_vid),1);

    % Make a sum map of all of the undivided frames to determine the ideal
    % cropping area.
    crop_mask = zeros(size(confocal_vid{1}));
    sum_map = zeros(size(confocal_vid{1}));

    for n=1:length(confocal_vid)
        frm_nonzeros = imclose((confocal_vid{n}>0), ones(5)); 
        sum_map = sum_map+frm_nonzeros;
    end
    % imagesc(sum_map); axis image;

    max_frm_mask = sum_map >= ( max(sum_map(:))*0.8 );
    [C, h, w, max_largest_rect] =FindLargestRectangles(max_frm_mask,[1 1 0], [300 150]);
    cropregion = regionprops(max_largest_rect,'BoundingBox');
    cropregion = floor(cropregion.BoundingBox);
    % Bound our crop region to where the sum_map actually exists
    cropregion(cropregion<1) = 1;
    if cropregion(4)>size(sum_map,1)
        cropregion(4) = size(sum_map,1);
    end
    if cropregion(3)>size(sum_map,2)
        cropregion(3) = size(sum_map,2);
    end    
    
    maxcropregion = [cropregion(1:2), cropregion(1)+cropregion(3), cropregion(2)+cropregion(4)];

    average_frm_mask = sum_map > ceil(mean(sum_map(:)));
    % Find the largest incribed rectangle in this mask.
    [C, h, w, largest_rect] =FindLargestRectangles(average_frm_mask,[1 1 0], [300 150]);

    % Find the coordinates for each corner of the rectangle, and
    % return them
    cropregion = regionprops(largest_rect,'BoundingBox');
    cropregion = round(cropregion.BoundingBox);

    cropregion = [cropregion(1:2), cropregion(1)+cropregion(3), cropregion(2)+cropregion(4)];

    % Bound our crop region to where the sum_map actually exists
    cropregion(cropregion<1) = 1;
    if cropregion(4)>size(sum_map,1)
        cropregion(4) = size(sum_map,1);
    end
    if cropregion(3)>size(sum_map,2)
        cropregion(3) = size(sum_map,2);
    end

    sum_map_crop = sum_map(cropregion(2):cropregion(4),cropregion(1):cropregion(3));

    figure(1); imagesc( sum_map_crop ); axis image; title('Ideal cropped sum map');


    %% Register these frames together, removing residual rotation.

    im_only_vid = cell(length(confocal_vid),1);
    im_only_vid_ref = cell(length(confocal_vid),1);
    split_im_only_vid = cell(length(confocal_vid),1);

    for n=1:length(confocal_vid)

        frm_nonzeros = imclose((confocal_vid{n}>0), ones(11));

        masked_frm = frm_nonzeros.*largest_rect;
        cropped_masked_frm = masked_frm( cropregion(2):cropregion(4),cropregion(1):cropregion(3) );


        im_only_vid{n} = confocal_vid{n}( cropregion(2):cropregion(4),cropregion(1):cropregion(3) );
        if loadsplit
            split_im_only_vid{n} = split_vid{n}( cropregion(2):cropregion(4),cropregion(1):cropregion(3) );
        end
        reg_only_vid{n} = confocal_vid{n}( maxcropregion(2):maxcropregion(4),maxcropregion(1):maxcropregion(3) );
    end


    for n=1:length(im_only_vid)
        prc_filled(n) = sum( sum( imclose((im_only_vid{n}>0), ones(11)) ) )./ (size(im_only_vid{n},1)*size(im_only_vid{n},2));
    end

    max_fill = max(prc_filled);
    for n=1:length(im_only_vid) % Find the first frame that has the highest fill %.
        if prc_filled(n) == max_fill
            ref_ind =n;
            break;
        end
    end

    [optimizer, metric]  = imregconfig('monomodal');
    % optimizer.GradientMagnitudeTolerance = 1e-4;
    % optimizer.MinimumStepLength = 1e-5;
    % optimizer.MaximumStepLength = 0.04;
    % optimizer.MaximumIterations = 100;

    tic;

    % Register the image stack forward. It is more stable if we align to the 
    % frame with the largest percentage of the cropped region covered.
    parfor n=1:length(reg_only_vid)

        % Register using the cropped frame
        forward_reg_tform{n}=imregtform(reg_only_vid{n}, reg_only_vid{ref_ind},'affine',...
                                optimizer, metric,'PyramidLevels',1, 'InitialTransformation', affine2d());%,'DisplayOptimization',true);

        tforms(:,:,n) = forward_reg_tform{n}.T;
    end
    toc;

    mean_tforms = mean(tforms,3);

    % Determine the typical distance from the mean xform.
    for t=1:size(tforms,3)
        tfo = tforms(:,:,t);
        frob_dist(t) = abs(tfo(:)'*mean_tforms(:));
    end

    mean_frob_dist = mean(frob_dist);
    std_frob_dist = std(frob_dist);


    max_frob_dist = 3*std_frob_dist+mean_frob_dist;

    %%
    reg_confocal_vid = cell(length(confocal_vid),1);
    reg_split_vid = cell(length(confocal_vid),1);

    for f=1:length(im_only_vid)    

        if f~=ref_ind

            % We would NOT expect large changes here- so if there is a big
            % jump, use the surrounding transforms to force some stablility on
            % the registration.
            tfo = tforms(:,:,f);
            t=1;
            theend = false;
            while abs(tfo(:)'*mean_tforms(:)) > max_frob_dist
                tfo = tforms(:,:,f-t);
                t=t+1;

                if (f-t) == 0
                   theend = true;
                   break;
                end
            end

            if theend % If we reached the end in that direction, check the other direction
                t=1;
                while abs(tfo(:)'*mean_tforms(:)) > frob_dist  && (f+t) < length(im_only_vid)
    %                 f
                    tfo = tforms(:,:,f+t);
                    t=t+1;
                end
            end

            reg_confocal_vid{f}= imwarp(im_only_vid{f}, affine2d(tfo),'OutputView', imref2d(size(im_only_vid{1})) );
            
            if loadsplit
                reg_split_vid{f}= imwarp(split_im_only_vid{f}, affine2d(tfo),'OutputView', imref2d(size(im_only_vid{1})) );
            end
    %         if f == 42
    %             figure(1); imagesc(im_only_vid{f}); axis image; colormap gray;
    %             figure(2); imagesc(reg_confocal_vid{f}); axis image; colormap gray;
    %             drawnow;
    %         end        


        else
            reg_confocal_vid{f}= im_only_vid{f};
            if loadsplit
                reg_split_vid{f}= split_im_only_vid{f};
            end
        end
    end


    %% Remake the sum map with our rotated data.

    % Make a sum map of all of the undivided frames to determine the ideal
    % cropping area.
    crop_mask = zeros(size(reg_confocal_vid{1}));
    sum_map = zeros(size(reg_confocal_vid{1}));

    for f=1:length(reg_confocal_vid)
        frm_nonzeros = imclose((reg_confocal_vid{f}>0), ones(5)); 
        sum_map = sum_map+frm_nonzeros;
    end


    average_frm_mask = sum_map>mean2(sum_map);
    % Find the largest incribed rectangle in this mask.
    [C, h, w, largest_rect] =FindLargestRectangles(average_frm_mask,[1 1 0], [300 150]);

    % Find the coordinates for each corner of the rectangle, and
    % return them
    cropregion = regionprops(largest_rect,'ConvexHull');
    cropregion = cropregion.ConvexHull;
    % TL TR BR BL
    cropregion = [floor(min(cropregion)); [ceil(max(cropregion(:,1))) floor( min(cropregion(:,2))) ]; ...
                  ceil(max(cropregion));  [floor(min(cropregion(:,1))) ceil( max(cropregion(:,2)))] ];
    cropregion(cropregion==0) = 1;
    cropregion( cropregion(:,1)>size(reg_confocal_vid{1},2),1 ) = size(reg_confocal_vid{1},2);
    cropregion( cropregion(:,2)>size(reg_confocal_vid{1},1),2 ) = size(reg_confocal_vid{1},1);

    figure(2); imagesc(sum_map); axis image; title('Ideal cropped sum map- registered data');




    %% Determine the divided files' degree of overlap with the ideal cropping region, add them to the stack if possible.
    % crop_mask( cropregion(2,2):cropregion(3,2), cropregion(1,1):cropregion(2,1) ) = 1;
    % 
    % div_confocal_vid = confocal_vid(div_frms);
    % div_acc_frame_list = acc_frame_list(div_frms);
    % if loadsplit
    %     div_split_vid = split_vid(div_frms);
    % end
    % 
    % for f=1:length(div_confocal_vid)
    %     frm_nonzeros = (confocal_vid{f}>0);
    %     
    %     frm_nonzeros = imclose(frm_nonzeros, ones(5)); % There will always be noise in a simple threshold- get rid of it.
    %     
    %     frm_cc = bwconncomp(frm_nonzeros);    
    %     frm_rp{f} = regionprops(frm_cc,'Area','BoundingBox');
    %     
    % end


    %% Output the cropped frames
    confocal_vid_out = uint8( zeros( size(reg_confocal_vid{1},1), size(reg_confocal_vid{1},2), length(confocal_vid) ));

    if loadsplit
        split_vid_out = uint8( zeros( size(reg_split_vid{1},1), size(reg_split_vid{1},2), length(confocal_vid) ));
    end

    for i=1:length(confocal_vid)

        confocal_vid_out(:,:,i) = uint8( reg_confocal_vid{i} );

        if loadsplit
            split_vid_out(:,:,i) = uint8( reg_split_vid{i} );
        end

    end

    outfolder = 'region_cropped';
    mkdir(pathname, outfolder);

    frmcount = ['_n' num2str(size(confocal_vid_out,3))];

    dlmwrite( fullfile(pathname, outfolder,[acceptable_frame_fname_out frmcount '_acceptable_frames.csv']),acc_frame_list);
    confocal_vidobj = VideoWriter( fullfile(pathname, outfolder, [confocal_fname_out frmcount '.avi']), 'Grayscale AVI' );

%     if loadsplit
%         split_vidobj = VideoWriter( fullfile(pathname, outfolder, [split_fname_out frmcount '.avi']), 'Grayscale AVI' );
%     end


    open(confocal_vidobj);
    writeVideo(confocal_vidobj,confocal_vid_out);
    close(confocal_vidobj);

%     if loadsplit
%         open(split_vidobj);
%         writeVideo(split_vidobj,split_vid_out);
%         close(split_vidobj);
%     end

    % Write the average images.
    imwrite(uint8(sum(confocal_vid_out,3)./sum_map), fullfile(pathname, outfolder, [confocal_fname_out frmcount '_AVG.tif']) );

    if loadsplit
        imwrite(uint8(sum(split_vid_out,3)./sum_map), fullfile(pathname, outfolder, [split_fname_out frmcount '_AVG.tif']) );
    end
    
    close all;
end
