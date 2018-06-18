function [ reference_frames, totalNumFrames ] = extract_candidate_reference_frames( fPathName, desinusoid_matrix, STRIP_SIZE, BAD_STRIP_THRESHOLD, NUM_FRM_PER_GROUP )
%EXTRACT_CANDIDATE_REFERENCE_FRAMES Extracts a list of candidate reference
% frames from a input video.
%   [reference_frames] = EXTRACT_CANDIDATE_REFERENCE_FRAMES( fName, pName,desinusoid_matrix, STRIP_SIZE, BAD_STRIP_THRESHOLD, NUM_FRM_PER_GROUP )


    vidReader = VideoReader( fPathName );
    
    i=1;
    while(hasFrame(vidReader))
        image_stack(:,:,i) = uint8(double(readFrame(vidReader))* desinusoid_matrix);
        
        frame_mean(i) = mean2(image_stack(:,:,i));
        i = i+1;
    end
    numFrames = i-1;
    
    % Get some basic heuristics from each modality.
    mode_mean = mean(frame_mean(:));
    mode_dev = std(frame_mean(:));

    frame_contenders = (1:numFrames);
    totalNumFrames = numFrames;

    strip_inds = 0:STRIP_SIZE:size(image_stack(:,:,1),1);
    strip_inds(1) = 1;
    if strip_inds(end) ~= size(image_stack,1)
        if (size(image_stack,1)-strip_inds(end)) > STRIP_SIZE/2
            strip_inds = [strip_inds size(image_stack,1)];
        else
            strip_inds(end) = size(image_stack,1);
        end
    end
    num_strips = length(strip_inds)-1;
    
    %% Filter by image mean
    mean_contenders = false(1,numFrames);
    for f=1:numFrames        
        mean_contenders(f) =  (frame_mean(f) < mode_mean+2*mode_dev);% &&...
                              %(frame_mean(f,modalityInd) > mode_mean(modalityInd)-2*mode_dev(modalityInd));
    end
    
    frame_contenders = frame_contenders(mean_contenders);
    
    numFrames = length(frame_contenders);
    
    radon_bandwidth = zeros(numFrames,num_strips);
    
    disp(['Filtered by image mean... ' num2str(length(frame_contenders)) ' frames remain.']);
    
    
    if exist('parfor','builtin') == 5 % If we can multithread it, do it!
        %% Filter by Radon transform FWHM        
        parfor f=1:numFrames
            
            frame_ind=frame_contenders(f);
            
            for s=1:num_strips
    
                % Get the log power spectrum for us to play with
                pwr_spect = ( abs(fftshift(fft2(image_stack(strip_inds(s):strip_inds(s+1),:, frame_ind),512, 512))).^2);
                % From our padding, the center vertical frequency will be
                % garbage- remove it for our purposes.
                pwr_spect = log10(pwr_spect(:,[1:256 258:512]));
                
                pwr_spect(isnan(pwr_spect(:))) =0;
                pwr_spect(isinf(pwr_spect(:))) =0;
                
                % Threshold is set using the upper 2 std devs
                thresh_pwr_spect = ( pwr_spect>(mean(pwr_spect(:))+2*std(pwr_spect(:))) );
                
                radoned = radon( thresh_pwr_spect );                

                % Determine the minimum and maximum FWHM
                halfmax = repmat(max(radoned)./2,[size(radoned,1) 1]);
                fwhm = sum(radoned>halfmax);

                radon_bandwidth(f,s) = max(fwhm)-min(fwhm);

            end
        end        
        
        threshold = ceil(mean(radon_bandwidth(:))+ 2*std(radon_bandwidth(:)));
        % After thresholding and removal, update the contenders list. One
        % bad strip is not enough to kick out the frame.
        frame_contenders = frame_contenders(~(sum(radon_bandwidth > threshold,2)>BAD_STRIP_THRESHOLD) );        
        
        disp(['Filtered by Radon FWHM... ' num2str(length(frame_contenders)) ' frames remain.']);
        clear radon_bandwidth;
        
        %% Determine NCC between pairs of frames.
        
        contender_image_stack = image_stack(:,:, frame_contenders);
        clear image_stack;
        
        frm1 = double( contender_image_stack(:,:,1) );
        [m, n] = size(frm1);
        paddiffm = (size(frm1,1)*2)-1-m;
        paddiffn = (size(frm1,2)*2)-1-n;

        % Make a mask to remove any edge effects.
        [maskdistx, maskdisty] = meshgrid( 1:(size(frm1,2)*2)-1, 1:(size(frm1,1)*2)-1);
        
        maskdistx = maskdistx-(size(maskdistx,2)/2);
        maskdisty = maskdisty-(size(maskdisty,1)/2);
        
        xcorr_mask = sqrt(maskdistx.^2 + maskdisty.^2) < max(size(frm1))/2;
        % Determine the number of pixels that will overlap between the two
        % images at any given point.
        numberOfOverlapPixels = arfs_local_sum(ones(size(frm1)),size(frm1,1),size(frm1,2));
        
        padfrm1 = padarray(frm1,[paddiffm paddiffn],0,'post');
        
        fft_frm1 =  fft2( padfrm1 );
        
        seq_ncc = zeros(length(frame_contenders),1);
        seq_ncc_offset = zeros(length(frame_contenders)-1,2);
        frame_group = zeros(length(frame_contenders)-1,1);
        fft_ims = zeros(size(padfrm1,1), size(padfrm1,2), length(frame_contenders));
        
        fft_ims(:,:,1) = fft2( padfrm1 );
        
        group = 1;        
        frame_group(1) = group;
        seq_ncc(1) = NaN; % 0 to 1 doesn't align.
        seq_ncc_offset(1,:) = [NaN NaN];
        
        % This is NOT sped up by multithreading.
        % Flipping the 2nd frame to the first halves the number of dfts we
        % calculate.
        
        for f=2:length(frame_contenders)
            
            frm2 = double(contender_image_stack(:,:,f));
                         
            padfrm2 = padarray(frm2,[paddiffm paddiffn],0,'post'); 
            fft_ims(:,:,f) = fft2( padfrm2 );
            fft_frm2 = fft_ims(:,:,f);

            [seq_ncc(f), seq_ncc_offset(f,:), peakmaskedncc]  = auto_ref_ncc(frm1, fft_frm1, frm2, fft_frm2, numberOfOverlapPixels, xcorr_mask);
            
            if isnan(seq_ncc(f))
                group = group+1;
%                 imagesc(peakmaskedncc); axis image;
%                 pause;
            end
            
            frame_group(f) = group; % Frame 2's group.
            
            % Frame 2 is now frame 1
            fft_frm1 = fft_frm2;
            frm1 = frm2;            
        end
        
        disp(['Determined pairwise NCC. Grouping frames...']);
        
        clear frm1 frm2 padfrm1 padfrm2 fft_frm1 fft_frm2;
                
        num_groups = max(frame_group);
        seq_ncc_offset(:,1) = seq_ncc_offset(:,1)-size(contender_image_stack, 1);
        seq_ncc_offset(:,2) = seq_ncc_offset(:,2)-size(contender_image_stack, 2);
        %%
        % Find the edge frames in each group, and compare them to
        % non-adjacent groups (because we already know non-adjacents don't
        % align)
        
        for i=1:num_groups
            first_grp_inds = find(frame_group == i);
            first_grp_inds = first_grp_inds(~isnan(seq_ncc_offset(first_grp_inds,1)));
             
            first_grp_frms = double(contender_image_stack(:,:,first_grp_inds));
            first_grp_fft_frms = fft_ims(:,:,first_grp_inds);
            
            first_centered_offset = seq_ncc_offset(first_grp_inds,:);                     
            first_centered_offset = cumsum(first_centered_offset);
            first_centered_offset = first_centered_offset- ( ones(size(first_centered_offset,1),1)*mean(first_centered_offset,1) );
            
            [minoff, mininds] = min( first_centered_offset );
            [maxoff, maxinds] = max( first_centered_offset );
            
            first_extremeinds = unique([mininds, maxinds]);
            first_centered_offset = first_centered_offset(first_extremeinds,:);
            for j=1:num_groups
                second_grp_inds = find(frame_group == j);
                if j>i && ~isempty(second_grp_inds) %This prevents us from double checking an offset we already know.                     
                    second_grp_inds = second_grp_inds(~isnan(seq_ncc_offset(second_grp_inds,1)));
                    
                    second_grp_frms = double(contender_image_stack(:,:,second_grp_inds));
                    second_grp_fft_frms = fft_ims(:,:,second_grp_inds);
                    
                    second_centered_offset = seq_ncc_offset(second_grp_inds,:);                                        
                    second_centered_offset = cumsum(second_centered_offset);
                    second_centered_offset = second_centered_offset- ( ones(size(second_centered_offset,1),1)*mean(second_centered_offset,1) );

                    [minoff, mininds] = min( second_centered_offset );
                    [maxoff, maxinds] = max( second_centered_offset );
                    
                    second_extremeinds = unique([mininds, maxinds]);
                    second_centered_offset = second_centered_offset(second_extremeinds,:);
                    
                    ncc_ind_offset = zeros(length(first_extremeinds)*length(second_extremeinds),2);
                    
                    for f=1:length(first_extremeinds)                        
                        for s=1:length(second_extremeinds)
                            ind = (f-1)*length(second_extremeinds)+s;
                            [this_ncc, ncc_ind_offset(ind ,:), peakmaskedncc]  = auto_ref_ncc(first_grp_frms(:,:,first_extremeinds(f)), first_grp_fft_frms(:,:,first_extremeinds(f)),... 
                                                                                              second_grp_frms(:,:,second_extremeinds(s)), second_grp_fft_frms(:,:,second_extremeinds(s)),...
                                                                                              numberOfOverlapPixels, xcorr_mask);
                              ncc_ind_offset(ind,1) = ncc_ind_offset(ind,1)-size(contender_image_stack, 1);
                              ncc_ind_offset(ind,2) = ncc_ind_offset(ind,2)-size(contender_image_stack, 2);
                                                            
%                               if any(~isnan(ncc_ind_offset(s,:)))
%                                   ncc_ind_offset
%                               figure(1); imagesc(peakmaskedncc); axis image;
%                               transim = imtranslate(second_grp_frms(:,:,second_extremeinds(s)),[ncc_ind_offset(ind,2), ncc_ind_offset(ind,1)],'FillValues',0);
%                               figure(s+1); imshowpair(transim,first_grp_frms(:,:,first_extremeinds(f)));
%                               pause;
%                               end
                            % Make all alignments relative to the first
                            % group.
                            ncc_ind_offset(ind,:) = second_centered_offset(s,:)-ncc_ind_offset(ind,:)-first_centered_offset(f,:);
                            
                        end
                    end                                                          
                    % Remove the NaNs (failures to align).
                    ncc_ind_offset = ncc_ind_offset(~isnan(ncc_ind_offset(:,1)),:);                                        
                    

                    % If the offset isn't empty after that check, then
                    % find the translations closest together, take their
                    % median, and set this group's indices to the first group.
                    if any(~isempty(ncc_ind_offset)) && all(median(ncc_ind_offset) < (size(first_grp_frms(:,:,1))/3))
                        
                        disp(['Within distance cutoff threshold: Reassigning group ' num2str(j) ' to ' num2str(i) '...']);
                        frame_group(frame_group == j) = i;
                    end
                    
                end
            end            
        end

        
        groups = unique(frame_group);
        groups_b4 = groups;
        to_remove = false(length(frame_group),1);
        for g=1:length(groups)
            to_remove(groups(g) == frame_group) = (sum(groups(g) == frame_group)<NUM_FRM_PER_GROUP);
        end
        
        
        % If the images are part of a group that is too small, then remove
        % them.
        frame_group = frame_group(~to_remove);
        frame_contenders = frame_contenders(~to_remove);
        contender_image_stack= contender_image_stack(:,:,~to_remove);
%         fft_ims= fft_ims(:,:,~to_remove);
        clear fft_ims;
        seq_ncc = seq_ncc(~to_remove);
        seq_ncc_offset = seq_ncc_offset(~to_remove,:);     
        
        groups = unique(frame_group);
        
        fprintf('Removed %d groups from contention because they contained less than %d frames.\n', ...
                 length(groups_b4)-length(groups), NUM_FRM_PER_GROUP)
        
        %% Filter by neighbor NCC

        % When calculating the ncc threshold for a given group, only use sequential frames (further separated in time isn't fair).
        sequential_frames = [true, diff(frame_contenders) == 1]';
        ncc_threshold = median(seq_ncc(sequential_frames & ~isnan(seq_ncc) ));

%         histogram(seq_ncc,20); hold on; plot([ncc_threshold ncc_threshold],[0 10],'r'); hold off;

        rem_voting = nan(size(frame_contenders));
        voting_capacity = zeros(size(frame_contenders));  

        % Move a sliding window along the ncc values and determine which
        % frames have poor NCC with their neighbors.
        % Reminder: Each NCC value is a comparison between two frames
        for f=1:length(seq_ncc)
            if (f-1 ~= 0) && (frame_contenders(f)-frame_contenders(f-1) == 1) %If the frames are sequential (are separated by 1 frame), compare them.

                if isnan(seq_ncc(f-1))% (NaN is always beneath threshold)
                    below_threshold_vote = true;
                else
                    below_threshold_vote = seq_ncc(f-1) < ncc_threshold;
                end
                
                if isnan(rem_voting(f-1))
                    rem_voting(f-1) = 0;
                end
                if isnan(rem_voting(f))
                    rem_voting(f) = 0; 
                end

                voting_capacity(f-1)=voting_capacity(f-1)+1;
                voting_capacity(f)= voting_capacity(f)+1;
                % Add weighted votes between the frame of interest and previous frame                
                rem_voting(f-1) = rem_voting(f-1) + below_threshold_vote;
                rem_voting(f) = rem_voting(f) + below_threshold_vote;
            else
%                 frame_contenders(f)
            end

            if (f+1 < length(frame_contenders)) && (frame_contenders(f+1)-frame_contenders(f) == 1) % Between f and f+1
                
                if isnan(seq_ncc(f))% (NaN is always beneath threshold)
                    below_threshold_vote = true;
                else
                    below_threshold_vote = seq_ncc(f) < ncc_threshold;
                end

                if isnan(rem_voting(f))
                    rem_voting(f) = 0; 
                end
                if isnan(rem_voting(f+1))
                    rem_voting(f+1) = 0; 
                end

                voting_capacity(f) = voting_capacity(f)+1;
                voting_capacity(f+1)= voting_capacity(f+1)+1;
                % Add weighted votes between the frame of interest and next frame
                rem_voting(f) = rem_voting(f) + below_threshold_vote;
                rem_voting(f+1) = rem_voting(f+1) + below_threshold_vote; 
            else
%                 frame_conttenders(f)
            end            
        end

        rem_voting(isnan(rem_voting)) =0;
        
        to_remove = rem_voting == voting_capacity;
        to_retain = true(length(to_remove)-1,1);
        for v=2 :length(to_remove)
            if to_remove(v-1) || to_remove(v) %If there's a vote of no confidence in either frame involved in making the decision, drop it.
                to_retain(v-1) = false;
            end
        end

        % If you recieved all the votes you could to get removed, then you get dropped.
        frame_group = frame_group(to_retain);
        frame_contenders = frame_contenders(to_retain);
        contender_image_stack= contender_image_stack(:,:,to_retain);
%         fft_ims= fft_ims(:,:,to_retain);
        seq_ncc = seq_ncc(to_retain);
        seq_ncc_offset = seq_ncc_offset(to_retain,:);            

        %%
        % Do some final cutoffs until we have only a few possible
        % reference frames.
        % Only take frames with a minimal amount of movement between frames,
        % if possible.
        average_offset = mean(abs( seq_ncc_offset(~isnan(seq_ncc_offset(:,1)),:) ));
        std_offset = std(abs(seq_ncc_offset(~isnan(seq_ncc_offset(:,1)),:) ));
        keep_list = [];
        for g=1:length(groups)
            
            groupind = find(groups(g) == frame_group);
            
            % If the absolute average offset is below the global threshold,
            % then include the frames.
            to_retain = sum( abs(seq_ncc_offset(groupind,:)) < repmat(average_offset, size(seq_ncc_offset(groupind,:),1),1), 2) == 2;
            
            % If by doing this we completely eradicate this group,
            % reconsider the removal by relaxing the offset to 1 std over
            % the mean.
            if sum(to_retain) == 0
                to_retain = sum( abs(seq_ncc_offset(groupind,:)) < repmat(average_offset+std_offset, size(seq_ncc_offset(groupind,:),1),1), 2) == 2;
                
                % If we still don't have anything, then remove the frames
                % from consideration anyway as they're probably unreliable.                                
            end
            % Keep track of which indexes are worth keeping.
            keep_list = [keep_list; groupind(to_retain)];
        end
        %%
        
        frame_group = frame_group(keep_list);
        frame_contenders = frame_contenders(keep_list);
        contender_image_stack= contender_image_stack(:,:,keep_list);
%         fft_ims= fft_ims(:,:,keep_list);
        seq_ncc = seq_ncc(keep_list);
        seq_ncc_offset = seq_ncc_offset(keep_list,:);   

        groups = unique(frame_group);
        
        [seq_ncc, sortinds] =sort(seq_ncc,1,'descend');
        
        % Sort the images based on their sequential NCC
        frame_group = frame_group(sortinds);
        frame_contenders = frame_contenders(sortinds);
        contender_image_stack= contender_image_stack(:,:,sortinds);
%         fft_ims= fft_ims(:,:,sortinds);
        
        seq_ncc_offset = seq_ncc_offset(sortinds,:);
        
        ref_size =0;
        for g=1:length(groups)
            groupind = find(groups(g) == frame_group);
            ref_size = max([length(groupind) ref_size]);
        end
        
        reference_frames = -ones(ref_size, length(groups));
        
        for g=1:length(groups)
            groupind = find(groups(g) == frame_group);
            
%             vidObj = VideoWriter(['Group_' num2str(g) '.avi']);
%             open(vidObj);
            
            for v=1:length(groupind)
%                 writeVideo(vidObj, image_stack(:,:,frame_contenders(groupind(v))) );
                reference_frames(v,g) = frame_contenders(groupind(v));
            end
            
%             close(vidObj);
        end
        
    else
        % NOP for now.
    end

end

