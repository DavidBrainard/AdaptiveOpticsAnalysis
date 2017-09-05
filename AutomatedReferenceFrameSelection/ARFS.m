% Rob Cooper 06-30-2017
%
% This script is an implementation of the algorithm outlined by Salmon et
% al 2017: "An Automated Reference Frame Selection (ARFS) Algorithm for 
% Cone Imaging with Adaptive Optics Scanning Light Ophthalmoscopy"
%

clear;
close all;

STRIP_SIZE = 40;

locInd=1; % Location index

% For Debug
mov_path = {pwd}; %,...
            %pwd,...
            %pwd};
stack_fname = {'NC_11028_20160601_OD_confocal_0011_desinusoided.avi'};
% {'NC_11049_20170629_confocal_OD_0000_0003_desinusoided.avi',...
%                'NC_11049_20170629_avg_OD_0000_desinusoided.avi',...
%                'NC_11049_20170629_split_det_OD_0000_desinusoided.avi'};

for modalityInd=1 : 1%length(stack_fname) 
    vidReader = VideoReader( fullfile(mov_path{locInd,modalityInd}, stack_fname{locInd, modalityInd}) );
    
    i=1;
    while(hasFrame(vidReader))
        image_stack(:,:,i,modalityInd) = uint8(readFrame(vidReader));        
        frame_mean(i,modalityInd) = mean2(image_stack(:,:,i,modalityInd));
        i = i+1;
    end
    numFrames = i-1;
    
    % Get some basic heuristics from each modality.
    mode_mean(modalityInd) = mean(frame_mean(:,modalityInd));
    mode_dev(modalityInd) = std(frame_mean(:,modalityInd));

    frame_contenders(:,modalityInd) = (1:numFrames);
    

    strip_inds = 0:STRIP_SIZE:size(image_stack(:,:,1, 1),2);
    strip_inds(1) = 1;
    if strip_inds(end) ~= size(image_stack,2)
        if (size(image_stack,2)-strip_inds(end)) > STRIP_SIZE/2
            strip_inds = [strip_inds size(image_stack,2)];
        else
            strip_inds(end) = size(image_stack,2);
        end
    end
    num_strips = length(strip_inds)-1;
    
    %% Filter by image mean
    mean_contenders = false(1,numFrames);
    for f=1:numFrames        
        mean_contenders(f) =  (frame_mean(f,modalityInd) < mode_mean(modalityInd)+2*mode_dev(modalityInd)) &&...
                              (frame_mean(f,modalityInd) > mode_mean(modalityInd)-2*mode_dev(modalityInd));
    end
    
    frame_contenders = frame_contenders(mean_contenders,modalityInd);
    
    numFrames = length(frame_contenders);
    
    radon_bandwidth = zeros(numFrames,num_strips);
    %%
    
    if exist('parfor','builtin') == 5 % If we can multithread it, do it!
        %% Filter by Radon transform FWHM
        parfor f=1:numFrames
            
            frame_ind=frame_contenders(f,modalityInd);
            
            for s=1:num_strips
    
                % Get the log power spectrum for us to play with
                pwr_spect = ( abs(fftshift(fft2(image_stack(strip_inds(s):strip_inds(s+1),:, frame_ind),512, 512))).^2);
                % From our padding, the center vertical frequency will be
                % garbage- remove it for our purposes.
                pwr_spect = log10(pwr_spect(:,[1:256 258:512]));
                
                % Threshold is set using the upper 2 std devs
                thresh_pwr_spect = ( pwr_spect>(mean(pwr_spect(:))+2*std(pwr_spect(:))) );
                
                radoned = radon( thresh_pwr_spect );                

                % Determine the minimum and maximum FWHM
%                 tic;
                halfmax = repmat(max(radoned)./2,[727 1]);
%                 toc;
                fwhm = sum(radoned>halfmax);
                
                % 
                radon_bandwidth(f,s) = max(fwhm)-min(fwhm);

            end
        end        
        
        threshold = mean(radon_bandwidth(:))+ 2*std(radon_bandwidth(:));
        % After thresholding and removal, update the contenders list.
        frame_contenders = frame_contenders(~any(radon_bandwidth > threshold,2));        
        
        contender_image_stack = image_stack(:,:, frame_contenders(:, modalityInd), modalityInd);
        
        frm1 = double( contender_image_stack(:,:,1) );
        [m, n] = size(frm1);
        paddiffm = (size(frm1,1)*2)-1-m;
        paddiffn = (size(frm1,2)*2)-1-n;

        % Make a mask to remove any edge effects.
        [maskdistx, maskdisty] = meshgrid( 1:(size(frm1,2)*2)-1, 1:(size(frm1,1)*2)-1);
        
        maskdistx = maskdistx-(size(maskdistx,2)/2)-STRIP_SIZE;
        maskdisty = maskdisty-(size(maskdistx,1)/2)-STRIP_SIZE;
        
        xcorr_mask = sqrt(maskdistx.^2 + maskdisty.^2) <400;
        % Determine the number of pixels that will overlap between the two
        % images at any given point.
        numberOfOverlapPixels = arfs_local_sum(ones(size(frm1)),size(frm1,1),size(frm1,2));
        
        padfrm1 = padarray(frm1,[paddiffm paddiffn],0,'post');
        
        fft_frm1 =  fft2( padfrm1 );
        
        ncc = zeros(length(frame_contenders),1);
        ncc_offset = zeros(length(frame_contenders)-1,2);
        frame_group = zeros(length(frame_contenders)-1,1);
        fft_ims = zeros(size(padfrm1,1), size(padfrm1,2), length(frame_contenders));
        
        fft_ims(:,:,1) = fft2( padfrm1 );
        
        group = 1;        
        frame_group(1) = group;
        ncc(1) = NaN; % 0 to 1 doesn't align.
        ncc_offset(1,:) = [NaN NaN];
        
        % This is NOT sped up by multithreading.
        % Flipping the 2nd frame to the first halves the number of dfts we
        % calculate.
        tic;
        for f=2:length(frame_contenders)
            
            frm2 = double(contender_image_stack(:,:,f));
                         
            padfrm2 = padarray(frm2,[paddiffm paddiffn],0,'post'); 
            fft_ims(:,:,f) = fft2( padfrm2 );
            fft_frm2 = fft_ims(:,:,f);

            [ncc(f), ncc_offset(f,:), peakmaskedncc]  = auto_ref_ncc(frm1, fft_frm1, frm2, fft_frm2, numberOfOverlapPixels, xcorr_mask);
            
            if isnan(ncc(f))
                group = group+1;
            end
            
            frame_group(f) = group; % Frame 2's group.
            
            % Frame 2 is now frame 1
            fft_frm1 = fft_frm2;
            frm1 = frm2;            
        end
        toc;
                
        num_groups = max(frame_group);
        ncc_offset(:,1) = ncc_offset(:,1)-size(contender_image_stack, 1);
        ncc_offset(:,2) = ncc_offset(:,2)-size(contender_image_stack, 2);
        %%
        % Find the edge frames in each group, and compare them to
        % non-adjacent groups (because we already know non-adjacents don't
        % align)
        for i=1:num_groups
            first_grp_inds = find(frame_group == i);
            first_grp_inds = first_grp_inds(~isnan(ncc_offset(first_grp_inds,1)));
             
            first_grp_frms = double(contender_image_stack(:,:,first_grp_inds));
            first_grp_fft_frms = fft_ims(:,:,first_grp_inds);
            
            first_centered_offset = ncc_offset(first_grp_inds,:);                     
            first_centered_offset = cumsum(first_centered_offset);
            first_centered_offset = first_centered_offset- ( ones(size(first_centered_offset,1),1)*mean(first_centered_offset,1) );
            
            [minoff, mininds] = min( first_centered_offset );
            [maxoff, maxinds] = max( first_centered_offset );
            
            first_extremeinds = unique([mininds, maxinds]);
            first_centered_offset = first_centered_offset(first_extremeinds,:);
            for j=1:num_groups
                if j>i %This prevents us from double checking an offset we already know. 
                    second_grp_inds = find(frame_group == j);
                    second_grp_inds = second_grp_inds(~isnan(ncc_offset(second_grp_inds,1)));
                    
                    second_grp_frms = double(contender_image_stack(:,:,second_grp_inds));
                    second_grp_fft_frms = fft_ims(:,:,second_grp_inds);
                    
                    second_centered_offset = ncc_offset(second_grp_inds,:);                                        
                    second_centered_offset = cumsum(second_centered_offset);
                    second_centered_offset = second_centered_offset- ( ones(size(second_centered_offset,1),1)*mean(second_centered_offset,1) );

                    [minoff, mininds] = min( second_centered_offset );
                    [maxoff, maxinds] = max( second_centered_offset );
                    
                    second_extremeinds = unique([mininds, maxinds]);
                    second_centered_offset = second_centered_offset(second_extremeinds,:);
                    
                    ncc_ind_offset = zeros(length(first_extremeinds)*length(second_extremeinds),2);
                    
                    for f=1:length(first_extremeinds)                        
                        for s=1:length(second_extremeinds)
                            ind = ((f-1)*(length(first_extremeinds))+s);
                            [this_ncc, ncc_ind_offset(ind ,:), peakmaskedncc]  = auto_ref_ncc(first_grp_frms(:,:,first_extremeinds(f)), first_grp_fft_frms(:,:,first_extremeinds(f)),... 
                                                                                      second_grp_frms(:,:,second_extremeinds(s)), second_grp_fft_frms(:,:,second_extremeinds(s)),...
                                                                                      numberOfOverlapPixels, xcorr_mask);
                              ncc_ind_offset(ind,1) = ncc_ind_offset(ind,1)-size(contender_image_stack, 1);
                              ncc_ind_offset(ind,2) = ncc_ind_offset(ind,2)-size(contender_image_stack, 2);
                              
                              
%                               if any(~isnan(ncc_ind_offset(s,:)))
%                               figure(1); imagesc(peakmaskedncc); axis image;
%                               transim = imtranslate(second_grp_frms(:,:,second_extremeinds(s)),[ncc_ind_offset(s,2), ncc_ind_offset(s,1)],'FillValues',0);
%                               figure(s+1); imshowpair(transim,first_grp_frms(:,:,first_extremeinds(f)));
%                               end
                            % Make all alignments relative to the first
                            % group.
                            ncc_ind_offset(ind,:) = second_centered_offset(s,:)-ncc_ind_offset(ind,:)-first_centered_offset(f,:);
                        end
                        ncc_ind_offset
                        
                    end
                    
                    pause;
                    
                    % Remove the NaNs (failures to align).
                    ncc_ind_offset = ncc_ind_offset(~isnan(ncc_ind_offset(:,1)),:);

                    % If the offset isn't empty after that check, then
                    % find the translations closest together, average
                    % them, and add this group's indices to the first group.
                    if ~isempty(ncc_ind_offset)
                        
                        
                    end
                    
                end
            end
        end
        %%
%         for f=2:length(ncc_offset)
%             if ~isnan(ncc_offset(f,1)) 
%             transim = imtranslate(contender_image_stack(:,:,f),[ncc_offset(f,2), ncc_offset(f,1)],'FillValues',0);        
%             figure(2);imshowpair(transim,contender_image_stack(:,:,f-1));
%             title(num2str([ncc_offset(f,2), ncc_offset(f,1)]));
%             end
%             
%         end
        
        
        %% Filter by neighbor NCC
        % When calculating the ncc threshold, only use sequential frames (further separated in time isn't fair).
        sequential_frames = [true; diff( frame_contenders ) == 1];
        ncc_threshold = median(ncc(sequential_frames & ~isnan(ncc) ));
        
        seq_ncc = ncc(sequential_frames);
        
        hist(seq_ncc,20); hold on; plot([ncc_threshold ncc_threshold],[0 10],'r'); hold off;
                
         
        rem_voting = nan(size(frame_contenders));
        voting_capacity = zeros(size(frame_contenders));               
        % Move a sliding window along the ncc values and determine which
        % frames have poor NCC with their neighbors.
        % Reminder: Each NCC value is a comparison between two frames
        for f=1:length(ncc)
            if (f-1 ~= 0) && (frame_contenders(f)-frame_contenders(f-1) == 1) % Between f-1 and f
                
                atthresh = ncc(f-1) < ncc_threshold;
                
                if isnan(rem_voting(f-1))
                    rem_voting(f-1) = 0;
                end
                if isnan(rem_voting(f))
                    rem_voting(f) = 0; 
                end
                
                voting_capacity(f-1)=voting_capacity(f-1)+1;
                voting_capacity(f)= voting_capacity(f)+1;
                % Add weighted votes between the frame of interest and previous frame
                % (NaN is always beneath threshold)
                rem_voting(f-1) = rem_voting(f-1) + ((ncc(f-1) < ncc_threshold) );
                rem_voting(f) = rem_voting(f) + ((ncc(f-1) < ncc_threshold) );
            else
%                 frame_contenders(f)
            end
            
            if (f+1 < length(frame_contenders)) && (frame_contenders(f+1)-frame_contenders(f) == 1) % Between f and f+1
                atthresh = ncc(f) < ncc_threshold;
                
                if isnan(rem_voting(f))
                    rem_voting(f) = 0; 
                end
                if isnan(rem_voting(f+1))
                    rem_voting(f+1) = 0; 
                end
                
                voting_capacity(f) = voting_capacity(f)+1;
                voting_capacity(f+1)= voting_capacity(f+1)+1;
                % Add weighted votes between the frame of interest and next frame
                rem_voting(f) = rem_voting(f) + (ncc(f) < ncc_threshold);
                rem_voting(f+1) = rem_voting(f+1) + (ncc(f) < ncc_threshold); 
            else
%                 frame_conttenders(f)
            end            
        end
        
        % VALIDATE THIS!
        to_remove = rem_voting == voting_capacity;
        remove_by_vote_paired = true(length(to_remove)-1,1);
        for v=2 :length(to_remove)
            if to_remove(v-1) || to_remove(v)
                remove_by_vote_paired(v-1) = false;
            end
        end
        
%         figure; plot(absolute_offsets(:,1), absolute_offsets(:,2),'.'); title('Before removal');
        
%         % If you recieved all the votes you could to get removed, then you get dropped.
%         frame_contenders = frame_contenders(rem_voting ~= voting_capacity);
%         ncc = ncc(remove_by_vote_paired);
%         ncc_offset = ncc_offset(remove_by_vote_paired,:);
%         absolute_offsets = absolute_offsets(remove_by_vote_paired,:);
        
%         figure; plot(absolute_offsets(:,1), absolute_offsets(:,2),'.'); title('After removal');
%         % If we can't be sure you're a good frame, then you get dropped.
%         voting_capacity = voting_capacity(rem_voting ~= voting_capacity);
%         frame_contenders = frame_contenders(voting_capacity ~= 2);
%         ncc = ncc(rem_voting ~= voting_capacity);
%         ncc_offset = ncc_offset(voting_capacity ~= 2,:);
        

        
        
        
        
        
        
        
        
        for f=1:length(frame_contenders)-1
            frame_ind=frame_contenders(f,modalityInd);
            frm1 = double(image_stack(:,:,frame_ind,modalityInd));
            frm2 = double(image_stack(:,:,frame_ind+1,modalityInd));
            xcored(f) = std2( ( ((frm2-mean2(frm2))./std2(frm2)) - ((frm1-mean2(frm1))./std2(frm1))) );
        end
        hist(xcored,20);
    else
        
    end
    toc;
end


    


                