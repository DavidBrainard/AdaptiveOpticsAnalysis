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
mov_path = {pwd,...
            pwd,...
            pwd};
stack_fname = {'NC_11049_20170629_confocal_OD_0000_0003_desinusoided.avi',...
               'NC_11049_20170629_avg_OD_0000_desinusoided.avi',...
               'NC_11049_20170629_split_det_OD_0000_desinusoided.avi'};

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
    if strip_inds(end) ~= size(image_stack(1:40,:,93,modalityInd),2)
        strip_inds = [strip_inds size(image_stack(:,:,1,modalityInd),2)];
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
        
        maskdistx = maskdistx-(size(maskdistx,2)/2);
        maskdisty = maskdisty-(size(maskdistx,1)/2);
        
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
        ncc_offset(:,1) = ncc_offset(:,1)-size(contender_image_stack(:,:,1), 1);
        ncc_offset(:,2) = ncc_offset(:,2)-size(contender_image_stack(:,:,1), 2);
        
        % Find the center-most frame in each group, and compare them to
        % non-adjacent groups (because we already know non-adjacents don't
        % align)
        for i=1:num_groups
            first_grp_inds = (frame_group == i);
            first_grp_frms = contender_image_stack(:,:,first_grp_inds);
            first_grp_fft_frms = fft_ims(:,:,first_grp_inds);
            
            first_group_offset = ncc_offset(first_grp_inds,:);            
            first_group_offset(isnan(first_group_offset(:,1)),:) = [];
            
            first_group_offset = cumsum(first_group_offset);
            
            [minoff, mininds] = min( first_group_offset );
            [maxoff, maxinds] = max( first_group_offset );
            
            extremeinds = unique([mininds, maxinds]);
            
            for j=1:numgroups
                if i~=j
                    second_grp_inds = (frame_group == i);
                    second_grp_frms = fft_ims(:,:,second_grp_inds);
                    second_grp_fft_frms = fft_ims(:,:,second_grp_inds);
                    
                    [this_ncc, ncc_ind_offset, peakmaskedncc]  = auto_ref_ncc(frm1, fft_frm1, frm2, fft_frm2, numberOfOverlapPixels, xcorr_mask);
                    
                end
            end
        end
        
        
        %% Filter by neighbor NCC
        % When calculating the ncc threshold, only use sequential frames (further separated in time isn't fair).
        sequential_frames = diff( frame_contenders ) == 1;
        ncc_threshold = median(ncc(sequential_frames & ~isnan(ncc) ));
        
        seq_ncc = ncc(sequential_frames);
        
        hist(seq_ncc,20); hold on; plot([ncc_threshold ncc_threshold],[0 10],'r'); hold off;
                
         
        
%         absolute_offsets = cumsum(ncc_offset);
        
        rem_voting = nan(size(frame_contenders));
        voting_capacity = zeros(size(frame_contenders));               
        % Move a sliding window along the ncc values and determine which
        % frames consistently have poor NCC with their surrounding
        % neighbors.
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
        

        
        
        for f=2:length(ncc_offset)            
            transim = imtranslate(contender_image_stack(:,:,f),[ncc_offset(f-1,2), ncc_offset(f-1,1)],'FillValues',0);        
            figure(2);imshowpair(transim,contender_image_stack(:,:,f-1));
            title(num2str(this_ncc(f-1)));
            pause;
        end
        
        
        
        
        
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


    


                