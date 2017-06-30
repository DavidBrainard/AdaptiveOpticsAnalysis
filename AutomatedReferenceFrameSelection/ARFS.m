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
mov_path = {'C:\Users\Gabriel Cavanagh\Box Sync\General_Postdoc_Work\Ref_Frame_Finder', 'C:\Users\Gabriel Cavanagh\Box Sync\General_Postdoc_Work\Ref_Frame_Finder', 'C:\Users\Gabriel Cavanagh\Box Sync\General_Postdoc_Work\Ref_Frame_Finder'};
stack_fname = {'NC_11049_20170629_confocal_OD_0000_desinusoided.avi', 'NC_11049_20170629_avg_OD_0000_desinusoided.avi', 'NC_11049_20170629_split_det_OD_0000_desinusoided.avi'};

for modalityInd=1 : 1%length(stack_fname) 
    vidReader = VideoReader( fullfile(mov_path{locInd,modalityInd}, stack_fname{locInd, modalityInd}) );
    
    i=1;
    while(hasFrame(vidReader))
        image_stack(:,:,i,modalityInd) = uint8(readFrame(vidReader));
        i = i+1;
    end
    numFrames = i-1;


    frame_contenders(:,modalityInd) = 1:i;
    vectorablestk = double(image_stack(:,:,:,modalityInd));

    % Get some basic heuristics from each modality
    mode_mean(modalityInd) = mean(vectorablestk(:));
    mode_dev(modalityInd) = std(vectorablestk(:));

    strip_inds = 0:STRIP_SIZE:size(image_stack(:,:,1, 1),2);
    strip_inds(1) = 1;
    if strip_inds(end) ~= size(image_stack(1:40,:,93,modalityInd),2)
        strip_inds = [strip_inds size(image_stack(:,:,1,modalityInd),2)];
    end
    num_strips = length(strip_inds)-1;
    
    % Debug examples
        disted=image_stack(1:40,:,93,modalityInd);
        undisted=image_stack(100:140,:,93,modalityInd);        
        
    centeraffinityfunc = sqrt(abs( (1:727)-round(727/2)));
    centeraffinityfunc = 1-(centeraffinityfunc./max(centeraffinityfunc));/.
    centeraffinityfunc = repmat(centeraffinityfunc',[1 180]);
    
%     houghs = zeros(1445,180,150,num_strips);
    strip_affinity = zeros(150,num_strips);
    %%
    tic;
    if exist('parfor','builtin') == 5 % If we can multithread it, do it!
        parfor f=1:numFrames
            
            for s=1:num_strips
    
                % Get the log power spectrum for us to play with
                pwr_spect = ( abs(fftshift(fft2(image_stack(strip_inds(s):strip_inds(s+1),:, f),512, 512))).^2);
                % From our padding, the center vertical frequency will be
                % garbage- remove it for our purposes.
                pwr_spect = log10(pwr_spect(:,[1:256 258:512]));
                
                % Threshold is set using the upper 2 std devs
                thresh_pwr_spect = ( pwr_spect>(mean(pwr_spect(:))+2*std(pwr_spect(:))) );
                
                radoned = radon( thresh_pwr_spect );                
                figure; imagesc(radoned);
                
                % Determine the minimum and maximum FWHM
                tic;
                halfmax = repmat(max(radoned)./2,[727 1]);
                toc;
                fwhm = sum(radoned>halfmax);
                
                % Multiply the energies by the radon xform to get the
                % average central affinity
%                 affinity=radoned.*centeraffinityfunc;
                strip_affinity(f,s) = max(fwhm)-min(fwhm);
                         
            end 
        end
        
    else
        
    end
    toc;
end

                

                