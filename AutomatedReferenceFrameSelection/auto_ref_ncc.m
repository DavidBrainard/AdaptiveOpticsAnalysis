function [ ncc, ncc_offset, ncc_map ] = auto_ref_ncc( frm1, fft_frm1, frm2, fft_frm2, overlapPixelMap, xcorr_mask )

    if ~exist('overlapPixelMap','var')
        % Determine the number of pixels that will overlap between the two
        % EQUAL-SIZED images at any given point.
        overlapPixelMap = arfs_local_sum(ones(size(frm1)),size(frm1,1),size(frm1,2));
    end

    if ~exist('xcorr_mask','var')
        % Make a mask to remove any edge effects.
        [maskdistx, maskdisty] = meshgrid( 1:(size(frm1,2)*2)-1, 1:(size(frm1,1)*2)-1);

        maskdistx = maskdistx-(size(maskdistx,2)/2);
        maskdisty = maskdisty-(size(maskdistx,1)/2);

        xcorr_mask = sqrt(maskdistx.^2 + maskdisty.^2) <400;
    end

    [m, n] = size(frm1);

    % Denominator for NCC (Using the local sum method used by
    % MATLAB and J.P Lewis.            
    local_sum_F = arfs_local_sum(frm1,m,n);
    local_sum_F2 = arfs_local_sum(frm1.*frm1,m,n);

    rotfrm2 = rot90(frm2,2); % Need to rotate because we'll be conjugate in the fft (double-flipped)
    local_sum_T = arfs_local_sum(rotfrm2,m,n);
    local_sum_T2 = arfs_local_sum(rotfrm2.*rotfrm2,m,n);

    % f_uv, assumes that the two images are the same size (they should
    % be for this application)- otherwise we'd need to calculate
    % overlap, ala normxcorr2_general by Dirk Padfield
    denom_F = max( local_sum_F2- ((local_sum_F.^2)./overlapPixelMap), 0);
    % t
    denom_T = max( local_sum_T2- ((local_sum_T.^2)./overlapPixelMap), 0);

    % Numerator for NCC
    numerator = fftshift( real(ifft2( fft_frm1 .* conj( fft_frm2 ) ) ));

    numerator = numerator - local_sum_F.*local_sum_T./overlapPixelMap;

    denom = sqrt(denom_F.*denom_T);

    ncc_frm = (numerator./denom);
    ncc_frm(isnan(ncc_frm)) =0;
    %             Too slow; only use in the future if we need the performance             
    %             stddev_ncc_frm = stdfilt(ncc_frm,ones(11)).*xcorr_mask;
    ncc_map = (ncc_frm.*xcorr_mask);
    ncc_map(isnan(ncc_map)) =0;

    [ncc, ncc_ind]  = max(ncc_map(:));

    [roff, coff]= ind2sub(size(ncc_map),ncc_ind);            

    % Method from the Salmon paper
    peakmaskedncc = ncc_map;
%     peakmaskedncc( (roff-10):(roff+10),(coff-10):(coff+10)) = 0; % Mask out the found peak
    peakmaskedncc( peakmaskedncc > (ncc*(1-1/2.71)) )  = 1;
    peakmaskedncc( peakmaskedncc < (ncc*(1-1/2.71)) )  = 0;

    comps = bwconncomp(imdilate(peakmaskedncc,ones(11)));
    
    if comps.NumObjects > 1 % There's multiple peaks; so massage the data a little more until to remove noise.
               
        [horzgrad, vertgrad]=gradient(ncc_frm);
        grad = (abs(horzgrad).*abs(vertgrad)).*xcorr_mask;
        
        blurry_grad = imgaussfilt(grad,10);
        
        % Remove garbage data
        blurry_grad(isnan(blurry_grad)) = mean(blurry_grad(~isnan(blurry_grad)));
        blurry_grad(isinf(blurry_grad)) = mean(blurry_grad(~isinf(blurry_grad)));
        
        [blurryncc, blurry_ncc_ind]  = max(blurry_grad(:));
        
        [roff, coff]= ind2sub(size(ncc_map),blurry_ncc_ind);  
        
        
        blurry_grad( blurry_grad > (blurryncc*(1-1/2.71)) )  = 1;
        blurry_grad( blurry_grad < (blurryncc*(1-1/2.71)) )  = 0;
        % After thresholding, should only be one component.
        comps = bwconncomp(blurry_grad);
        
        if comps.NumObjects > 1 % They don't align; make sure we keep track of that.
            ncc = NaN;
            ncc_offset = [NaN NaN];
        else
            ncc = ncc_map(blurry_ncc_ind);
            ncc_offset = [roff, coff];
        end
    else                
        ncc_offset = [roff, coff];
    end

end

