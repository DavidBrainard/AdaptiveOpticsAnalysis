function [ sdImageMasked  ] = tam_etal_capillary_func( pceptVideo )
% Robert F Cooper
% 01-08-2016
%
% This function incorporates an algorithm based on Johnny Tam's motion 
% contrast enhancements for an input avi dataset filename.

numFrames = size(pceptVideo,3);

for i=1:numFrames
    pceptVideo(:,:,i) = medfilt2(pceptVideo(:,:,i), [5,5]);
end

divVideo = zeros(size(pceptVideo,1),size(pceptVideo,2),numFrames-1);

% Create the division video
% figure(1); title('Division video');
for i=1:numFrames-1
    
    % Find the regions in each frame that do not have data, make a mask to
    % ignore in all frames
    ignoremask = (pceptVideo(:,:,i) == 0) | (pceptVideo(:,:,i+1) ==0);
    
    divVideo(:,:,i) = ~ignoremask.*(pceptVideo(:,:,i)./pceptVideo(:,:,i+1));
    
    im = divVideo(:,:,i);

    im(isnan(im)) = 0;
    
    divVideo(:,:,i) = im;

%     imagesc( divVideo(:,:,i) ); axis image; colormap gray;
%     pause(0.2);
end

mfDivVideo = zeros(size(pceptVideo,1),size(pceptVideo,2),size(divVideo,3)-1);

for i=1:size(divVideo,3)-1
    
    % Find the regions in each frame that do not have data, make a mask to
    % ignore in all frames
    ignoremask = (divVideo(:,:,i) == 0) | (divVideo(:,:,i+1) ==0);
    ignoremask = imdilate(ignoremask, strel('disk',5) );
    
    mfDivVideo(:,:,i) = ~ignoremask.*(divVideo(:,:,i)+divVideo(:,:,i+1))/2 ;
    
    im = mfDivVideo(:,:,i);
    
    im(isnan(im)) = 0;
    
    mfDivVideo(:,:,i) = im;
    
end

% meanMfDiv = mean(mfDivVideo,3); % Arithmetic std deviation

% Geometric std deviation
% meanMfDiv = nthroot(prod(fDivProdVideo,3), size(mfDivVideo,3)); 
% Make everything a 1 so that it does not contribute to the multiplication
mfDivProdIm = ones(size(pceptVideo,1),size(pceptVideo,2));
for i=1:size(mfDivVideo,3)
    
    im = mfDivVideo(:,:,i);
    im(im==0) = 1;
    
    mfDivProdIm = mfDivProdIm.*im;
end
meanMfDiv = nthroot(mfDivProdIm, size(mfDivVideo,3));
logMeanMfDiv = log(meanMfDiv);

sdImage = zeros(size(pceptVideo,1),size(pceptVideo,2));
mfDivMask = zeros(size(pceptVideo,1),size(pceptVideo,2));

for i=1:size(mfDivVideo,3)
   
%     sdImage = sdImage + (mfDivVideo(:,:,i)-meanMfDiv).^2;

    mfDivImage = mfDivVideo(:,:,i);
    %Only take the log over the areas that are not 0.
    
    mfDivMask = mfDivMask + 1*(mfDivImage ~= 0);
    mfDivImage( mfDivImage ~= 0 ) = log(mfDivImage( mfDivImage ~= 0 ));

    sdImage = sdImage + (mfDivImage -logMeanMfDiv ).^2;
%     figure(1); imagesc(sdImage); colormap gray;
%     figure(2); imagesc(mfDivImage); colormap gray;
end

% mfDivMask( mfDivMask < 0.33*numFrames ) = 1;

sdImage = exp( sqrt( sdImage./(mfDivMask-1) ) );
% sdImage = sdImage./(size(mfDivVideo,3)-1);

% Squish undefined parts of the image by making them the minimum value
undefsd = (isnan(sdImage) | isinf(sdImage) | (sdImage ==0));


sdImage( undefsd ) = min( sdImage( ~undefsd ) );

% Gaussian filter the shit out of it because we a) don't have many frames
% to play with and b) don't care about the background noise
sdImage = imgaussfilt(sdImage,10,'FilterSize',51);

% Linearly contrast stretch the image
sdImagemin = min(sdImage(:));
sdImageminsub = sdImage-sdImagemin;
sdImagestretched = 255*sdImageminsub./max(sdImageminsub(:));

% figure(1); imagesc(adapthisteq(sdImagestretched,'ClipLimit',0.005)); colormap gray; axis image; title('SD Image');

% figure(1); imagesc( sdImagestretched ); colormap gray; axis image; title('SD Image');

% imwrite(uint8(sdImagestretched),parula(256),'NC_11049_20160209_OD_confocal_0003_ref_65_affine_crop_TAM_capillaries.tif');

threshold = mean(sdImagestretched(:))  + std(sdImagestretched(:))/3;
% figure(2); imagesc( sdImagestretched > threshold ); colormap gray; axis image; title('Vessel Mask');

sdImageMasked = (sdImagestretched > threshold);

end

