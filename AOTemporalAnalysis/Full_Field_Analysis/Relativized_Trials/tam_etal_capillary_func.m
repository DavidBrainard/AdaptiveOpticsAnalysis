function [ sdImageMasked, varImage  ] = tam_etal_capillary_func( pceptVideo )
% Robert F Cooper
% 01-08-2016
%
% This function incorporates an algorithm based on Johnny Tam's motion 
% contrast enhancements for an input avi dataset filename.
%
% @param pceptVideo:
%
% A 3D stack containing the video frames from which the motion contrast will
% be generated

numFrames = size(pceptVideo,3);
ignoremask = zeros(size(pceptVideo));

for i=1:numFrames-1  
    % Find the regions in each frame that do not have data, make a mask to
    % ignore in all frames
    ignoremask(:,:,i) = imdilate( (pceptVideo(:,:,i) == 0) | (pceptVideo(:,:,i+1) ==0), strel('disk',5));
end

for i=1:numFrames    
%     pceptVideo(:,:,i) = medfilt2(pceptVideo(:,:,i), [5,5]);
    pceptVideo(:,:,i) = imgaussfilt(pceptVideo(:,:,i),6,'FilterSize',31);     
end

divVideo = zeros(size(pceptVideo,1),size(pceptVideo,2),numFrames-1);

% Create the division video
% figure(1); title('Division video');
for i=1:numFrames-1
    

    divVideo(:,:,i) = ~ignoremask(:,:,i).*(pceptVideo(:,:,i)./pceptVideo(:,:,i+1));
    
    im = divVideo(:,:,i);

    im(isnan(im)) = 0;
    
    divVideo(:,:,i) = im;

%     imagesc( divVideo(:,:,i) ); axis image; colormap gray; title(num2str(i));
%     drawnow;
%     pause(1);
end
clear pceptVideo;

mfDivVideo = zeros(size(divVideo,1), size(divVideo,2), size(divVideo,3)-1);

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
clear divVideo;

% Geometric std deviation
% meanMfDiv = nthroot(prod(fDivProdVideo,3), size(mfDivVideo,3)); 
% Make everything a 1 so that it does not contribute to the multiplication
mfDivProdIm = ones(size(mfDivVideo,1),size(mfDivVideo,2));
for i=1:size(mfDivVideo,3)
    
    im = mfDivVideo(:,:,i);
    im(im==0) = 1;
    
    mfDivProdIm = mfDivProdIm.*im;
end
meanMfDiv = nthroot(mfDivProdIm, size(mfDivVideo,3));
logMeanMfDiv = log(meanMfDiv);

varImage = zeros(size(mfDivVideo,1),size(mfDivVideo,2));
mfDivMask = zeros(size(mfDivVideo,1),size(mfDivVideo,2));

for i=1:size(mfDivVideo,3)
   
%     sdImage = sdImage + (mfDivVideo(:,:,i)-meanMfDiv).^2;

    mfDivImage = mfDivVideo(:,:,i);
    %Only take the log over the areas that are not 0.
    
    mfDivMask = mfDivMask + 1*(mfDivImage ~= 0);
    mfDivImage( mfDivImage ~= 0 ) = log(mfDivImage( mfDivImage ~= 0 ));

    varImage = varImage + (mfDivImage -logMeanMfDiv ).^2;
%     figure(1); imagesc(sdImage); colormap gray;
%     figure(2); imagesc(mfDivImage); colormap gray;
end

% mfDivMask( mfDivMask < 0.33*numFrames ) = 1;

sdImage = exp( sqrt( varImage./(mfDivMask-1) ) );

varImage = varImage./(mfDivMask-1);

% Squish undefined parts of the image by making them 0
undefsd = (isnan(sdImage) | isinf(sdImage) );

sdImage( undefsd ) = 1;

meaningfulregion = sdImage~=1;

% edgemask = 20;
% varImage(1:edgemask,:)=min(varImage(:));
% varImage(:,1:edgemask)=min(varImage(:));
% varImage(end-edgemask+1:end,:)=min(varImage(:));
% varImage(:,end-edgemask+1:end)=min(varImage(:));

% Gaussian filter the shit out of it because we a) don't have many frames
% to play with and b) don't care about the background noise
sdImage = imgaussfilt(sdImage,10,'FilterSize',51);

% figure(1); imagesc(varImage); colormap gray;

% centerVarImage = varImage(edgemask+1:end-edgemask+1,edgemask+1:end-edgemask+1);

% Linearly contrast stretch the non-zero parts of the image.
sdImagemin = min(sdImage(meaningfulregion));
sdImageminsub = sdImage-sdImagemin;
sdImagestretched = 255*sdImageminsub./max(sdImageminsub(:));

% imagesc( sdImagestretched ); colormap gray; axis image; title('SD Image');


threshold = mean(sdImage(:))  + std(sdImage(:));
% figure(2); imagesc( sdImage > threshold ); colormap gray; axis image; title('Vessel Mask');

sdImageMasked = (sdImage > threshold);

end

