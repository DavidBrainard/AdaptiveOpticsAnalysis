% AlignAOSLOImagesConv
%
% Some playing around with AO SLO frames provided by Jessica Morgan.
%
% 7/10/12  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Load in the frames
movieName = 'movie_4214';
aosloImagesObj = VideoReader([movieName '.avi']);
vidFrames = double(read(aosloImagesObj))/255;
nFrames = size(vidFrames,4);

%% Parameters
fullImageRows = size(vidFrames,1);
fullImageCols = size(vidFrames,2);
alignImageRows = round(3*fullImageRows/4);
alignImageCols = round(3*fullImageCols/4);
alignImageRowOffset = round(fullImageCols/8);
alignImageColOffset = round(fullImageCols/8);

%% Generate reference image
refImage = squeeze(vidFrames(:,:,1,round(nFrames/2)));
figure; clf;

%% Align the frames
refImageAlignPart = refImage(alignImageRowOffset: alignImageRowOffset+alignImageRows-1, ...
        alignImageColOffset:alignImageColOffset+alignImageCols-1);
[refRows,refCols] = size(refImageAlignPart);
alignedImages = zeros(alignImageRows,alignImageCols,nFrames);
unalignedImages = zeros(alignImageRows,alignImageCols,nFrames);

for j = 1:nFrames
    testImage = squeeze(vidFrames(:,:,1,j));
    
    % Figure out where test is relative to the reference, using cross-correlation
    % This is pretty right out of the documentation for normxcorr2
    testImageAlignPart = testImage(alignImageRowOffset:alignImageRowOffset+alignImageRows-1, ...
    	alignImageColOffset:alignImageColOffset+alignImageCols-1);
    unalignedImages(:,:,j) = testImageAlignPart;
    cc = normxcorr2(refImageAlignPart,testImage);
    [max_cc, imax] = max(abs(cc(:)));
    [ypeak, xpeak] = ind2sub(size(cc),imax(1));
    corr_offset(j,:) = [ (ypeak-refRows+1) (xpeak-refCols+1) ];
    alignedImages(:,:,j) = testImage(corr_offset(j,1)  + (0:alignImageRows-1) , ...
        corr_offset(j,2) + (0:alignImageCols- 1) );
    
    % Show what is happening
    subplot(1,3,1);
    imshow(squeeze(alignedImages(:,:,1)));
    subplot(1,3,2);
    imshow(testImageAlignPart);
    subplot(1,3,3);
    imshow(squeeze(alignedImages(:,:,j)));
    drawnow;
end

%% Average both the aligned and unligned versions
alignedAvgImage = mean(alignedImages,3);
alignedAvgImage = alignedAvgImage/max(alignedAvgImage(:));
unalignedAvgImage = mean(unalignedImages,3);
unalignedAvgImage = unalignedAvgImage/max(unalignedAvgImage(:));
figure; clf;
subplot(1,2,1);
imshow(unalignedAvgImage);
title('Average Unaligned');
subplot(1,2,2);
imshow(alignedAvgImage);
title('Average Aligned');
savefig('AveragedImagesConv.pdf',gcf,'pdf');

% Plot movement
figure; clf; hold on
plot(corr_offset(:,2)-alignImageColOffset,corr_offset(:,1)-alignImageRowOffset,'r');
plot(corr_offset(1,2)-alignImageColOffset,corr_offset(1,1)-alignImageRowOffset,'ro','MarkerFaceColor','r','MarkerSize',8);
axis([-30 30 -30 30]);
title('Eye Movements');
xlabel('Horizontal (pixels)');
ylabel('Vertical (pixels)');
savefig('EyeMovementsConv.pdf',gcf','pdf');

% Write a movie of the aligned frames
framesPerSecond = 15;
outFilename = [movieName '_alignedConv.avi'];
%writerObj = VideoWriter(outFilename,'Uncompressed AVI');
writerObj = VideoWriter(outFilename,'Motion JPEG AVI');
writerObj.FrameRate = framesPerSecond;
open(writerObj);
for j = 1:nFrames
    writeVideo(writerObj,squeeze([unalignedImages(:,:,j) alignedImages(:,:,j)]));
end
close(writerObj);


