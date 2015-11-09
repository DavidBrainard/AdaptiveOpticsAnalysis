% AlignAOSLOImages
%
% Some playing around with AO SLO frames provided by Jessica Morgan.
% Try the Bergen method.
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

%% Make space
alignedImages = zeros(fullImageRows,fullImageCols,nFrames);
unalignedImages = zeros(fullImageRows,fullImageCols,nFrames);

for j = 1:nFrames
    testImage = squeeze(vidFrames(:,:,1,j));
    %M = estAffineMulti2(refImage,testImage);
    M = estAffineIter2(refImage,testImage,10);
    alignedImage = warpAffine2(testImage,M);
    
    unalignedImages(:,:,j) = testImage;
    alignedImages(:,:,j) = alignedImage;
    
    % Show what is happening
    subplot(1,3,1);
    imshow(squeeze(alignedImages(:,:,1)));
    subplot(1,3,2);
    imshow(squeeze(unalignedImages(:,:,j)));
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
savefig('AveragedImagesBergen.pdf',gcf,'pdf');

% Write a movie of the aligned frames
framesPerSecond = 15;
outFilename = [movieName '_alignedBergen.avi'];
%writerObj = VideoWriter(outFilename,'Uncompressed AVI');
writerObj = VideoWriter(outFilename,'Motion JPEG AVI');
writerObj.FrameRate = framesPerSecond;
open(writerObj);
for j = 1:nFrames
    writeVideo(writerObj,squeeze([unalignedImages(:,:,j) alignedImages(:,:,j)]));
end
close(writerObj);


