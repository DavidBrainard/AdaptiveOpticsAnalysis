% FindThatCone
%
% See if we can find where a particular cone is in an image
% by looking for a small piece of the image.
%
% 8/3/12  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Load in the frames
inMovieName = 'movie_4214';
outMovieBaseName = 'FindThatCone';
aosloImagesObj = VideoReader([inMovieName '.avi']);
vidFrames = double(read(aosloImagesObj))/255;
nFrames = size(vidFrames,4);

%% Parameters
transStr = 'translation';
coneLocRow = 233;
coneLocCol = 269;
linesBefore = 80;
lineHalfHeight = 15;
lineHalfWidth = 50;
lineRows = coneLocRow-lineHalfHeight:coneLocRow+lineHalfHeight;
lineCols = coneLocCol-lineHalfWidth:coneLocCol+lineHalfWidth;
showHalfHeight = 30;
showHalfWidth = 30;

%% Generate reference image and pull out the part we
% are going to use to convolve with
refImageIndex = 1;
refImage = squeeze(vidFrames(:,:,1,refImageIndex));
refPart = refImage(lineRows-linesBefore,lineCols);
showLineImage = refImage;
showLineImage(lineRows,lineCols) = 0;
showConeRows = coneLocRow-showHalfHeight:coneLocRow+showHalfHeight;
showConeCols = coneLocCol-showHalfWidth:coneLocCol+showHalfWidth;
showConeImage = refImage(showConeRows,showConeCols);
showConeImage(showHalfHeight+1,showHalfWidth+1) = 0;

%% Align the frames
[refRows,refCols] = size(refPart);
rawConeImages = zeros(2*showHalfHeight+1,2*showHalfWidth+1,nFrames);
foundConeImages = zeros(2*showHalfHeight+1,2*showHalfWidth+1,nFrames);

%% Try to find the cone in each frame
figure; clf;
for j = 1:nFrames
    testImage = squeeze(vidFrames(:,:,1,j));
    
    % Figure out where reference signature is in test.
    % This is pretty right out of the documentation for normxcorr2
    cc = normxcorr2(refPart,testImage);
    [max_cc, imax] = max(abs(cc(:)));
    [ypeak, xpeak] = ind2sub(size(cc),imax(1));
    corr_offset(j,:) = [ (ypeak-refRows+1) (xpeak-refCols+1) ];
    testConeLocRow = corr_offset(j,1) + lineHalfHeight + linesBefore;
    testConeLocCol = corr_offset(j,2) + lineHalfWidth;
    
    % Pull out raw ref location in test images
    rawConeImages(:,:,j) = testImage(showConeRows,showConeCols);
    rawConeImages(showHalfHeight+1,showHalfWidth+1,j) = 0;
    subplot(1,2,1); imshow(rawConeImages(:,:,j));
    title(sprintf('Raw image %d, ref is %d',j,refImageIndex));

    % Pull out the area round the found cone.
    testConeRows = testConeLocRow-showHalfHeight:testConeLocRow+showHalfHeight;
    testConeCols = testConeLocCol-showHalfWidth:testConeLocCol+showHalfWidth;
    foundConeImages(:,:,j) = testImage(testConeRows,testConeCols);
    foundConeImages(showHalfHeight+1,showHalfWidth+1,j) = 0;
    subplot(1,2,2); imshow(foundConeImages(:,:,j));
    title(sprintf('Cone centered image %d, ref is %d',j,refImageIndex));
    
    % Take a look as it goes by
    drawnow;
    %pause;
end

% Write a movie of the aligned frames
framesPerSecond = 15;
outFilename = sprintf([outMovieBaseName '_%s_%d_%d_%d.avi'],transStr,2*lineHalfHeight+1,2*lineHalfWidth+1,linesBefore);
%writerObj = VideoWriter(outFilename,'Uncompressed AVI');
writerObj = VideoWriter(outFilename,'Motion JPEG AVI');
writerObj.FrameRate = framesPerSecond;
open(writerObj);
for j = 1:nFrames
    writeVideo(writerObj,squeeze([rawConeImages(:,:,j) foundConeImages(:,:,j)]));
end
close(writerObj);


