clear;
close all;


% Robert F Cooper
% Created 11-16-2017
% This script creates a super-video that is used to get an excellently
% accurate capillary map.


mov_path = pwd;


fNames = read_folder_contents(mov_path,'avi');

varImages = [];
%% Load the files
for f=1:length(fNames)

    if strcmp('_piped.avi', fNames{f}(end-length('_piped.avi')+1:end))
        fNames{f}
        temporal_stack_reader = VideoReader( fullfile(mov_path,fNames{f}) );

        temporal_stack=zeros(temporal_stack_reader.Height, temporal_stack_reader.Width,... 
                             round(temporal_stack_reader.Duration*temporal_stack_reader.FrameRate));
        i=1;
        while(hasFrame(temporal_stack_reader))
            temporal_stack(:,:,i) = double(readFrame(temporal_stack_reader));
            i=i+1;
        end
        
        [ ~, thisVarImage ] = tam_etal_capillary_func( temporal_stack );
        
        varImages = cat(3,varImages,thisVarImage);
    end
end

%%
undefsd = (isnan(varImages) | isinf(varImages) );

varImages( undefsd ) = 1;
gausfiltImages = zeros(size(varImages));
sum_map = zeros(size(varImages,1),size(varImages,2));

for i=1:size(varImages,3)
    
    sum_map = sum_map+(varImages(:,:,i)>0);
    gausfiltImages(:,:,i) = imgaussfilt(varImages(:,:,i),15,'FilterSize',71);
end

sdImage = exp( real(sqrt( sum(gausfiltImages,3)./(sum_map-1))) );

edgemask = 15;
sdImage(1:edgemask,:)=min(sdImage(:));
sdImage(:,1:edgemask)=min(sdImage(:));
sdImage(end-edgemask+1:end,:)=min(sdImage(:));
sdImage(:,end-edgemask+1:end)=min(sdImage(:));

notmin = sdImage~=min(sdImage(:));

sdImagemin = min(sdImage(notmin));
sdImageminsub = sdImage-sdImagemin;
sdImagestretched = 255*sdImageminsub./max(sdImageminsub(:));

% sdImagestretched= double(adapthisteq(uint8(sdImagestretched), 'NumTiles', [4 4], 'Distribution','uniform'));

figure(1); imagesc( sdImagestretched ); colormap gray; axis image;


threshold = mean(sdImagestretched(notmin))  + std(sdImagestretched(notmin));

imageMask = imclose(sdImagestretched>threshold, strel('disk',7));
imageMask(1:edgemask,:)=true;
imageMask(:,1:edgemask)=true;
imageMask(end-edgemask+1:end,:)=true;
imageMask(:,end-edgemask+1:end)=true;


figure(2); imagesc( imageMask ); colormap gray; axis image;