clear;
close all;


% Robert F Cooper
% Created 11-16-2017
% This script creates a super-video that is used to get an excellently
% accurate capillary map.

if ~exist('contains','builtin')
    contains = @(t,p)~isempty(strfind(t,p));
end
mov_path = pwd;


fNames = read_folder_contents(mov_path,'avi');

for f=1:length(fNames)
    if strcmp('_piped.avi', fNames{f}(end-length('_piped.avi')+1:end)) && ...
       contains(fNames{f}, 'confocal')
        temporal_stack_reader = VideoReader( fullfile(mov_path,fNames{f}) );
        varImages = nan(temporal_stack_reader.Height, temporal_stack_reader.Width, length(fNames));
        break;
    end
end


%% Load the files
parfor f=1:length(fNames)

    if strcmp('_piped.avi', fNames{f}(end-length('_piped.avi')+1:end)) && ...
       contains(fNames{f}, 'confocal')
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
        
        varImages(:,:,f) = thisVarImage;
    end
end
%%

tokeep = true(length(fNames),1);
for f=1:length(fNames)
    naner= isnan(varImages(:,:,f));
    tokeep(f) = ~all(naner(:));
end

varImages = varImages(:,:,tokeep);

%%
undefsd = (isnan(varImages) | isinf(varImages) );

edgemask = 15;
varImages( undefsd ) = 0;
varImages(1:edgemask,:,:) = 0;
varImages(:,1:edgemask,:) = 0;
varImages(end-edgemask-1:end,:,:) = 0;
varImages(:,end-edgemask-1:end,:) = 0;

gausfiltImages = zeros(size(varImages));
sum_map = zeros(size(varImages,1),size(varImages,2));

for i=1:size(varImages,3)
    
    nonzmask = imerode((varImages(:,:,i)>0),ones(7)); % Try and get rid of any edge artifacts from each map.
    
    sum_map = sum_map+nonzmask;
    gausfiltImages(:,:,i) = nonzmask.*imgaussfilt(varImages(:,:,i),15,'FilterSize',71);
    
%     figure(1); imagesc(gausfiltImages(:,:,i)); colormap gray;
%     pause;
end

sdImage = exp( real(sqrt( sum(gausfiltImages,3)./(sum_map-1))) );
%%
sdImage(isinf(sdImage))=1;
sdImage(isnan(sdImage))=1;


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

edges = 0:2:255;
binned = discretize(sdImagestretched(notmin),edges);

edges(mode(binned))

threshold = median(sdImagestretched(notmin))  + 0.8.*std(sdImagestretched(notmin));

capillary_mask = imclose(sdImagestretched>threshold, strel('disk',7));
capillary_mask(1:edgemask,:)=true;
capillary_mask(:,1:edgemask)=true;
capillary_mask(end-edgemask+1:end,:)=true;
capillary_mask(:,end-edgemask+1:end)=true;


figure(2); imagesc( capillary_mask ); colormap gray; axis image;
save('ALL_TRIALS_cap_map.mat','capillary_mask');
