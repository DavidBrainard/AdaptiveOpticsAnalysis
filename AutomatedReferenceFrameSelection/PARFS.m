% Rob Cooper 06-30-2017
%
% This script is an moderately modified implementation of the algorithm outlined by Salmon et
% al 2017: "An Automated Reference Frame Selection (ARFS) Algorithm for 
% Cone Imaging with Adaptive Optics Scanning Light Ophthalmoscopy".
%
% At UPenn, we have decided to call it "PARFS: Pretty Accurate Reference
% Frame Selection"

clear;
close all;


NUM_REF_OUTPUT = 3;
MODALITY_WEIGHTS = [1/3 1/3 1/3]; % The weights applied to each modality. Adjust these values if you want the first modality (say, confocal) to carry more weight in the reference frame choice.
STRIP_SIZE = 40; % The size of the strip at which we'll analyze the distortion
BAD_STRIP_THRESHOLD = 0; % Having more bad strips than this will result in a frame's removal from consideration.
MIN_NUM_FRAMES_PER_GROUP=20; % A group must have more than this number of frames otherwise it will be dropped from consideration

locInd=1; % Location index

[stack_fname, mov_path] = uigetfile(fullfile(pwd,'*.avi'), 'MultiSelect', 'on');

% stack_fname = stack_fname;
% mov_path = mov_path;

for m=1 : length(stack_fname) 
    tic;
    
    refs{m} = extract_candidate_reference_frames(stack_fname{m}, mov_path, STRIP_SIZE, BAD_STRIP_THRESHOLD, MIN_NUM_FRAMES_PER_GROUP)
    
    toc;
end


%% Look for correspondence between all of the modalities.

intersected = refs{1};
for m=2 : length(stack_fname)    
    intersected = union(intersected, refs{m});    
end

for r=1:length(intersected)
    of_interest = intersected(r);
    
    whichind = length(stack_fname)*ones(length(stack_fname),1); % Weight heavily against a reference frame if it doesn't show in all modalities.
    for m=1 : length(stack_fname)
        rank = find(refs{m}==of_interest);
        if ~isempty(rank)
            whichind(m) = rank*MODALITY_WEIGHTS(m);
        end
    end
    average_rank(r) = sum(whichind);
end

[rankings, rankinds ] = sort(average_rank,2,'ascend');

intersected = intersected(rankinds);

intersected(1:3)