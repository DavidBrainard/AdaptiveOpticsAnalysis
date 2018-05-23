clear;
close all;


%% Filename determination and handling
[confocal_fname, pathname] = uigetfile('*.tif', 'Select the trial average images you wish to relativize:', 'MultiSelect','on');

if ~iscell(confocal_fname)
    confocal_fname={confocal_fname};
end

for i=1:length(confocal_fname)
    
    trial_im{i} = imread(fullfile(pathname,confocal_fname{i}));
    imsize(i,:) = size(trial_im{i});
end

%%
cropsize = min(imsize);

cropped_ims = zeros(cropsize(1),cropsize(2),length(confocal_fname));

for i=1:length(confocal_fname)
    
    residualleft = ceil((imsize(i,:)-cropsize)/2);
    residualright = floor((imsize(i,:)-cropsize)/2);
    
    world_ref{i} = imref2d(size(cropped_ims), [1+residualleft(2), (imsize(i,2)-residualright(2))], ...
                                              [1+residualleft(1), (imsize(i,1)-residualright(1))] );
    cropped_ims(:,:,i) = trial_im{i}( 1+residualleft(1):(imsize(i,1)-residualright(1)),...
                                      1+residualleft(2):(imsize(i,2)-residualright(2)) );
end

[optimizer, metric]  = imregconfig('multimodal');
[monooptimizer, monometric]  = imregconfig('monomodal');

tforms = cell(length(trial_im),length(trial_im));
% reg_ims = cell(length(trial_im),length(trial_im));
% Find the transform from each image to every other
% for i=1:length(trial_im)
ref_im = 1;
i=ref_im;
    for j=1:length(trial_im)
        if i~=j
            tic;
            tforms{i,j} = imregtform(trial_im{j},trial_im{ref_im},... % First get close
                                     'rigid',optimizer,metric, 'PyramidLevels',3);
            tforms{i,j} = imregtform(trial_im{j},trial_im{ref_im},... % Then tweak for affine
                                     'affine',monooptimizer,monometric, 'PyramidLevels',1,'InitialTransformation',tforms{i,j});
            toc;

%             tforms{i,j}.T
%             reg_ims(:,:,j) = imwarp(trial_im{j}, imref2d(size(trial_im{j})), tforms{1,j},'OutputView', imref2d(size(trial_im{1})) );
%             figure(1);imshowpair(trial_im{1},reg_ims(:,:,j));

        end
    end
% end

%%
stk_name = 'Relativized_Trial_Stack.tif';
reg_ims = repmat(trial_im{ref_im},[1 1 length(trial_im)]);
imwrite(reg_ims(:,:,ref_im),stk_name)

for j=1:length(trial_im)
    if j ~= ref_im
        reg_ims(:,:,j) = imwarp(trial_im{j}, imref2d(size(trial_im{j})), tforms{ref_im,j}, 'OutputView', imref2d(size(trial_im{ref_im})) );
        imwrite(reg_ims(:,:,j),stk_name,'WriteMode','append')
        figure(1); imshowpair(trial_im{ref_im}, reg_ims(:,:,j)); pause(.5)
        
    end
end


% Make a sum map, and create an average image.
sum_map = zeros(size(reg_ims(:,:,1)));

for f=1:size(reg_ims,3)
    frm_nonzeros = (reg_ims(:,:,f)>0); 
    sum_map = sum_map+frm_nonzeros;
end

% matchexp = '_\d*_ref_\d*_'; %String match that is expected to show up in the filename of each image. E.g. '_0018_ref_7_'
% avgname = regexprep(confocal_fname{1},matchexp,'_ALL_TRIALS_');
% matchexp = '_n\d*_';
% avgname = regexprep(avgname,matchexp,'_');

imwrite(uint8(sum(double(reg_ims),3)./sum_map), fullfile(pathname,[confocal_fname{i}(1:end-8) '_ALL_TRIALS.tif']));


%% Transform all of the videos associated with each avg image
for i=1:length(confocal_fname)
    mov_name_in = fullfile(pathname,[confocal_fname{i}(1:end-8) '.avi']);
    mov_name_out = fullfile(pathname,[confocal_fname{i}(1:end-8) '_piped.avi']);
    
    confocal_vidin = VideoReader( mov_name_in );
    confocal_vidout = VideoWriter( mov_name_out, 'Grayscale AVI' );
    
    open(confocal_vidout);    
    while hasFrame(confocal_vidin)
        frm_in = readFrame(confocal_vidin);
        
        if ~isempty(tforms{ref_im,i})
            writeVideo(confocal_vidout, imwarp(frm_in, imref2d(size(frm_in)), tforms{ref_im,i},...
                                                       'OutputView', imref2d(size(trial_im{ref_im})) ));
        else
            writeVideo(confocal_vidout, frm_in);
        end
    end
    close(confocal_vidout);
    
end


