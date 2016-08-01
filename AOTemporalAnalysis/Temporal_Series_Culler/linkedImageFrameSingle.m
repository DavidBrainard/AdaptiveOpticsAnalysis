function [ selected, cropregion ] = linkedImageFrameSingle( images, titles )
% FUNCTION [ selectionid ] = linkedImageFrameSingle( images, titles )
% Robert Cooper 05-10-2016
%   This frame allows a user to select a single image from a group of
%   images. This updated version
%     Copyright (C) 2016  Robert F Cooper
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%   Inputs:
%       @images: This is a cell array of 3d arrays which will be displayed
%       in the subplots.
%
%       @titles: This is a cell array containing matched titles (via index)
%                to the images
%
%   Outputs:
%       @selected: These are the images that will be transferred
%
%

multihand = figure('name',['Press Enter To Include, or Tab to exclude.']);
% update = @updateContents;
set(multihand,'KeyPressFcn',@keypress);
set(multihand,'HitTest','off');



%% Find the image with the largest bounds; scale the others off of it.
imsizes = cell2mat(cellfun(@size,images,'UniformOutput',0));

[maxval ind] = max(imsizes(:,2));

sliderwidth=0.02;

imfigs = imagesc(images{1}); colormap gray; axis image; axis off; title(titles{1},'Interpreter', 'none'); hold on;
imrect = rectangle('Position',[1 1 size(images{1},2)-1 size(images{1},1)-1],'LineWidth',3,'EdgeColor','g');
maskpoly = plot(-1,-1,'LineWidth',2,'Color','b'); hold off;
set(imfigs,'ButtonDownFcn',@butdwn);

numims = length(images);

if numims > 1
    sliderui = uicontrol(multihand, 'Style', 'slider','Min',1,'Max', numims,'Value', 1,'SliderStep',[ 1/(numims-1) 1/(numims-1) ],'Units','normalized','Position',[1-sliderwidth 0 sliderwidth 1 ], 'Callback',{@sliderMotion});
end
%% Create the image viewer
selectionid = 1; %length(images);
selected = ones(length(images),1);

mask = ones(size(images{1}));
summask = zeros(size(images{1}));
maskslices = cell(length(images),1);

% Make all of the masks.
for m=1:length(selected)           
    % Make a convex hull mask from each image
    out = regionprops(images{m}>1,'ConvexHull');
    out = {out.ConvexHull};
    out = cell2mat(out');

    if size(out,1) ~= 0 && size(out,2) ~= 0
        convregion = convhull(out(:,1), out(:,2));
        convregion = [out(convregion,1), out(convregion,2)];

        maskslices{m} = poly2mask(convregion(:,1),convregion(:,2),...
                            size(images{m},1), size(images{m},2));
        summask = summask + maskslices{m};
    end    
end

% figure(499); imagesc(summask);

suggestMask()

uiwait(multihand);


    % Run this whenever you want tto update the contents of the figure
    function [selection]=updateContents(newimages, newtitles)
%         % Empty out the figure
%         clf(multihand);
% 
%         images = newimages;
%         titles = newtitles;
%         
%         imsizes = cell2mat(cellfun(@size,images,'UniformOutput',0));
% 
%         [maxval ind] = max(imsizes(:,2));
% 
%         sliderwidth=0.02;
% %         ratios = (imsizes(:,2)./maxval);
% 
%         figure(multihand); % Make that handle current, for the subplots.
%         
%         imfigs = imagesc(images{1}); colormap gray; axis image; axis off; title(titles{1},'Interpreter', 'none')
%         imrect = rectangle('Position',[1 1 size(images{1},2)-1 size(images{1},1)-1],'LineWidth',4,'EdgeColor','w');
%         set(imfigs,'ButtonDownFcn',@butdwn);
%         
% 
%         numims = length(images);
% 
%         if length(plotloc) > 1
%             sliderui = uicontrol(multihand, 'Style', 'slider','Min',1,'Max', numims,'Value', numims,'SliderStep',[ 1/(numims-1) 1/(numims-1) ],'Units','normalized','Position',[1-sliderwidth 0 sliderwidth 1 ], 'Callback',{@sliderMotion});
%         end
% %         sliderui = uicontrol(multihand, 'Style', 'slider','Min',1,'Max', length(images),'Value', length(images),'Units','SliderStep',[1/length(images) 1/length(images)],'normalized','Position',[1-sliderwidth 0 sliderwidth 1 ], 'Callback',{@sliderMotion});
%         selectionid = length(images);
%         selected = -ones(length(images),1);
%         set(imrect,'EdgeColor','w');
%         
% 
%         uiwait(multihand);
%         selection = selected;
    end

    function keypress(src, evt)
        
%         get(gcf,'CurrentCharacter')
        
        if strcmp(get(gcf,'CurrentCharacter'),char(13)) % Character 13 is the enter character...
%             acceptlayer()

            % Find the largest incribed rectangle in our current mask.
            [C, h, w, mask] =FindLargestRectangles(mask,[1 1 0], [300 150]);
            
            % Find the coordinates for each corner of the rectangle, and
            % return them
            cropregion = regionprops(mask,'ConvexHull');
            cropregion = cropregion.ConvexHull;
            % TL TR BR BL
            cropregion = [floor(min(cropregion)); [ceil(max(cropregion(:,1))) floor( min(cropregion(:,2))) ]; ...
                          ceil(max(cropregion));  [floor(min(cropregion(:,1))) ceil( max(cropregion(:,2)))] ];
            cropregion(cropregion==0) = 1;
            cropregion( cropregion(3,1)>size(images{1},2) ) = size(images{1},2);
            cropregion( cropregion(3,2)>size(images{1},1) ) = size(images{1},1);
            selected = logical(selected);
            uiresume;
        elseif strcmp(get(gcf,'CurrentCharacter'),'r') % r for reject
            rejectlayer()
        elseif strcmp(evt.Key,'uparrow')
            if (selectionid+1) <= length(images)
                butdwn([], selectionid+1)
            end
        elseif strcmp(evt.Key,'downarrow')
            if (selectionid-1) > 0
                butdwn([], selectionid-1)
            end
        end
        
    end

    function suggestMask()
        mask = ones(size(images{1}));
        
        for i=1:length(maskslices)
            if selected(i) == 1 && ~isempty(maskslices{i})
                mask = mask.*maskslices{i};
                
            end
        end
        
%         figure(500); imagesc(mask); colormap gray;
        if any( mask(:) ~=0 )
            
            maskhull = regionprops(mask,'ConvexHull');

            set(maskpoly,'XData', maskhull.ConvexHull(:,1), 'YData', maskhull.ConvexHull(:,2));
            drawnow;
            % Use for final masking?
        end
    end

    function acceptlayer()
               
        selected(selectionid) = 1;
        
        
        if (selectionid+1) <= length(images)
            butdwn([], selectionid+1)
        end
        
%         if( all(selected > -1) )            
%             uiresume;
%         end
    end

    function rejectlayer()
                     
        selected(selectionid) = 0;
        
        suggestMask()
        if (selectionid+1) <= length(images)            
            butdwn([], selectionid+1)
        end
  
%         if( all(selected > -1) )            
%             uiresume;
%         end
    end

    function butdwn(src, evt)
       
        if ~isempty(src)
            selectionid = find(imfigs == src);
        else
            selectionid = evt;
        end

        % Update the figure positions
        %If it isn't selected in some form (rejected/accepted), then
        % paint it white
        if selected(selectionid) == -1
            set(imrect,'EdgeColor','w');
        elseif selected(selectionid) == 1
            set(imrect,'EdgeColor','g');
        elseif selected(selectionid) == 0
            set(imrect,'EdgeColor','r');
        end

        set(sliderui,'Value', selectionid )
        
        set(imfigs,'CData',images{selectionid});
        title( [titles{selectionid} ' Total Selected: ' num2str(sum(selected)) ]);
        set(imfigs,'ButtonDownFcn',@butdwn);
        
    end

    function sliderMotion(src, evt)
       
        neworigin = get(src,'Value');

        neworigin = round(neworigin);
        selectionid = neworigin;
        
        
        if selected(selectionid) == -1
            set(imrect,'EdgeColor','w');
        elseif selected(selectionid) == 1
            set(imrect,'EdgeColor','g');
        elseif selected(selectionid) == 0
            set(imrect,'EdgeColor','r');
        end
        
        set(imfigs,'CData',images{selectionid});
        set(imfigs,'ButtonDownFcn',@butdwn);
        title(titles{selectionid});
        set(sliderui,'Value', neworigin )
    end
end