function [ update, selected, multihand ] = linkedImageFrameSingle( images, titles )
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
update = @updateContents;
set(multihand,'KeyPressFcn',@keypress);
set(multihand,'HitTest','off');



%% Find the image with the largest bounds; scale the others off of it.
imsizes = cell2mat(cellfun(@size,images,'UniformOutput',0));

[maxval ind] = max(imsizes(:,2));

sliderwidth=0.02;

imfigs = imagesc(images{1}); colormap gray; axis image; axis off; title(titles{1},'Interpreter', 'none'); hold on;
imrect = rectangle('Position',[1 1 size(images{1},2)-1 size(images{1},1)-1],'LineWidth',4,'EdgeColor','w');
maskpoly = plot(-1,-1,'LineWidth',2,'Color','b'); hold off;
set(imfigs,'ButtonDownFcn',@butdwn);

numims = length(images);

if numims > 1
    sliderui = uicontrol(multihand, 'Style', 'slider','Min',1,'Max', numims,'Value', numims,'SliderStep',[ 1/(numims-1) 1/(numims-1) ],'Units','normalized','Position',[1-sliderwidth 0 sliderwidth 1 ], 'Callback',{@sliderMotion});
end
%% Create the image viewer
selectionid = length(images);
selected = -ones(length(images),1);

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
            acceptlayer()
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
        for i=1:length(selected)           
            if selected(i) == 1
                
                % Make a convex hull mask from the accepted image
                out = regionprops(images{i}>1,'ConvexHull');
                tmpmask = poly2mask(out.ConvexHull(:,1),out.ConvexHull(:,2), size(images{i},1), size(images{i},2));
                
                
                mask = mask.*tmpmask;
            end            
        end
        maskhull = regionprops(mask,'ConvexHull');
        
        set(maskpoly,'XData', maskhull.ConvexHull(:,1));
        set(maskpoly,'YData', maskhull.ConvexHull(:,2));
        % Use for final masking?
%         [C, h, w, m] =FindLargestRectangles(tmpmask,[1 1 0], [300 150]);
                
%         figure(500); imagesc(mask); colormap gray;
    end

    function acceptlayer()
               
        selected(selectionid) = 1;
        suggestMask()
        
        if (selectionid+1) <= length(images)
            butdwn([], selectionid+1)
        end
        
        if( all(selected > -1) )            
            uiresume;
        end
    end

    function rejectlayer()
                     
        selected(selectionid) = 0;

        if (selectionid-1) > 0
            butdwn([], selectionid-1)
        end
  
        if( all(selected > -1) )            
            uiresume;
        end
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

        set(sliderui,'Value', neworigin )
    end
end