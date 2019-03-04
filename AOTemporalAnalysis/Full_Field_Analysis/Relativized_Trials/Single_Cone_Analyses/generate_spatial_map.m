function [] = generate_spatial_map(single_cone_response, coords, valid_cells, names, suffix, saveplots, tohighlight)
%generate_spatial_map(single_cone_response, coords, valid_cells, names)
%
% 2018-08-29 Robert F Cooper
%
% This function outputs a map of the cone responses.
%

if ~exist('tohighlight','var')
    tohighlight = zeros(size(valid_cells));
end

for j=1:size(single_cone_response,2)
    
    upper_thresh = quantile(single_cone_response(:,j),0.95); 
    lower_thresh = quantile(single_cone_response(:,j),0.05);

    thismap = parula( ((upper_thresh-lower_thresh)*100)+2); 

    figure; clf;%imagesc(ref_image); hold on; colormap gray;
    axis image; hold on;

    percentmax = zeros(size(coords,1));

    [V,C] = voronoin(coords,{'QJ'});

    for i=1:size(coords,1)

        if valid_cells(i)
            percentmax(i) = single_cone_response(i,j);

            if percentmax(i) > upper_thresh
                percentmax(i) = upper_thresh;
            elseif percentmax(i) < lower_thresh
                percentmax(i) = lower_thresh;
            end

            thiscolorind = round((percentmax(i)-lower_thresh)*100)+1;

            vertices = V(C{i},:);

            if ~isnan(thiscolorind) && all(vertices(:,1)<max(coords(:,1))) && all(vertices(:,2)<max(coords(:,1))) ... % [xmin xmax ymin ymax] 
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0)
                if tohighlight(i)
    %             plot(coords(i,1), coords(i,2),'.','Color', thismap(thiscolorind,:), 'MarkerSize', 15 );
                    patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', thismap(thiscolorind,:), 'EdgeColor','none');
                else
                    patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', thismap(thiscolorind,:));
                end
            end
        end
    end
    colorbar
    axis([min(coords(:,1)) max(coords(:,1)) min(coords(:,2)) max(coords(:,2)) ])
    caxis([lower_thresh upper_thresh])
    set(gca,'Color','k'); 
    title(strrep(names{j}(1:end-4),'_','\_'))
    hold off; drawnow;
    if saveplots
        set(gcf, 'Renderer', 'painters');
        saveas(gcf, [names{j}(1:end-4) suffix '_spatial_map.png']);
    end
end

end

