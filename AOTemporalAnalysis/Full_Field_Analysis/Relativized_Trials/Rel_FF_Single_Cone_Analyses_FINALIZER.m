
clear;
% close all;

load('allnW.mat');

allfits = [fitAmp_0nW fitAmp_50nW fitAmp_450nW];

valid = all(~isnan(allfits),2);

intensities = repmat( [0 log10(50) log10(450)],[size(allfits,1) 1]);


%% Individual Histograms
% To compare against the previously observed population-based dose-response,
% we will create histograms of the amplitudes from all cones for each 
% stimulus intensity, and calculate descriptive statistics from both of the
% histograms, such as their mean, median, and what percentage of their
% values are above 0. 

figure(1); histogram(fitAmp_0nW(valid),'BinWidth',0.1,'Normalization','probability'); 
axis([0 8 0 0.5]);
title(['\bf0nW: \rmMean: ' num2str(mean(fitAmp_0nW(valid)))...
       ' Median: ' num2str(median(fitAmp_0nW(valid))) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '0nW_histo.png');

figure(2); histogram(fitAmp_50nW(valid),'BinWidth',0.1,'Normalization','probability');
axis([0 8 0 0.15]); 
title(['\bf50nW: \rmMean: ' num2str(mean(fitAmp_50nW(valid)))...
       ' Median: ' num2str(median(fitAmp_50nW(valid))) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '50nW_histo.png');

figure(3); histogram(fitAmp_450nW(valid),'BinWidth',0.1,'Normalization','probability');
axis([0 8 0 0.15]); 
title(['\bf450nW: \rmMean: ' num2str(mean(fitAmp_450nW(valid)))...
       ' Median: ' num2str(median(fitAmp_450nW(valid))) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '450nW_histo.png');
%% Individual Spatal maps

for j=1:size(allfits,2)
    
    upper_thresh = 3; 
    lower_thresh = 0;

    thismap = parula((upper_thresh*100)+1); 

    figure(3+j); %imagesc(ref_image); hold on; colormap gray;
    axis image; hold on;

    percentmax = zeros(size(allcoords,1));

    [V,C] = voronoin(allcoords,{'QJ'});

    for i=1:size(allcoords,1)

        if valid(i)
            percentmax(i) = allfits(i,j);

            if percentmax(i) > upper_thresh
                percentmax(i) = upper_thresh;
            elseif percentmax(i) < lower_thresh
                percentmax(i) = lower_thresh;
            end

            thiscolorind = round(percentmax(i)*100)+1;

            vertices = V(C{i},:);

            if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0) %&& slopes(i)<0.1                
    %             plot(allcoords(i,1),allcoords(i,2),'.','Color', thismap(thiscolorind,:), 'MarkerSize', 15 );
                patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', thismap(thiscolorind,:));

            end
        end
    end
    colorbar
    axis([min(allcoords(:,1)) max(allcoords(:,1)) min(allcoords(:,2)) max(allcoords(:,2)) ])
    caxis([lower_thresh upper_thresh])
    set(gca,'Color','k'); 
    title(['Spatial map ' num2str(j-1)])
    hold off; drawnow;
    saveas(gcf, ['spatial_map_' num2str(j-1) '.png']);
end

%% Determine each cone's slope

slopes = intensities'\allfits';

slopes=slopes(1,:);


%% Slope plots

figure(3+j+1); histogram(slopes(valid),'BinWidth',0.1,'Normalization','probability');
axis([0 3 0 0.15]); 
title(['\bf amplitude vs log-intensity slope: \rmMean: ' num2str(mean(slopes(valid)))...
       ' Median: ' num2str(median(slopes(valid))) ]);
xlabel('Slope');
ylabel('Probability');
saveas(gcf, 'slope_histo.png');

%% Slope spatial plot

upper_thresh = 1.25; 
lower_thresh = 0;

thismap = parula((upper_thresh*100)+1); 

figure(3+j+2); %imagesc(ref_image); hold on; colormap gray;
axis image; hold on;

percentmax = zeros(size(allcoords,1));

[V,C] = voronoin(allcoords,{'QJ'});

for i=1:size(allcoords,1)
    
    if valid(i)
        percentmax(i) = slopes(i);
        
        if percentmax(i) > upper_thresh
            percentmax(i) = upper_thresh;
        elseif percentmax(i) < lower_thresh
            percentmax(i) = lower_thresh;
        end
        
        thiscolorind = round(percentmax(i)*100)+1;
        
        vertices = V(C{i},:);
        
        if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0) %&& slopes(i)<0.1                
%             plot(allcoords(i,1),allcoords(i,2),'.','Color', thismap(thiscolorind,:), 'MarkerSize', 15 );
            patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', thismap(thiscolorind,:));

        end
    end
end
colorbar
axis([min(allcoords(:,1)) max(allcoords(:,1)) min(allcoords(:,2)) max(allcoords(:,2)) ])
caxis([lower_thresh upper_thresh])
set(gca,'Color','k'); 
title('Slope spatial map')
hold off; drawnow;
saveas(gcf, 'spatial_map_slopes.png');

