
clear;
% close all;

load('450nW');
fitAmp_550nm = fitAmp; 
fitMean_550nm = fitMean;
fitAngle_550nm = fitAngle;

load('19000nW.mat');
fitAmp_675nm = fitAmp; 
fitMean_675nm = fitMean;
fitAngle_675nm = fitAngle;


% Just Amp
% allfits = [fitAmp_550nm fitAmp_675nm];

% Euclidian Distance
% allfits = [sqrt(fitAmp_550nm.^2   + fitMean_550nm.^2)...
%            sqrt(fitAmp_675nm.^2  + fitMean_675nm.^2)];

 
allfits = [ (fitAmp_550nm   + abs(fitMean_550nm)) ...
            (fitAmp_675nm  + abs(fitMean_675nm)) ];

valid = all(~isnan(allfits),2);

wavelengths = repmat( [550 675],[size(allfits,1) 1]);

slopes = allfits(:,1)-allfits(:,2);%diff(allfits,[],2);

upper_slope_thresh = quantile(slopes(valid),0.95); 
lower_slope_thresh = quantile(slopes(valid),0.05);



%% Individual Histograms
% To compare against the previously observed population-based dose-response,
% we will create histograms of the amplitudes from all cones for each 
% stimulus intensity, and calculate descriptive statistics from both of the
% histograms, such as their mean, median, and what percentage of their
% values are above 0. 

figure(1); histogram(allfits(valid, 1),'BinWidth',0.1,'Normalization','probability'); 
axis([-1 6 0 0.2]);
title(['\bf550nm: \rmMean: ' num2str(mean(allfits(valid, 1))) ...
       ' Median: ' num2str(median(allfits(valid, 1))) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '550nm_histo.png');

response_threshold = quantile(allfits(valid, 1),0.95)

figure(2); histogram(allfits(valid, 2),'BinWidth',0.1,'Normalization','probability');
axis([-1 6 0 0.12]); 
title(['\bf675nm: \rmMean: ' num2str(mean(allfits(valid, 2))) ...
       ' Median: ' num2str(median(allfits(valid, 2))) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '675nm_histo.png');


%% Individual Spatal maps

upper_thresh = quantile(allfits(:),0.95); 
lower_thresh = quantile(allfits(:),0.05);

for j=1:size(allfits,2)
    
    

    thismap = parula( ((upper_thresh-lower_thresh)*100)+2); 

    figure(3+j); clf;%imagesc(ref_image); hold on; colormap gray;
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

            thiscolorind = round((percentmax(i)-lower_thresh)*100)+1;

            vertices = V(C{i},:);

            if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0) %&& allfits(i,j)<lower_thresh %&& slopes(i)<lower_slope_thresh              
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




%% Slope plots

figure(3+j+1); histogram(slopes(valid),'BinWidth',0.1,'Normalization','probability');
% axis([-.5 3 0 0.15]); 
title(['\bf amplitude vs log-intensity slope: \rmMean: ' num2str(mean(slopes(valid)))...
       ' Median: ' num2str(median(slopes(valid))) ]);
xlabel('Slope');
ylabel('Probability');
saveas(gcf, 'slope_histo.png');

%% Slope spatial plot


thismap = parula(((upper_slope_thresh-lower_slope_thresh)*100)+2); 

figure(3+j+2); clf;%imagesc(ref_image); hold on; colormap gray;
axis image; hold on;

percentmax = zeros(size(allcoords,1));

[V,C] = voronoin(allcoords,{'QJ'});

for i=1:size(allcoords,1)
    
    if valid(i)
        percentmax(i) = slopes(i);
        
        if percentmax(i) > upper_slope_thresh
            percentmax(i) = upper_slope_thresh;
        elseif percentmax(i) < lower_slope_thresh
            percentmax(i) = lower_slope_thresh;
        end
        
        thiscolorind = round((percentmax(i)-lower_slope_thresh)*100)+1;
        
        vertices = V(C{i},:);
        
        if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0) && all(allfits(i,:)<lower_thresh)...
%                                 && slopes(i)<-1
%             plot(allcoords(i,1),allcoords(i,2),'.','Color', thismap(thiscolorind,:), 'MarkerSize', 15 );
            patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', thismap(thiscolorind,:));

        end
    end
end
colorbar
axis([min(allcoords(:,1)) max(allcoords(:,1)) min(allcoords(:,2)) max(allcoords(:,2)) ])
caxis([lower_slope_thresh upper_slope_thresh])
set(gca,'Color','k'); 
title('Slope spatial map')
hold off; drawnow;
saveas(gcf, 'spatial_map_slopes.png');

%% Vs plot

figure(8); clf; hold on;
plot(allfits(valid,1), allfits(valid,2),'k.');
xlabel('550nm Response');
ylabel('675nm Response');
title('550nm vs 675nm responses');hold off;
axis square; grid on;
saveas(gcf, '550_vs_675_response.png');

figure(9); clf; hold on;
compass( 100*(fitAmp_550nm(valid)-fitAmp_675nm(valid)), 100*(abs(fitMean_550nm(valid))-abs(fitMean_675nm(valid))) );
diffangle = atan2d(abs(fitMean_550nm(valid))-abs(fitMean_675nm(valid)), fitAmp_550nm(valid)-fitAmp_675nm(valid));
rose(diffangle*2*pi/360);
legend('Individual differences*100','Radial histogram')
title('550nm vs 675nm mean/stddev responses');hold off; axis square;
saveas(gcf, '550_vs_675_compass_rose_response.png');


%% Boxplot of the amplitudes from each intensity.

figure(10);
boxplot(allfits,'notch','on','Labels', {'675nm','550nm'});
xlabel('Stimulus wavelength');
ylabel('Stimulus amplitude');
title('Stimulus amplitudes for each stimulus wavelength')
saveas(gcf, 'allamps_boxplot.png');


