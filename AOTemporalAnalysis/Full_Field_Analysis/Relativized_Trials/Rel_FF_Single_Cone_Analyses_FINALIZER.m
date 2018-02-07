
clear;
% close all;

load('0nW.mat');
fitAmp_0nW = fitAmp; 
fitMean_0nW = fitMean;
fitAngle_0nW = fitAngle;

load('50nW.mat');
fitAmp_50nW = fitAmp; 
fitMean_50nW = fitMean;
fitAngle_50nW = fitAngle;

load('450nW.mat');
fitAmp_450nW = fitAmp; 
fitMean_450nW = fitMean;
fitAngle_450nW = fitAngle;

% Just Amp
% allfits = [fitAmp_0nW fitAmp_50nW fitAmp_450nW];

allfits = [ (fitAmp_0nW   + abs(fitMean_0nW)) ...
            (fitAmp_50nW  + abs(fitMean_50nW)) ... 
            (fitAmp_450nW + abs(fitMean_450nW)) ];

valid = all(~isnan(allfits),2);

intensities = repmat( [0 log10(50) log10(450)],[size(allfits,1) 1]);

diffamp = diff(allfits,[],2);
bigdiff = allfits(:,3)-allfits(:,1);

total1=sum(~isnan(diffamp(:,1)));
total2=sum(~isnan(diffamp(:,2)));
total3=sum(~isnan(bigdiff));

zero_to_fiftynW = 100*sum(sign(diffamp( ~isnan(diffamp(:,1)) ,1)) == 1)./total1
fifty_to_fourfiftynW = 100*sum(sign(diffamp( ~isnan(diffamp(:,2)) ,2)) == 1)./total2
zero_to_fourfiftnW = 100*sum(sign(bigdiff( ~isnan(bigdiff))) == 1)./total3


%% Determine each cone's slope

slopes = [ones(3,1) intensities']\allfits';

intercepts=slopes(1,:);
slopes=slopes(2,:);

upper_slope_thresh = quantile(slopes,0.95); 
lower_slope_thresh = quantile(slopes,0.05);

%% Individual Histograms
% To compare against the previously observed population-based dose-response,
% we will create histograms of the amplitudes from all cones for each 
% stimulus intensity, and calculate descriptive statistics from both of the
% histograms, such as their mean, median, and what percentage of their
% values are above 0. 

figure(1); histogram(allfits(valid, 1),'BinWidth',0.1,'Normalization','probability'); 
axis([-1 6 0 0.2]);
title(['\bf0nW: \rmMean: ' num2str(mean(allfits(valid, 1))) ...
       ' Median: ' num2str(median(allfits(valid, 1))) ...
       ' %>0: ' num2str( 100*sum(allfits(valid, 1)>0)./ sum(valid) ) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '0nW_histo.png');

response_threshold = quantile(allfits(valid, 1),0.95)

figure(2); histogram(allfits(valid, 2),'BinWidth',0.1,'Normalization','probability');
axis([-1 6 0 0.12]); 
title(['\bf50nW: \rmMean: ' num2str(mean(allfits(valid, 2))) ...
       ' Median: ' num2str(median(allfits(valid, 2))) ...
       ' %>95th: ' num2str( 100*sum(allfits(valid, 2)>response_threshold)./ sum(valid) ) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '50nW_histo.png');


figure(3); histogram(allfits(valid, 3),'BinWidth',0.1,'Normalization','probability');
axis([-1 6 0 0.12]); 
title(['\bf450nW: \rmMean: ' num2str(mean(allfits(valid, 3))) ...
       ' Median: ' num2str(median(allfits(valid, 3))) ...
       ' %>95th: ' num2str( 100*sum(allfits(valid, 3)>response_threshold)./ sum(valid) ) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '450nW_histo.png');


disp([ num2str(100*sum(allfits(valid, 2)>response_threshold) / sum(valid)) '% of 50nW cone responses are over the 95th percentile of the 0nW condition'])
disp([ num2str(100*sum(allfits(valid, 3)>response_threshold) / sum(valid)) '% of 450nW cone responses are over the 95th percentile of the 0nW condition'])

%% Individual Spatal maps

for j=1:size(allfits,2)
    
    upper_thresh = quantile(allfits(:),0.95); 
    lower_thresh = quantile(allfits(:),0.05);

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
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0) 
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

figure(7); histogram(slopes(valid),'BinWidth',0.1,'Normalization','probability');
axis([-.5 3 0 0.15]); 
title(['\bf amplitude vs log-intensity slope: \rmMean: ' num2str(mean(slopes(valid)))...
       ' Median: ' num2str(median(slopes(valid))) ]);
xlabel('Slope');
ylabel('Probability');
saveas(gcf, 'slope_histo.png');

%% Slope spatial plot

thismap = parula(((upper_slope_thresh-lower_slope_thresh)*100)+2); 

figure(8); clf;%imagesc(ref_image); hold on; colormap gray;
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
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0)
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

%% Individual Low-response maps
    
for j=1:size(allfits,2)
    
    upper_thresh = quantile(allfits(:),0.95); 
    lower_thresh = quantile(allfits(:),0.05);

    thismap = parula( ((upper_thresh-lower_thresh)*100)+2); 

    figure(8+j); clf;%imagesc(ref_image); hold on; colormap gray;
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
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0) && allfits(i,j)<=lower_thresh
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
    saveas(gcf, ['low_response_map_' num2str(j-1) '.png']);
end

%% Individual Increase maps


lower_fourfifty_thresh = quantile(allfits(:,3)-allfits(:,1),0.05);
lower_fifty_thresh = quantile(allfits(:,2)-allfits(:,1),0.05);

lowestfourfifty = (allfits(:,3)-allfits(:,1)) < lower_fourfifty_thresh;
lowestfifty = (allfits(:,2)-allfits(:,1)) < lower_fifty_thresh;

figure(12); clf;%imagesc(ref_image); hold on; colormap gray;
axis image; hold on;

percentmax = zeros(size(allcoords,1));

[V,C] = voronoin(allcoords,{'QJ'});

for i=1:size(allcoords,1)

    if valid(i)

        vertices = V(C{i},:);

        if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0) 

            if lowestfifty(i) && lowestfourfifty(i)
                patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', 'r' );
            elseif lowestfifty(i)
                patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', 'b' );
            elseif lowestfourfifty(i)
                patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', 'g' );
            end

        end
    end
end

axis([min(allcoords(:,1)) max(allcoords(:,1)) min(allcoords(:,2)) max(allcoords(:,2)) ])
set(gca,'Color', 'k'); 
title(['Increasing map '])
hold off; drawnow;
saveas(gcf, ['increase_map.png']);

%% Plot each relationship on the plot
figure(13); hold on;
plot(fitAmp_50nW,abs(fitMean_50nW),'g.',...
     fitAmp_450nW,abs(fitMean_450nW),'r.',...
     fitAmp_0nW,abs(fitMean_0nW),'b.');


%% Boxplot of the amplitudes from each intensity.

figure(14);
boxplot(allfits,'notch','on','Labels', {'0nW','50nW','450nW'});
xlabel('Stimulus irradiance');
ylabel('Stimulus amplitude');
title('Stimulus amplitudes for each stimulus irradiance')
saveas(gcf, 'allamps_boxplot.png');

%% Vs plots
figure(15); clf; hold on;
plot(allfits(valid,1), allfits(valid,3),'k.');
xlabel('0nW Response');
ylabel('450nW Response');
title('0nW vs 450nW responses');hold off;
axis square; grid on;
saveas(gcf, '0_vs_450_response.png');

figure(16); clf; hold on;
compass( 100*(fitAmp_450nW(valid)-fitAmp_0nW(valid)), 100*(abs(fitMean_450nW(valid))-abs(fitMean_0nW(valid))) );
diffangle = atan2d(abs(fitMean_450nW(valid))-abs(fitMean_0nW(valid)), fitAmp_450nW(valid)-fitAmp_0nW(valid));
rose(diffangle*2*pi/360);
legend('Individual differences*100','Radial histogram')
title('0nW vs 450nW mean/stddev responses');hold off; axis square;
saveas(gcf, '0_vs_450_compass_rose_response.png');

