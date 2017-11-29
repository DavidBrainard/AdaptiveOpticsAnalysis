
clear;
% close all;

load('0nW.mat');
load('50nW.mat');
load('450nW.mat');

allfits = [fitAmp_0nW fitAmp_50nW fitAmp_450nW];

valid = all(~isnan(allfits),2);

intensities = repmat( [0 log10(50) log10(450)],[size(allfits,1) 1]);

diffamp = diff(allfits,[],2);
bigdiff = fitAmp_450nW-fitAmp_0nW;

total1=sum(~isnan(diffamp(:,1)));
total2=sum(~isnan(diffamp(:,2)));
total3=sum(~isnan(bigdiff));

zero_to_fiftynW = 100*sum(sign(diffamp( ~isnan(diffamp(:,1)) ,1)) == 1)./total1
fifty_to_fourfiftynW = 100*sum(sign(diffamp( ~isnan(diffamp(:,2)) ,2)) == 1)./total2
zero_to_fourfiftnW = 100*sum(sign(bigdiff( ~isnan(bigdiff))) == 1)./total3


%% Boxplot of the amplitudes from each intensity.

figure(10);
boxplot(allfits,'notch','on','Labels', {'0nW','50nW','450nW'});
xlabel('Stimulus irradiance');
ylabel('Stimulus amplitude');
title('Stimulus amplitudes for each stimulus irradiance')
saveas(gcf, 'allamps_boxplot.png');

%% Determine each cone's slope

slopes = intensities'\allfits';

slopes=slopes(1,:);

%% Individual Histograms
% To compare against the previously observed population-based dose-response,
% we will create histograms of the amplitudes from all cones for each 
% stimulus intensity, and calculate descriptive statistics from both of the
% histograms, such as their mean, median, and what percentage of their
% values are above 0. 

figure(1); histogram(fitAmp_0nW(valid),'BinWidth',0.1,'Normalization','probability'); 
axis([-1 6 0 0.2]);
title(['\bf0nW: \rmMean: ' num2str(mean(fitAmp_0nW(valid)))...
       ' Median: ' num2str(median(fitAmp_0nW(valid))) ...
       ' %>0: ' num2str( 100*sum(fitAmp_0nW(valid)>0)./ sum(valid) ) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '0nW_histo.png');

response_threshold = quantile(fitAmp_0nW(valid),0.95)

figure(2); histogram(fitAmp_50nW(valid),'BinWidth',0.1,'Normalization','probability');
axis([-1 6 0 0.12]); 
title(['\bf50nW: \rmMean: ' num2str(mean(fitAmp_50nW(valid)))...
       ' Median: ' num2str(median(fitAmp_50nW(valid))) ...
       ' %>95th: ' num2str( 100*sum(fitAmp_50nW(valid)>response_threshold)./ sum(valid) ) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '50nW_histo.png');


figure(3); histogram(fitAmp_450nW(valid),'BinWidth',0.1,'Normalization','probability');
axis([-1 6 0 0.12]); 
title(['\bf450nW: \rmMean: ' num2str(mean(fitAmp_450nW(valid)))...
       ' Median: ' num2str(median(fitAmp_450nW(valid))) ...
       ' %>95th: ' num2str( 100*sum(fitAmp_450nW(valid)>response_threshold)./ sum(valid) ) ]);
xlabel('Amplitude (std devs)');
ylabel('Probability');
saveas(gcf, '450nW_histo.png');


disp([ num2str(100*sum(fitAmp_50nW(valid)>response_threshold) / sum(valid)) '% of 50nW cone responses are over the 95th percentile of the 0nW condition'])
disp([ num2str(100*sum(fitAmp_450nW(valid)>response_threshold) / sum(valid)) '% of 450nW cone responses are over the 95th percentile of the 0nW condition'])

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
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0) %&& slopes(i)<0.05             
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
%     saveas(gcf, ['spatial_map_' num2str(j-1) '.png']);
end




%% Slope plots

figure(3+j+1); histogram(slopes(valid),'BinWidth',0.1,'Normalization','probability');
axis([-.5 3 0 0.15]); 
title(['\bf amplitude vs log-intensity slope: \rmMean: ' num2str(mean(slopes(valid)))...
       ' Median: ' num2str(median(slopes(valid))) ]);
xlabel('Slope');
ylabel('Probability');
saveas(gcf, 'slope_histo.png');

%% Slope spatial plot

upper_thresh = quantile(slopes,0.95); 
lower_thresh = quantile(slopes,0.05);

thismap = parula(((upper_thresh-lower_thresh)*100)+2); 

figure(3+j+2); clf;%imagesc(ref_image); hold on; colormap gray;
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
        
        thiscolorind = round((percentmax(i)-lower_thresh)*100)+1;
        
        vertices = V(C{i},:);
        
        if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0) %&& slopes(i)<0.05                
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

