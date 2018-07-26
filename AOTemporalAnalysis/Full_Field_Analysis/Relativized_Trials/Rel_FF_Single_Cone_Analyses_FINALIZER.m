
clear;
% close all;

load('0nW.mat');
fitAmp_0nW = AmpResp;
fitMedian_0nW = MedianResp;
fitTTP_0nW = TTPResp;
Prestim_0nW = mean(stim_prestim_means,2,'omitnan');

load('50nW.mat');
fitAmp_50nW = AmpResp;
fitMedian_50nW = MedianResp;
fitTTP_50nW = TTPResp;
Prestim_50nW = mean(stim_prestim_means,2,'omitnan');

load('450nW.mat');
fitAmp_450nW = AmpResp;
fitMedian_450nW = MedianResp;
fitTTP_450nW = TTPResp;
Prestim_450nW = mean(stim_prestim_means,2,'omitnan');

% Just Amp
allfits = [ (fitAmp_0nW   + abs(fitMedian_0nW)) ...
            (fitAmp_50nW  + abs(fitMedian_50nW)) ... 
            (fitAmp_450nW + abs(fitMedian_450nW)) ];
    

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

zero_to_50_inc = (allfits(valid,2)-allfits(valid,1)>0);
fifty_to_450_inc = (allfits(valid,3)-allfits(valid,2)>0);

percent_chance_to_always_increase = 100*sum(zero_to_50_inc & fifty_to_450_inc) / sum(valid)


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

% figure(1); histogram(allfits(valid, 1),'BinWidth',0.1,'Normalization','probability'); 
% axis([-1 6 0 0.2]);
% title(['\bf0nW: \rmMean: ' num2str(mean(allfits(valid, 1))) ...
%        ' Median: ' num2str(median(allfits(valid, 1))) ...
%        ' %>0: ' num2str( 100*sum(allfits(valid, 1)>0)./ sum(valid) ) ]);
% xlabel('Amplitude (std devs)');
% ylabel('Probability');
% saveas(gcf, '0nW_histo.png');
% 
% response_threshold = quantile(allfits(valid, 1),0.95)
% 
% figure(2); histogram(allfits(valid, 2),'BinWidth',0.1,'Normalization','probability');
% axis([-1 6 0 0.12]); 
% title(['\bf50nW: \rmMean: ' num2str(mean(allfits(valid, 2))) ...
%        ' Median: ' num2str(median(allfits(valid, 2))) ...
%        ' %>95th: ' num2str( 100*sum(allfits(valid, 2)>response_threshold)./ sum(valid) ) ]);
% xlabel('Amplitude (std devs)');
% ylabel('Probability');
% saveas(gcf, '50nW_histo.png');
% 
% 
% figure(3); histogram(allfits(valid, 3),'BinWidth',0.1,'Normalization','probability');
% axis([-1 6 0 0.12]); 
% title(['\bf450nW: \rmMean: ' num2str(mean(allfits(valid, 3))) ...
%        ' Median: ' num2str(median(allfits(valid, 3))) ...
%        ' %>95th: ' num2str( 100*sum(allfits(valid, 3)>response_threshold)./ sum(valid) ) ]);
% xlabel('Amplitude (std devs)');
% ylabel('Probability');
% saveas(gcf, '450nW_histo.png');


% disp([ num2str(100*sum(allfits(valid, 2)>response_threshold) / sum(valid)) '% of 50nW cone responses are over the 95th percentile of the 0nW condition'])
% disp([ num2str(100*sum(allfits(valid, 3)>response_threshold) / sum(valid)) '% of 450nW cone responses are over the 95th percentile of the 0nW condition'])

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
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0) && allfits(i,j)>=upper_thresh
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
plot(fitAmp_50nW,abs(fitMedian_50nW),'g.',...
     fitAmp_450nW,abs(fitMedian_450nW),'r.',...
     fitAmp_0nW,abs(fitMedian_0nW),'b.');
xlabel('Std dev reponse');
ylabel('Absolute Mean reponse');
saveas(gcf, ['comparative_responses.png']); 


%% Boxplot of the amplitudes from each intensity.

figure(14);
boxplot(allfits,'notch','on','Labels', {'0nW','50nW','450nW'});
xlabel('Stimulus irradiance');
ylabel('Stimulus amplitude');
title('Stimulus amplitudes for each stimulus irradiance')
saveas(gcf, 'allamps_boxplot.png');

%% Histograms of the amplitudes from each intensity.

[~, edges]=histcounts(allfits(valid,1),10);

binwidth = diff(edges);

figure(34);clf; hold on;
histogram(allfits(valid,1),'BinWidth',binwidth(1));
histogram(allfits(valid,2),'BinWidth',binwidth(1));
histogram(allfits(valid,3),'BinWidth',binwidth(1));
xlabel('Aggregate Response');
ylabel('Number of Cones');
axis square;
% axis([0 4 -1 10])
saveas(gcf, 'allamps_histogramsplot.svg');


%% Vs plots
figure(15); clf; hold on;
plot(allfits(valid,1), allfits(valid,3),'k.');
% errorbar(allfits(valid,1), allfits(valid,3), ...
%          allfitampserr(valid,3), allfitampserr(valid,3), ...
%          allfitmedianserr(valid,3), allfitmedianserr(valid,3),'ko');
plot([-20 160], [-20 160],'k');
xlabel('0nW Response');
ylabel('450nW Response');
title(['450nW vs 0nW responses: ' num2str(zero_to_fourfiftnW) '% increased.']);hold off;
% axis([-1 8 -1 8])
axis([-20 160 -20 160])
axis square; grid on;
saveas(gcf, '450_vs_0nW_response.png');

figure(16); clf; hold on;
compass( 100*(fitAmp_450nW(valid)-fitAmp_0nW(valid)), 100*(abs(fitMedian_450nW(valid))-abs(fitMedian_0nW(valid))) );
diffangle = atan2d(abs(fitMedian_450nW(valid))-abs(fitMedian_0nW(valid)), fitAmp_450nW(valid)-fitAmp_0nW(valid));
rose(diffangle*2*pi/360);
legend('Individual differences*100','Radial histogram')
title('0nW vs 450nW mean/stddev responses');hold off; 
% axis([-100 700 -100 700])
axis square;
saveas(gcf, '0_vs_450_compass_rose_response.png');

figure(17); clf; hold on;
compass( 100*(fitAmp_450nW(valid)-fitAmp_50nW(valid)), 100*(abs(fitMedian_450nW(valid))-abs(fitMedian_50nW(valid))) );
diffangle = atan2d(abs(fitMedian_450nW(valid))-abs(fitMedian_50nW(valid)), fitAmp_450nW(valid)-fitAmp_50nW(valid));
rose(diffangle*2*pi/360);
legend('Individual differences*100','Radial histogram')
title('50nW vs 450nW mean/stddev responses');hold off; 
% axis([-200 400 -200 400])
axis square;
saveas(gcf, '50_vs_450_compass_rose_response.png');

figure(18); clf; hold on;
plot(allfits(valid,2), allfits(valid,3),'k.');
xlabel('50nW Response');
ylabel('450nW Response');

title('50nW vs 450nW responses');hold off;
% axis([-1 8 -1 8])
axis square; grid on;
saveas(gcf, '50_vs_450_response.png');

figure(19); clf; hold on;
plot(allfits(valid,1), allfits(valid,3),'.');
plot(allfits(valid,2), allfits(valid,3),'.');

legend('450nW vs 0nW','450nW vs 50nW');
title('450nW vs 0nW and 450nW vs 50nW responses');hold off;
% axis([-1 8 -1 8])
axis square; grid on;
saveas(gcf, '0_vs_450_n_50_vs_450_response.png');
saveas(gcf, '0_vs_450_n_50_vs_450_response.svg');

%%
figure(20); clf;
histogram(allfitampserr(valid,:)); hold on; histogram(allfitmedianserr(valid,:));
title('Bootstrapped Error- at least 20 trials'); legend('Amplitude Error','Median Error')
saveas(gcf, 'Bootstrapped_Error.png');

% figure(19);clf;
% quiver(fitAmp_0nW(valid),abs(fitMean_0nW(valid)), fitAmp_50nW(valid)-fitAmp_0nW(valid), abs(fitMean_50nW(valid))-abs(fitMean_0nW(valid)),0 );
% hold on;
% quiver(fitAmp_50nW(valid),abs(fitMean_50nW(valid)), fitAmp_450nW(valid)-fitAmp_50nW(valid), abs(fitMean_450nW(valid))-abs(fitMean_50nW(valid)) );

