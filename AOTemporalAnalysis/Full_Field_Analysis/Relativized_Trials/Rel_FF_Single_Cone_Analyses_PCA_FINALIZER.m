
clear;
% close all;

load('450nW.mat');


CRITICAL_REGION = 66:100;
% 450nW 
stddev_coeff_450nW = std_dev_coeff;
median_coeff_450nW = median_coeff;

stddev_stim_450nW = sqrt(stim_cell_var(:,CRITICAL_REGION));
stddev_control_450nW = sqrt(control_cell_var(:,CRITICAL_REGION));
median_stim_450nW = stim_cell_median(:,CRITICAL_REGION);
median_control_450nW = control_cell_median(:,CRITICAL_REGION);

valid_450nW = valid;

norm_nonan_ref = stddev_stim_450nW -stddev_control_450nW;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*stddev_coeff_450nW(:,1);
Stddev_450nW = sum(projected_ref,2);

norm_nonan_ref = median_stim_450nW;% -median_control_450nW;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*median_coeff_450nW(:,1);
Median_450nW = sum(projected_ref,2);

% 50nW
load('50nW.mat');
valid_50nW = valid;

norm_nonan_ref = sqrt(stim_cell_var(:,CRITICAL_REGION)) -sqrt(control_cell_var(:,CRITICAL_REGION));
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*stddev_coeff_450nW(:,1);
Stddev_50nW = sum(projected_ref,2);

norm_nonan_ref = stim_cell_median(:,CRITICAL_REGION);% -control_cell_median(:,CRITICAL_REGION);
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*median_coeff_450nW(:,1);
Median_50nW = sum(projected_ref,2);

% 0nW
load('0nW.mat');
valid_0nW = valid;

stddev_stim_0nW = sqrt(stim_cell_var(:,CRITICAL_REGION));
stddev_control_0nW = sqrt(control_cell_var(:,CRITICAL_REGION));
median_stim_0nW = stim_cell_median(:,CRITICAL_REGION);
median_control_0nW = control_cell_median(:,CRITICAL_REGION);

norm_nonan_ref = stddev_stim_0nW -stddev_control_0nW;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*stddev_coeff_450nW(:,1);
Stddev_0nW = sum(projected_ref,2);

norm_nonan_ref = median_stim_0nW;% -median_control_0nW;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*median_coeff_450nW(:,1);
Median_0nW = sum(projected_ref,2);

norm_nonan_ref = stddev_control_0nW;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*stddev_coeff_450nW(:,1);
tmp = sum(projected_ref,2);


% Create our response measure
allfits = [ (Stddev_0nW   + abs(Median_0nW)) ...
            (Stddev_50nW  + abs(Median_50nW)) ... 
            (Stddev_450nW + abs(Median_450nW)) ];
% allfits = [ (Stddev_0nW  ) ...
%             (Stddev_50nW  ) ... 
%             (Stddev_450nW ) ];        
       
valid = valid_0nW & valid_50nW & valid_450nW;

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



%% Determine each cone's slope

slopes = [ones(3,1) intensities']\allfits';

intercepts=slopes(1,:);
slopes=slopes(2,:);

upper_slope_thresh = quantile(slopes,0.95); 
lower_slope_thresh = quantile(slopes,0.02);

%% Reprojection - Stddev
figure(1); clf; 
delete('Stddev_low_responders.tif');

timeFromStim = (0:size(stddev_stim_450nW,2)-1)/16.66666;
large_n_valid = valid;% & ((allfits(:,3)+.1)<=allfits(:,1));
%(allfits(:,1)>0.8)); %  
% small_n_valid = valid & ((allfits(:,1)>.8) & (allfits(:,3)<4) ); 


for i=1:size(allcoords,1)
    if large_n_valid(i) 
        clf;
        warning off;
        for num_coeffs=0:3
            if num_coeffs > 0
                linreger = ([stddev_coeff_450nW(:,1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\stddev_stim_450nW');
                score_stim_450 = linreger(1:num_coeffs,:)';
                mu_stim_450 = linreger( num_coeffs+1,:)';

                linreger = ([stddev_coeff_450nW(:,1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\stddev_control_450nW');
                score_cont_450 = linreger(1:num_coeffs,:)';
                mu_cont_450 = linreger( num_coeffs+1,:)';
                
                linreger = ([stddev_coeff_450nW(:, 1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\(stddev_stim_450nW-stddev_control_450nW)');
                score_sub_450 = linreger(1:num_coeffs,:)';
                mu_sub_450 = linreger(num_coeffs+1, :)';

                linreger = ([stddev_coeff_450nW(:,1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\stddev_stim_0nW');
                score_stim_0 = linreger(1:num_coeffs,:)';
                mu_stim_0 = linreger( num_coeffs+1,:)';

                linreger = ([stddev_coeff_450nW(:,1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\stddev_control_0nW');
                score_cont_0 = linreger(1:num_coeffs,:)';
                mu_cont_0 = linreger( num_coeffs+1,:)';
                
                linreger = ([stddev_coeff_450nW(:, 1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\(stddev_stim_0nW-stddev_control_0nW)');
                score_sub_0 = linreger(1:num_coeffs,:)';
                mu_sub_0 = linreger(num_coeffs+1, :)';
                
            end
            

            % 450nW stimulus
            subplot(2,3,1); hold on;
            if num_coeffs > 0
                reprojected = score_stim_450(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_stim_450(i));
            else
                plot(timeFromStim, stddev_stim_450nW(i,:)); 
            end
            
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['450nW Stimulus (Response value: ' num2str(allfits(i,3)) ')']); axis([0 max(timeFromStim) 0 8]);        
            hold off;
            % 450nW control
            subplot(2,3,2); hold on;
            if num_coeffs > 0
                reprojected = score_cont_450(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_cont_450(i));
            else
                plot(timeFromStim, stddev_control_450nW(i,:)); 
            end
            
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['450nW Control']); axis([0 max(timeFromStim) 0 8]);
            hold off;
            
            % 450nW subbed
            subplot(2,3,3); hold on;
            if num_coeffs > 0
                reprojected = score_sub_450(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_sub_450(i));
            else
                subbed = (stddev_stim_450nW-stddev_control_450nW);
                plot(timeFromStim, subbed(i,:)); 
            end
            
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['450nW Subbed (Projected value: ' num2str(Stddev_450nW(i)) ')']); axis([0 max(timeFromStim) -1 7]);
            hold off;
            
            % 0nW "stimulus"
            subplot(2,3,4); hold on;
            if num_coeffs > 0
                reprojected = score_stim_0(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_stim_0(i));
            else
                plot(timeFromStim, stddev_stim_0nW(i,:)); 
            end
            
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['0nW "Stimulus" (Response value: ' num2str(allfits(i,1)) ')']); axis([0 max(timeFromStim) 0 6]);
            hold off;
            % 0nW "control"
            subplot(2,3,5); hold on;
            if num_coeffs > 0
                reprojected = score_cont_0(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_cont_0(i));
            else
                plot(timeFromStim, stddev_control_0nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['0nW "Control"']); axis([0 max(timeFromStim) 0 6]);        
            hold off;

            % 0nW subbed
            subplot(2,3,6); hold on;
            if num_coeffs > 0
                reprojected = score_sub_0(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_sub_0(i));
            else
                subbed = (stddev_stim_0nW-stddev_control_0nW);
                plot(timeFromStim, subbed(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['0nW Subbed (Projected value: ' num2str(Stddev_0nW(i)) ')']); axis([0 max(timeFromStim) -1 5]);        
            hold off;
            

%             figure(2); clf; imagesc(ref_image); hold on; text(allcoords(i,1),allcoords(i,2),['\leftarrow Cell#' num2str(i)])
            
        end
        
        drawnow;
%     pause;
%         imwrite(frame2im(getframe(gcf)), 'Stddev_low_responders.tif','WriteMode','append');
    end
    
    warning on;
end


%% Reprojection - Median
figure(1); clf;

% delete('Median_low_responders.tif');

for i=1:size(allcoords,1)
    if large_n_valid(i)  && i==78
        clf;
        warning off;
        for num_coeffs=0:3
            if num_coeffs > 0
                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\median_stim_450nW');
                score_stim_450 = linreger(1:num_coeffs,:)';
                mu_stim_450 = linreger(num_coeffs+1, :)';

                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\median_control_450nW');
                score_cont_450 = linreger(1:num_coeffs,:)';
                mu_cont_450 = linreger(num_coeffs+1, :)';
                
                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\(median_stim_450nW-median_control_450nW)');
                score_sub_450 = linreger(1:num_coeffs,:)';
                mu_sub_450 = linreger(num_coeffs+1, :)';

                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\median_stim_0nW');
                score_stim_0 = linreger(1:num_coeffs,:)';
                mu_stim_0 = linreger( num_coeffs+1, :)';

                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\median_control_0nW');
                score_cont_0 = linreger(1:num_coeffs,:)';
                mu_cont_0 = linreger(num_coeffs+1, :)';
                
                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\(median_stim_0nW-median_control_0nW)');
                score_sub_0 = linreger(1:num_coeffs,:)';
                mu_sub_0 = linreger(num_coeffs+1, :)';
            end

            % 450nW stimulus
            subplot(2,3,1); hold on;
            if num_coeffs > 0
                reprojected = score_stim_450(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_stim_450(i));
            else
                plot(timeFromStim, median_stim_450nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['450nW Stimulus (Response value: ' num2str(allfits(i,3)) ')']); axis([0 max(timeFromStim) -3 3]);        
            hold off;
            
            % 450nW control
            subplot(2,3,2); hold on;
            if num_coeffs > 0
                reprojected = score_cont_450(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_cont_450(i));
            else
                plot(timeFromStim, median_control_450nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['450nW Control']); axis([0 max(timeFromStim) -5 5]);
            hold off;
            
            % 450nW subbed
            subplot(2,3,3); hold on;
            if num_coeffs > 0
                reprojected = score_sub_450(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_sub_450(i));
            else
                subbed = (median_stim_450nW-median_control_450nW);
                plot(timeFromStim, subbed(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['450nW Subtracted (Projected value: ' num2str(abs(Median_450nW(i))) ')']); axis([0 max(timeFromStim) -5 5]);
            hold off;
            
            % 0nW "stimulus"
            subplot(2,3,4); hold on;
            if num_coeffs > 0
                reprojected = score_stim_0(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_stim_0(i));
            else
                plot(timeFromStim, median_stim_0nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['0nW "Stimulus" (Response value: ' num2str(allfits(i,1)) ')']); axis([0 max(timeFromStim) -3 3]);
            hold off;
            % 0nW "control"
            subplot(2,3,5); hold on;
            if num_coeffs > 0
                reprojected = score_cont_0(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_cont_0(i));
            else
                plot(timeFromStim, median_control_0nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['0nW "Control"']); axis([0 max(timeFromStim) -3 3]);        
            hold off;
            
            % 0nW subbed
            subplot(2,3,6); hold on;
            if num_coeffs > 0
                reprojected = score_sub_0(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_sub_0(i));
            else
                subbed = (median_stim_0nW-median_control_0nW);
                plot(timeFromStim, subbed(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['0nW Subtracted (Projected value: ' num2str(abs(Median_0nW(i))) ')']); axis([0 max(timeFromStim) -3 3]);
            hold off;
        end
%         pause;
        drawnow;
%         imwrite(frame2im(getframe(gcf)), 'Median_low_responders.tif','WriteMode','append');
    end
    warning on;
end

%% Individual Spatal maps

for j=1:size(allfits,2)
    
    upper_thresh = quantile(allfits(:),0.95); 
    lower_thresh = quantile(allfits(:),0.025);

    thismap = parula( ((upper_thresh-lower_thresh)*100)+2); 

    figure(3+j); clf;%imagesc(ref_image); hold on; colormap gray;
    axis image; hold on;

    percentmax = zeros(size(allcoords,1));

    [V,C] = voronoin(allcoords,{'QJ'});

    for i=1:size(allcoords,1)

        if large_n_valid(i)
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
    saveas(gcf, ['spatial_map_' num2str(j-1) '_lowresponse.png']);
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

figure(8); clf;imagesc(ref_image); hold on; colormap gray;
axis image; hold on;

percentmax = zeros(size(allcoords,1));

[V,C] = voronoin(allcoords,{'QJ'});
thesecoords=[];
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
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0) && slopes(i)<lower_slope_thresh
%             plot(allcoords(i,1),allcoords(i,2),'.','Color', thismap(thiscolorind,:), 'MarkerSize', 15 );
            patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', thismap(thiscolorind,:));

            thesecoords = [thesecoords; allcoords(i,:)];
        end
    end
end
colorbar
axis([min(allcoords(:,1)) max(allcoords(:,1)) min(allcoords(:,2)) max(allcoords(:,2)) ])
% caxis([lower_slope_thresh upper_slope_thresh])
set(gca,'Color','k'); 
title('Slope spatial map')
hold off; drawnow;
saveas(gcf, 'spatial_map_slopes.png');

%% Individual Increase maps - change to be based on profiles, 
% where any increase over 2sd kicks it out from the group


lower_fourfifty_thresh = quantile(allfits(:,3),0.05)
lower_fifty_thresh = quantile(allfits(:,2),0.05)

lowestfourfifty = (allfits(:,3)) < lower_fourfifty_thresh;
lowestfifty = (allfits(:,2)) < lower_fifty_thresh;

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
subplot(1,3,1);
plot(Stddev_0nW,abs(Median_0nW),'k.');
axis([-2 7 0 7]); title('0nW'); ylabel('Absolute median response'); xlabel('Std dev response');
subplot(1,3,2);
plot(Stddev_50nW,abs(Median_50nW),'k.');
axis([-2 7 0 7]); title('50nW'); ylabel('Absolute median response'); xlabel('Std dev response');
subplot(1,3,3);
plot(Stddev_450nW,abs(Median_450nW),'k.');
axis([-2 7 0 7]); title('450nW'); ylabel('Absolute median response'); xlabel('Std dev response');

xlabel('Std dev reponse');
ylabel('Absolute Median reponse');
saveas(gcf, ['comparative_responses.png']); 


%% Boxplot of the amplitudes from each intensity.

figure(14);
boxplot(allfits,'notch','on','Labels', {'0nW','50nW','450nW'}, 'Symbol', 'k.');
xlabel('Stimulus irradiance');
ylabel('Stimulus amplitude');
title('Stimulus amplitudes for each stimulus irradiance');
axis([0 4 -1 10])
saveas(gcf, 'allamps_boxplot.svg');

%% Vs plots
figure(15); clf; hold on;
plot(allfits(valid,1), allfits(valid,3),'k.');

plot([-20 160], [-20 160],'k');
xlabel('0nW Response');
ylabel('450nW Response');
title(['450nW vs 0nW responses: ' num2str(zero_to_fourfiftnW) '% increased.']);hold off;
axis equal; axis([-.5 1.5 -1 9]); grid on;
saveas(gcf, '450_vs_0nW_response.png');


figure(18); clf; hold on;
plot(allfits(valid,2), allfits(valid,3),'k.');
plot([-20 160], [-20 160],'k');
xlabel('50nW Response');
ylabel('450nW Response');

title('50nW vs 450nW responses');hold off;

axis equal; axis([-1 8 -1 8]);grid on;
saveas(gcf, '50_vs_450_response.png');

figure(19); clf; hold on;
plot(allfits(valid,1), allfits(valid,3),'.');
plot(allfits(valid,1), allfits(valid,2),'.');
plot(allfits(valid,2), allfits(valid,3),'.');

legend('50nW vs 0nW','450nW vs 0nW','450nW vs 50nW');
title('450nW vs 0nW and 450nW vs 50nW responses');hold off;

axis equal; axis([-1 10 -1 10]); grid on;
saveas(gcf, 'allvs_response.png');
saveas(gcf, 'allvs_response.svg');

%% Error plots

% figure(20); clf;
% histogram(allfitampserr(valid,:)); hold on; histogram(allfitmedianserr(valid,:));
% title('Bootstrapped Error- at least 10 trials'); legend('Amplitude Error','Median Error')
% saveas(gcf, 'Bootstrapped_Error.png');

% figure(19);clf;
% quiver(fitAmp_0nW(valid),abs(fitMean_0nW(valid)), fitAmp_50nW(valid)-fitAmp_0nW(valid), abs(fitMean_50nW(valid))-abs(fitMean_0nW(valid)),0 );
% hold on;
% quiver(fitAmp_50nW(valid),abs(fitMean_50nW(valid)), fitAmp_450nW(valid)-fitAmp_50nW(valid), abs(fitMean_450nW(valid))-abs(fitMean_50nW(valid)) );

function [ stddev, median, valid_cells, proj_stddev_profile, proj_median_profile ] = parse_and_project( mat_to_load, ref_stddev_coeff, ref_median_coeff )
% Robert F cooper 04-10-2018
%   This function loads in a given MATLAB mat file, and parses/processes
%   each of the mat-specific filenames into a format more easily handled by
%   the script.

load(mat_to_load);

% mu = mean(critical_period_std_dev_sub,1,'omitnan');
% norm_nonan_ref = bsxfun(@minus,critical_period_std_dev_sub,mu);
norm_nonan_ref = sqrt(stim_cell_var(:,critical_region))-sqrt(control_cell_var(:,critical_region));
projected_ref = norm_nonan_ref*ref_stddev_coeff;
stddev = projected_ref(:,1);

norm_nonan_ref = stim_cell_median(:,critical_region)-control_cell_median(:,critical_region);
projected_ref = norm_nonan_ref*ref_median_coeff;
median = projected_ref(:,1);

valid_cells = valid;

end
