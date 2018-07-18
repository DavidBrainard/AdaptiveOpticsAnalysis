
clear;
% close all;



NUM_COMPONENTS=3;
CRITICAL_REGION = 67:115;%66:100; %1:166; %
% 450nW
load('450nW.mat');
stddev_coeff_450nW = std_dev_coeff;
median_coeff_450nW = median_coeff;
std_explained_coeff_450nW = std_dev_explained;
std_explained_coeff_450nW = std_explained_coeff_450nW(1:NUM_COMPONENTS)';
med_explained_coeff_450nW = median_explained;
med_explained_coeff_450nW = med_explained_coeff_450nW(1:NUM_COMPONENTS)';

stddev_stim_450nW = sqrt(stim_cell_var(:,CRITICAL_REGION));
stddev_control_450nW = sqrt(control_cell_var(:,CRITICAL_REGION));
median_stim_450nW = stim_cell_median(:,CRITICAL_REGION);
median_control_450nW = control_cell_median(:,CRITICAL_REGION);

valid_450nW = valid;

% 50nW
load('50nW.mat');
valid_50nW = valid;

stddev_stim_50nW = sqrt(stim_cell_var(:,CRITICAL_REGION));
stddev_control_50nW = sqrt(control_cell_var(:,CRITICAL_REGION));
median_stim_50nW = stim_cell_median(:,CRITICAL_REGION);
median_control_50nW = control_cell_median(:,CRITICAL_REGION);

% 0nW
load('0nW.mat');
valid_0nW = valid;

stddev_stim_0nW = sqrt(stim_cell_var(:,CRITICAL_REGION));
stddev_control_0nW = sqrt(control_cell_var(:,CRITICAL_REGION));
median_stim_0nW = stim_cell_median(:,CRITICAL_REGION);
median_control_0nW = control_cell_median(:,CRITICAL_REGION);

valid = valid_0nW & valid_50nW & valid_450nW;

allcontrolstd = mean([stddev_control_450nW(valid,:);stddev_control_50nW(valid,:);stddev_control_0nW(valid,:);stddev_stim_0nW(valid,:)],'omitnan');
allcontrolmed = mean([median_control_450nW(valid,:);median_control_50nW(valid,:);median_control_0nW(valid,:);median_stim_0nW(valid,:)],'omitnan');

%plot(allcontrolstd); hold on; plot(allcontrolmed); hold off; axis([2 166 -2 4])
% title('All averaged control values');
% save('control_avgs.mat','allcontrolstd','allcontrolmed')

% Process things

% 450nW
norm_nonan_ref = stddev_stim_450nW -allcontrolstd;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*stddev_coeff_450nW(:,1:NUM_COMPONENTS);
Stddev_450nW = sum(projected_ref.*std_explained_coeff_450nW,2)./sum(std_explained_coeff_450nW);

norm_nonan_ref = median_stim_450nW -allcontrolmed;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*median_coeff_450nW(:,1:NUM_COMPONENTS);
Median_450nW = sum( (projected_ref.*med_explained_coeff_450nW) ,2)./sum(med_explained_coeff_450nW);

% 50nW
norm_nonan_ref = stddev_stim_50nW -allcontrolstd;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*stddev_coeff_450nW(:,1:NUM_COMPONENTS);
Stddev_50nW = sum(projected_ref.*std_explained_coeff_450nW,2)./sum(std_explained_coeff_450nW);

norm_nonan_ref = median_stim_50nW -allcontrolmed;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*median_coeff_450nW(:,1:NUM_COMPONENTS);
Median_50nW = sum(projected_ref.*med_explained_coeff_450nW,2)./sum(med_explained_coeff_450nW);

% 0nW
norm_nonan_ref = stddev_stim_0nW -allcontrolstd;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*stddev_coeff_450nW(:,1:NUM_COMPONENTS);
Stddev_0nW = sum(projected_ref.*std_explained_coeff_450nW,2)./sum(std_explained_coeff_450nW);

norm_nonan_ref = median_stim_0nW  -allcontrolmed;
projected_ref = (norm_nonan_ref-mean(norm_nonan_ref,2,'omitnan'))*median_coeff_450nW(:,1:NUM_COMPONENTS);
Median_0nW = sum(projected_ref.*med_explained_coeff_450nW,2)./sum(med_explained_coeff_450nW);


% Create our response measure
% allfits = [ sqrt(Stddev_0nW.^2   + abs(Median_0nW).^2 ) ...
%             sqrt(Stddev_50nW.^2  + abs(Median_50nW).^2 ) ... 
%             sqrt(Stddev_450nW.^2 + abs(Median_450nW).^2 ) ];      

allfits = [ (Stddev_0nW   + abs(Median_0nW)) ...
            (Stddev_50nW  + abs(Median_50nW) ) ... 
            (Stddev_450nW + abs(Median_450nW) ) ];  

% allfits = log10(allfits+1);

intensities = repmat( [0 log10(50) log10(450)],[size(allfits,1) 1]);

diffamp = diff(allfits(valid,:),[],2);
bigdiff = allfits(valid,3)-allfits(valid,1);

total1=sum(~isnan(diffamp(:,1)));
total2=sum(~isnan(diffamp(:,2)));
total3=sum(~isnan(bigdiff));

zero_to_fiftynW = 100*sum(sign(diffamp( ~isnan(diffamp(:,1)) ,1)) == 1)./total1
fifty_to_fourfiftynW = 100*sum(sign(diffamp( ~isnan(diffamp(:,2)) ,2)) == 1)./total2
zero_to_fourfiftnW = 100*sum(sign(bigdiff( ~isnan(bigdiff))) == 1)./total3

zero_to_50_inc = (allfits(valid,2)-allfits(valid,1)>0);
fifty_to_450_inc = (allfits(valid,3)-allfits(valid,2)>0);

lg_50=quantile(allfits(valid,2),0.95);
lg_450=quantile(allfits(valid,3),0.95);
large_n_valid = ((allfits(valid,2)>=lg_50) & (allfits(valid,3)>=lg_450) ); 

sm_50=quantile(allfits(valid,2),0.05);
sm_450=quantile(allfits(valid,3),0.05);

small_n_valid = ((allfits(valid,2)<=sm_50) | (allfits(valid,3)<=sm_450) ); 

%% Clumping analysis

% logResp = log(allfits+1);

rng('shuffle');

validcoords = allcoords(valid,:);

highestResp = (allfits(valid,3));
maxdist = max(max(triu(pdist2(validcoords, validcoords),1) ) );

edgevector = 0:(floor(maxdist))/100:floor(maxdist);
% edgevector = 0:(550)/100:550;

reflectance_edges = [min(highestResp) mean(highestResp,'omitnan')-std(highestResp,'omitnan')  mean(highestResp,'omitnan')+std(highestResp,'omitnan') max(highestResp)];
reflectance_classes = small_n_valid(valid); %discretize(highestResp,reflectance_edges); %small_n_valid

COI = 1;

% If discrete
withinDist = COI==reflectance_classes;

alldists = triu(pdist2(validcoords(withinDist,:), validcoords(withinDist,:)),1);

figure(3);
plot(allcoords(withinDist,1), allcoords(withinDist,2),'*')

alldists(alldists==0) =[];
histogram(alldists, edgevector);
base_values = histcounts(alldists,edgevector,'Normalization','cdf');


RUNS = 1000;
newvalues= zeros(RUNS, length(edgevector)-1);

parfor i=1:RUNS
   i
    randoResp = reflectance_classes( randperm(size(reflectance_classes,1)) );
    withinDist = COI==randoResp;

    alldists = triu(pdist2(validcoords(withinDist,:), validcoords(withinDist,:)),1);
    
    alldists(alldists==0) =[];
    [newvalues(i,:)] = histcounts(alldists,edgevector,'Normalization','cdf');

%  plot(validcoords(withinDist,1),validcoords(withinDist,2),'.')
%  pause;
end
runningMin = min(newvalues);
runningMax = max(newvalues);
runningMean = mean(newvalues);

figure(4);clf;
plot(edgevector, [0 runningMin],'k.-', edgevector, [0 runningMax],'k.-',edgevector,[0 runningMean],'k');
%     drawnow;
figure(1); clf;
plot(runningMean,base_values,'k', runningMin,runningMean,'k--', runningMax,runningMean,'k--');
axis equal; axis([0 1 0 1]);
saveas(gcf, 'upper_5_50n450_clumping.svg');

figure(2);clf;
plot(runningMean,base_values,'k', runningMin,runningMean,'k--', runningMax,runningMean,'k-.');
axis equal; axis([0 .05 0 .05]);
saveas(gcf, 'upper_5_50n450_clumping_zoomie.svg');



%% Determine each cone's slope

slopes = [ones(3,1) intensities']\allfits';

intercepts=slopes(1,:);
slopes=slopes(2,:);

upper_slope_thresh = quantile(slopes,0.95); 
lower_slope_thresh = quantile(slopes,0.05);

lowslopes= slopes<=0.05;

% overallslope=[allfits(valid,2) ones(size(allfits(valid,2),1),1)]\allfits(valid,3);
% x=-1:.1:10;
% plot(x, x.*overallslope(1)+overallslope(2),'k-.',allfits(valid,2), allfits(valid,3),'.')



%% Reprojection - Stddev
figure(1); clf; 
delete('Stddev_50greaterthan450_responders.tif');

timeFromStim = (0:size(stddev_stim_450nW,2)-1)/16.66666;




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

                linreger = ([stddev_coeff_450nW(:,1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\stddev_stim_50nW');
                score_stim_50 = linreger(1:num_coeffs,:)';
                mu_stim_50 = linreger( num_coeffs+1,:)';

                linreger = ([stddev_coeff_450nW(:,1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\stddev_control_50nW');
                score_cont_50 = linreger(1:num_coeffs,:)';
                mu_cont_50 = linreger( num_coeffs+1,:)';
                
                linreger = ([stddev_coeff_450nW(:, 1:num_coeffs) ones(size(stddev_coeff_450nW,1),1)]\(stddev_stim_50nW-stddev_control_50nW)');
                score_sub_50 = linreger(1:num_coeffs,:)';
                mu_sub_50 = linreger(num_coeffs+1, :)';
                
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
                reprojected = score_stim_50(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_stim_50(i));
            else
                plot(timeFromStim, stddev_stim_50nW(i,:)); 
            end
            
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['50nW Stimulus (Response value: ' num2str(allfits(i,2)) ')']); axis([0 max(timeFromStim) 0 6]);
            hold off;
            % 0nW "control"
            subplot(2,3,5); hold on;
            if num_coeffs > 0
                reprojected = score_cont_50(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_cont_50(i));
            else
                plot(timeFromStim, stddev_control_50nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['50nW Control']); axis([0 max(timeFromStim) 0 6]);        
            hold off;

            % 0nW subbed
            subplot(2,3,6); hold on;
            if num_coeffs > 0
                reprojected = score_sub_50(:,1:num_coeffs)*stddev_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_sub_50(i));
            else
                subbed = (stddev_stim_50nW-stddev_control_50nW);
                plot(timeFromStim, subbed(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Std dev Response');
            title(['50nW Subbed (Projected value: ' num2str(Stddev_50nW(i)) ')']); axis([0 max(timeFromStim) -1 5]);        
            hold off;
            

%             figure(2); clf; imagesc(ref_image); hold on; text(allcoords(i,1),allcoords(i,2),['\leftarrow Cell#' num2str(i)])
            
        end
        
        drawnow;
%     pause;
        imwrite(frame2im(getframe(gcf)), 'Stddev_50greaterthan450_responders.tif','WriteMode','append');
    end
    
    warning on;
end


%% Reprojection - Median
figure(1); clf;

delete('Median_50greaterthan450_responders.tif');

for i=1:size(allcoords,1)
    if large_n_valid(i)  
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

                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\median_stim_50nW');
                score_stim_50 = linreger(1:num_coeffs,:)';
                mu_stim_50 = linreger( num_coeffs+1, :)';

                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\median_control_50nW');
                score_cont_50 = linreger(1:num_coeffs,:)';
                mu_cont_50 = linreger(num_coeffs+1, :)';
                
                linreger = ([median_coeff_450nW(:, 1:num_coeffs) ones(size(median_coeff_450nW,1),1)]\(median_stim_50nW-median_control_50nW)');
                score_sub_50 = linreger(1:num_coeffs,:)';
                mu_sub_50 = linreger(num_coeffs+1, :)';
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
            title(['450nW Stimulus (Response value: ' num2str(allfits(i,3)) ')']); axis([0 max(timeFromStim) -5 5]);        
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
                reprojected = score_stim_50(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_stim_50(i));
            else
                plot(timeFromStim, median_stim_50nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['50nW Stimulus (Response value: ' num2str(allfits(i,2)) ')']); axis([0 max(timeFromStim) -3 3]);
            hold off;
            % 0nW "control"
            subplot(2,3,5); hold on;
            if num_coeffs > 0
                reprojected = score_cont_50(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_cont_50(i));
            else
                plot(timeFromStim, median_control_50nW(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['50nW Control']); axis([0 max(timeFromStim) -3 3]);        
            hold off;
            
            % 0nW subbed
            subplot(2,3,6); hold on;
            if num_coeffs > 0
                reprojected = score_sub_50(:,1:num_coeffs)*median_coeff_450nW(:,1:num_coeffs)';
                plot(timeFromStim,reprojected(i,:)+mu_sub_50(i));
            else
                subbed = (median_stim_50nW-median_control_50nW);
                plot(timeFromStim, subbed(i,:)); 
            end
            legend('Raw signal','1 Coefficient','2 Coefficients','3 Coefficients','Location','northwest');
            xlabel('Time from stimulus onset (s)'); ylabel('Median Response');
            title(['50nW Subtracted (Projected value: ' num2str(abs(Median_50nW(i))) ')']); axis([0 max(timeFromStim) -3 3]);
            hold off;
        end
%         pause;
        drawnow;
        imwrite(frame2im(getframe(gcf)), 'Median_50greaterthan450_responders.tif','WriteMode','append');
    end
    warning on;
end

%% Individual Spatal maps


for j=1:size(allfits,2)
    
    upper_thresh = quantile(allfits(:),0.95); 
    lower_thresh = quantile(allfits(:),0.05);

    thismap = parula( ((upper_thresh-lower_thresh)*100)+2); 

    figure(3+j); clf; %imagesc(ref_image); hold on; colormap gray;
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

%             if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
%                                     && all(vertices(:,1)>0) && all(vertices(:,2)>0) 
            if ~isnan(thiscolorind) && all(vertices(:,1)<395) && all(vertices(:,2)<450) ... % [xmin xmax ymin ymax] 
                                    && all(vertices(:,1)>110) && all(vertices(:,2)>165)
%             if ~isnan(thiscolorind) && all(vertices(:,1)<460) && all(vertices(:,2)<410) ... % [xmin xmax ymin ymax] 
%                                     && all(vertices(:,1)>100) && all(vertices(:,2)>50)                                
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
    set(gcf, 'Renderer', 'painters');
    saveas(gcf, ['spatial_map_' num2str(j-1) '_resp.svg']);
%      pause;
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
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0) %%&& slopes(i)<lower_slope_thresh

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
 set(gcf, 'Renderer', 'painters');
saveas(gcf, 'spatial_map_slopes.png');
% saveas(gcf, ['spatial_map_' num2str(j-1) '_lowslope.svg']);

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

        if  all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
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
figure(13); clf; hold on;
subplot(1,3,1);
plot(Stddev_0nW,abs(Median_0nW),'k.');
axis square; axis([-2 10 0 5]); grid on; title('0nW'); ylabel('Absolute median response'); xlabel('Std dev response');
subplot(1,3,2);
plot(Stddev_50nW,abs(Median_50nW),'k.');
axis square; axis([-2 10 0 5]); grid on; title('50nW'); ylabel('Absolute median response'); xlabel('Std dev response');
subplot(1,3,3);
plot(Stddev_450nW,abs(Median_450nW),'k.');
axis square; axis([-2 10 0 5]);  grid on; title('450nW'); ylabel('Absolute median response'); xlabel('Std dev response');


saveas(gcf, ['comparative_responses.png']); 
saveas(gcf, ['comparative_responses.svg']);

figure(23); clf; hold on;
errorbar(median(Stddev_0nW(valid)),median(abs(Median_0nW(valid))), median(abs(Median_0nW(valid)))-quantile(abs(Median_0nW(valid)),0.05),quantile(abs(Median_0nW(valid)),0.95)-median(abs(Median_0nW(valid))),...
                                                                   median(Stddev_0nW(valid))-quantile(Stddev_0nW(valid),0.05),quantile(Stddev_0nW(valid),0.95)-median(Stddev_0nW(valid)));
% axis square; axis([-2 10 0 5]); grid on; title('0nW'); ylabel('Absolute median response'); xlabel('Std dev response');
errorbar(median(Stddev_50nW(valid)),median(abs(Median_50nW(valid))), median(abs(Median_50nW(valid)))-quantile(abs(Median_50nW),0.05),quantile(abs(Median_50nW),0.95)-median(abs(Median_50nW(valid))),...
                                                                     median(Stddev_50nW(valid))-quantile(Stddev_50nW(valid),0.05),quantile(Stddev_50nW(valid),0.95)-median(Stddev_50nW(valid)));
% axis square; axis([-2 10 0 5]); grid on; title('50nW'); ylabel('Absolute median response'); xlabel('Std dev response');
errorbar(median(Stddev_450nW(valid)),median(abs(Median_450nW(valid))), median(abs(Median_450nW(valid)))-quantile(abs(Median_450nW(valid)),0.05),quantile(abs(Median_450nW(valid)),0.95)-median(abs(Median_450nW(valid))),...
                                                                       median(Stddev_450nW(valid))-quantile(Stddev_450nW(valid),0.05),quantile(Stddev_450nW(valid),0.95)-median(Stddev_450nW(valid)));
axis square; axis([-.5 2 0 1]);  grid on; title('Allerrors'); ylabel('Absolute median response'); xlabel('Std dev response');
saveas(gcf, ['comparative_mean_responses.svg']);

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
axis equal; %axis([-.5 1.5 -1 9]); grid on;
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

axis equal; %axis([-1 10 -1 10]); grid on;
saveas(gcf, 'allvs_response.png');
saveas(gcf, 'allvs_response.svg');

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

%% Error Spatal maps


load('450nW_bootstrapped.mat');

AvgResp_450 = Avg_Resp;
StdResp_450 = Std_Resp;
Avg_StddevResp_450 = Avg_StddevResp;
Std_StddevResp_450 = Std_StddevResp;
Avg_MedianResp_450 = Avg_MedianResp;
Std_MedianResp_450 = Std_MedianResp;

load('50nW_bootstrapped.mat');

AvgResp_50 = Avg_Resp;
StdResp_50 = Std_Resp;
Avg_StddevResp_50 = Avg_StddevResp;
Std_StddevResp_50 = Std_StddevResp;
Avg_MedianResp_50 = Avg_MedianResp;
Std_MedianResp_50 = Std_MedianResp;

load('0nW_bootstrapped.mat');

AvgResp_0 = Avg_Resp;
StdResp_0 = Std_Resp;
Avg_StddevResp_0 = Avg_StddevResp;
Std_StddevResp_0 = Std_StddevResp;
Avg_MedianResp_0 = Avg_MedianResp;
Std_MedianResp_0 = Std_MedianResp;


sm_50=quantile(AvgResp_50(valid),0.05);
sm_450=quantile(AvgResp_450(valid),0.05);

small_n_valid = ((AvgResp_50<=quantile(AvgResp_0(valid),0.95)) & (AvgResp_450<=quantile(AvgResp_0(valid),0.95)) );
small_n_valid = small_n_valid & ((StdResp_50<=max(StdResp_0(valid))) & (StdResp_450<=max(StdResp_0(valid))) );

lowrespond_percent=100*sum(small_n_valid)./size(allcoords,1)



figure(20); clf; %imagesc(ref_image); hold on; colormap gray;
axis image; hold on;

percentmax = zeros(size(allcoords,1));

[V,C] = voronoin(allcoords,{'QJ'});

upper_thresh = quantile(allfits(valid,3),0.95); 
lower_thresh = quantile(allfits(valid,1),0.05);
thismap = parula( ((upper_thresh-lower_thresh)*100)+2); 

for i=1:size(allcoords(valid),1)

    if small_n_valid(i)
        percentmax(i) = AvgResp_450(i);

        if percentmax(i) > upper_thresh
            percentmax(i) = upper_thresh;
        elseif percentmax(i) < lower_thresh
            percentmax(i) = lower_thresh;
        end

        thiscolorind = round((percentmax(i)-lower_thresh)*100)+1;

        vertices = V(C{i},:);

            if ~isnan(thiscolorind) && all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                    && all(vertices(:,1)>0) && all(vertices(:,2)>0) 
%         if ~isnan(thiscolorind) && all(vertices(:,1)<460) && all(vertices(:,2)<410) ... % [xmin xmax ymin ymax] 
%                                 && all(vertices(:,1)>100) && all(vertices(:,2)>50) 
                                
            patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', thismap(thiscolorind,:));

        end
    end
end
colorbar
axis([min(allcoords(:,1)) max(allcoords(:,1)) min(allcoords(:,2)) max(allcoords(:,2)) ])
caxis([lower_thresh upper_thresh])
set(gca,'Color','k'); 
title('Lower 5% of cones')
hold off; drawnow;
%     set(gcf, 'Renderer', 'painters');
%     saveas(gcf, ['spatial_map_' num2str(j-1) '_highresp.svg']);
%     pause;


figure;

plot(allfits(valid,1),AvgResp_0(valid),'.'); hold on;
plot(allfits(valid,2),AvgResp_50(valid),'.');
plot(allfits(valid,3),AvgResp_450(valid),'.');
plot([-1 12],[-1 12],'k')
xlabel('Unbootstrapped response')
ylabel('Bootstrapped response')
title('Bootstrapped vs unboostrapped responses');
axis equal; axis([-1 12 -1 12]);


