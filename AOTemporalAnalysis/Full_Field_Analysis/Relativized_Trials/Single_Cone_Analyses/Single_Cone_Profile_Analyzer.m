% 2018-08-29 Robert F Cooper
%
% This script analyzes the output from Rel_FF_Single_Cone_Analyses.
%


%% Load our data, and calculate some basic stats.
clear;
close all;

saveplots = false;
logmode = true;
DENSTOMETRY_THRESHOLD = 0.1;
RESPONSE_THRESHOLD = 0.3;

profile_dir = uigetdir(pwd);

single_cone_mat_files = read_folder_contents(profile_dir,'mat');

load(fullfile(profile_dir,single_cone_mat_files{1}),'allcoords');

stimAmp = nan(size(allcoords,1), length(single_cone_mat_files));
stimMedian = nan(size(allcoords,1), length(single_cone_mat_files));
stimTTP = nan(size(allcoords,1), length(single_cone_mat_files));
stimRespRange = nan(size(allcoords,1), length(single_cone_mat_files));
Prestim = nan(size(allcoords,1), length(single_cone_mat_files));

controlAmp = nan(size(allcoords,1), length(single_cone_mat_files));
controlMedian = nan(size(allcoords,1), length(single_cone_mat_files));

single_cone_response = nan(size(allcoords,1), length(single_cone_mat_files));
single_cone_control_response = nan(size(allcoords,1), length(single_cone_mat_files));
trial_validity = false(size(allcoords,1), length(single_cone_mat_files));

for i=1:length(single_cone_mat_files)
    single_cone_mat_files{i}
    load(fullfile(profile_dir,single_cone_mat_files{i}));
    
    stimAmp(:,i) = AmpResp;
    stimMedian(:,i) = MedianResp;
    stimTTP(:,i) = TTPResp;
%     stimRespRange(:,i) = stim_resp_range;
%     Prestim(:,i) = median(stim_prestim_means,2,'omitnan');
    
    controlAmp(:,i) = ControlAmpResp;
    controlMedian(:,i) = ControlMedianResp;    
    
    single_cone_response(:,i) = AmpResp;%+abs(MedianResp);
    single_cone_control_response(:,i) = ControlAmpResp;%;+abs(ControlMedianResp);
    trial_validity(:,i) = valid;
    
    if logmode 
        single_cone_response(:,i) = log10(single_cone_response(:,i)+1);
        single_cone_control_response(:,i) = log10(single_cone_control_response(:,i)+1);
    end
end

valid = all(~isnan(single_cone_response),2) & all(trial_validity,2);


return;

%% Individual Spatal maps

generate_spatial_map(single_cone_response, allcoords, valid, single_cone_mat_files,'_stimulus', saveplots);

generate_spatial_map(single_cone_control_response, allcoords, valid, single_cone_mat_files,'_control', saveplots);

%% Amplitude vs Median response
figure; hold on;
for i=1:length(single_cone_mat_files)
    plot(stimAmp(valid, i) ,abs(stimMedian(valid, i)),'.');
end
xlabel('Std dev reponse');
ylabel('Absolute Mean reponse');
if saveplots
    saveas(gcf, ['comparative_responses.png']); 
end

%% Bland-Altman plot

response_compare = [single_cone_response(:,1) single_cone_control_response(:,1)]
meanresp = mean(response_compare,2);
diffresp = diff(response_compare,[],2);

diffbias = mean(diffresp(:),'omitnan');
stddiff = 1.96*std(diffresp(:),'omitnan');
LOA = [diffbias-stddiff diffbias-stddiff;
       diffbias+stddiff diffbias+stddiff];
  
figure; plot(meanresp,diffresp,'.'); hold on;
currentsize = axis;
plot(currentsize(1:2), LOA(1,:),'r-.')
plot(currentsize(1:2), [diffbias diffbias],'b')
plot(currentsize(1:2), LOA(2,:),'r-.'); hold off;



%% Boxplot of the amplitudes from each intensity.
figure;
boxplot(single_cone_response,'notch','on');
ylabel('Stimulus amplitude');

if saveplots
    saveas(gcf, 'allresps_boxplot.png');
end

%% Histograms of the response from each mat file.

% [~, edges]=histcounts(single_cone_response(valid,1),10);
% 
% binwidth = diff(edges);

figure; hold on;
for i=1:length(single_cone_mat_files)
    if logmode
        histogram(single_cone_response(valid,i),'BinWidth',log10(1.2));
    else
        histogram(single_cone_response(valid,i),'BinWidth',0.2);
    end
end

xlabel('Aggregate Response');
ylabel('Number of Cones');

% axis([0 4 -1 10])
if saveplots
    saveas(gcf, 'allresps_histogramsplot.svg');
end

%% Vs plots

for i=1:length(single_cone_mat_files)
    figure; hold on;
    plot(single_cone_control_response(:,i),single_cone_response(:,i),'.');
    if logmode
        plot([-0.5 1.5],[-0.5 1.5],'k');
        axis equal;axis([-0.5 .5 -0.5 .5]); 
    else
        plot([-10 10],[-10 10],'k');
        axis equal; axis([-0.5 2 -0.5 15]); 
    end
    
    if saveplots
        saveas(gcf, [single_cone_mat_files{i}(1:end-4) '_VS_plot.png']);
    end
end

%% Display Cells under the 1:1 line

diffvalid = single_cone_response-single_cone_control_response;

lessthanvalid = diffvalid < abs(min(diffvalid))*0.8 & valid;

% lessthanvalid= (single_cone_response<single_cone_control_response) & valid;

% lessthanvalid = (single_cone_response<0.3) & valid & valid_densitometry;
% valid = valid & valid_densitometry;

for i=1:length(single_cone_mat_files)
    figure; hold on;    
    plot(single_cone_control_response(valid,i),single_cone_response(valid,i),'.');
    plot(single_cone_control_response(lessthanvalid(:,i)&valid,i),single_cone_response(lessthanvalid(:,i)&valid,i),'r.');
    plot([-10 10],[-10 10],'k');
         
    thelessthan{i} = find(lessthanvalid(:,i)==1);
    if logmode
        axis square;axis([-0.5 1.5 -0.5 1.5]); 
        xlabel('Log Control Response (mean control subtracted)')
        ylabel('Log Stimulus Response (mean control subtracted)');
    else
        axis equal;axis([-1 5 -1 20]);
        xlabel('Control Response (mean control subtracted)')
        ylabel('Stimulus Response (mean control subtracted)');
    end
    
    if saveplots
        saveas(gcf, [single_cone_mat_files{i}(1:end-4) '_VS_lessthan_plot.png']);
    end
    
    generate_spatial_map(single_cone_response(:,i), allcoords, lessthanvalid(:,i), single_cone_mat_files(i), '_VS', saveplots);
end

%% Angular coordinate histograms of VS plots.

for i=1:length(single_cone_mat_files)
    

    figure; hold on;
    histogram(single_cone_response(valid,i),40);
    title(strrep(single_cone_mat_files{i}(1:end-4),'_','\_'))
    if logmode
        xlabel('Log Stimulus response');
    else
        xlabel('Stimulus response');
    end
    ylabel('Number of cones');
    
    if saveplots
        saveas(gcf, [single_cone_mat_files{i}(1:end-4) '_resp_histogram.png']);
    end
end

%% Display results vs densitometry
% load('/local_data/Dropbox/General_Postdoc_Work/Dynamic_Densitometry/11049/Dynamic_Densitometry_combined_4_sec_545b25nm_3uW_20_single_cone_signals.mat')

lessthanvalid = (densitometry_fit_amplitude<=DENSTOMETRY_THRESHOLD) & valid & valid_densitometry;

% lessthanvalid = (single_cone_response<RESPONSE_THRESHOLD) & (densitometry_fit_amplitude<=DENSTOMETRY_THRESHOLD) & valid & valid_densitometry;
% lessthanvalid = all(lessthanvalid,2);

for i=1:length(single_cone_mat_files)
    figure; hold on;
    
    plot(single_cone_control_response(valid,i),single_cone_response(valid,i),'k.');
    plot(single_cone_control_response(lessthanvalid,i),single_cone_response(lessthanvalid,i),'r.');
    plot([-10 10],[-10 10],'k');
         
    thelessthan{i} = find(lessthanvalid==1);
    if logmode
        axis square;axis([-0.5 1.5 -0.5 1.5]); 
        xlabel('Log Control Response (mean control subtracted)')
        ylabel('Log Stimulus Response (mean control subtracted)');
    else
        axis equal;axis([-1 5 -1 20]);
        xlabel('Control Response (mean control subtracted)')
        ylabel('Stimulus Response (mean control subtracted)');
    end
    
    if saveplots
        saveas(gcf, [single_cone_mat_files{i}(1:end-4) '_VS_plot.png']);
    end
    
    generate_spatial_map(single_cone_response(:,i), allcoords, valid & valid_densitometry, single_cone_mat_files(i), '_VS', saveplots, ~lessthanvalid);
    
end

%% 

% lessthanvalid = (single_cone_response<RESPONSE_THRESHOLD) & (densitometry_fit_amplitude<=DENSTOMETRY_THRESHOLD) & valid & valid_densitometry;
 lessthanvalid = (densitometry_fit_amplitude<=DENSTOMETRY_THRESHOLD) & valid_densitometry;
numlowresp = 1;
figure; clf;
[V,C] = voronoin(allcoords,{'QJ'});
colors = ['rgbym'];
for j=1:length(single_cone_mat_files)
    
    for i=1:size(allcoords,1)
        
        vertices = V(C{i},:);

        if all(vertices(:,1)<max(allcoords(:,1))) && all(vertices(:,2)<max(allcoords(:,1))) ... % [xmin xmax ymin ymax] 
                                && all(vertices(:,1)>0) && all(vertices(:,2)>0) %&& ~isnan(single_cone_response(i,j))

            if all(lessthanvalid(i,:))
                patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', 'w');
                numlowresp = numlowresp + 1;
            elseif lessthanvalid(i,j)
                patch(V(C{i},1),V(C{i},2),ones(size(V(C{i},1))),'FaceColor', colors(j));
            end
        end
    end
end
axis image;
axis([0 max(allcoords(:,1)) 0 max(allcoords(:,2)) ])    
 set(gca,'Color','k');     

% if saveplots
%     saveas(gcf, [single_cone_mat_files{i}(1:end-4) '_agreement_plot.png']);
% end

%% Repeatability of timepoints 1 and 2.
lessthanvalid = (densitometry_fit_amplitude<=DENSTOMETRY_THRESHOLD) & valid & valid_densitometry;

lownotdens = find((lessthanvalid < (valid & valid_densitometry) ) & ...
             (single_cone_response(:,1)<0.45 & single_cone_response(:,2)<0.45));

figure; hold on;
plot(single_cone_response(:,1),single_cone_response(:,2),'k.');
plot([-10 10],[-10 10],'k');
axis equal; axis([-0.5 2 -0.5 15]);
axis([-0.5 2 -0.5 2]);
xlabel('Timepoint 1'); ylabel('Timepoint 2');
title('Responses between both time points.')
% plot(single_cone_response(lownotdens,1),single_cone_response(lownotdens,2),'b.');
plot(single_cone_response(lessthanvalid,1),single_cone_response(lessthanvalid,2),'r.');

if saveplots
    saveas(gcf, [single_cone_mat_files{1}(1:end-4) '_Repeat_plot_thresh_' num2str(DENSTOMETRY_THRESHOLD) '.png']);
end

figure; hold on;
% plot(mean(single_cone_response(lownotdens,:),2), densitometry_fit_amplitude(lownotdens),'b*')
plot(mean(single_cone_response(lessthanvalid,:),2), densitometry_fit_amplitude(lessthanvalid),'r*')
xlabel('Mean cone response (std dev + median)');
ylabel('Densitometry fit amplitude');

if saveplots
    saveas(gcf, [single_cone_mat_files{1}(1:end-4) '_dens_vs_mag_plot_thresh_' num2str(DENSTOMETRY_THRESHOLD) '.png']);
end

%% Prestimulus behavior vs normalized response range
figure; hold on;
for i=1:length(single_cone_mat_files)
    
    plot( Prestim(valid,i), stimRespRange(valid, i),'.');
        
end