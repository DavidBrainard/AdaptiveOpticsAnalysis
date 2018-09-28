% 2018-08-29 Robert F Cooper
%
% This script analyzes the output from Rel_FF_Single_Cone_Analyses.
%


%% Load our data, and calculate some basic stats.
clear;
close all;

saveplots = false;
logmode = false;


profile_dir = uigetdir(pwd);

single_cone_mat_files = read_folder_contents(profile_dir,'mat');

load(fullfile(profile_dir,single_cone_mat_files{1}),'allcoords');

stimAmp = nan(size(allcoords,1), length(single_cone_mat_files));
stimMedian = nan(size(allcoords,1), length(single_cone_mat_files));
stimTTP = nan(size(allcoords,1), length(single_cone_mat_files));
Prestim = nan(size(allcoords,1), length(single_cone_mat_files));

controlAmp = nan(size(allcoords,1), length(single_cone_mat_files));
controlMedian = nan(size(allcoords,1), length(single_cone_mat_files));

single_cone_response = nan(size(allcoords,1), length(single_cone_mat_files));
single_cone_control_response = nan(size(allcoords,1), length(single_cone_mat_files));

for i=1:length(single_cone_mat_files)
    single_cone_mat_files{i}
    load(fullfile(profile_dir,single_cone_mat_files{i}));
    
    stimAmp(:,i) = AmpResp;
    stimMedian(:,i) = MedianResp;
    stimTTP(:,i) = TTPResp;
    Prestim(:,i) = mean(stim_prestim_means,2,'omitnan');
    
    controlAmp(:,i) = ControlAmpResp;
    controlMedian(:,i) = ControlMedianResp;    
    
    single_cone_response(:,i) = AmpResp+abs(MedianResp);
    single_cone_control_response(:,i) = ControlAmpResp+abs(ControlMedianResp);
    
    if logmode 
        single_cone_response(:,i) = log10(single_cone_response(:,i)+1);
        single_cone_control_response(:,i) = log10(single_cone_control_response(:,i)+1);
    end
end

valid = all(~isnan(single_cone_response),2);


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
    histogram(single_cone_response(valid,i),20);
end

xlabel('Aggregate Response');
ylabel('Number of Cones');
axis square;
% axis([0 4 -1 10])
if saveplots
    saveas(gcf, 'allresps_histogramsplot.svg');
end

%% Vs plots

for i=1:length(single_cone_mat_files)
    figure; hold on;
    plot(single_cone_control_response(:,i),single_cone_response(:,i),'.');
    plot([-10 10],[-10 10],'k');
    axis square; axis([-0.5 1.5 -0.5 1.5]); 
end

%% Display Cells under the 1:1 line

% lessthanvalid= (single_cone_response<single_cone_control_response) & valid;
lessthanvalid = (single_cone_response<1) & (densitometry_fit_amplitude<0.1) & valid & valid_densitometry;
valid = valid & valid_densitometry;

for i=1:length(single_cone_mat_files)
    figure; hold on;    
    plot(single_cone_control_response(valid,i),single_cone_response(valid,i),'.');
    plot(single_cone_control_response(lessthanvalid(:,i),i),single_cone_response(lessthanvalid(:,i),i),'r.');
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
        saveas(gcf, [single_cone_mat_files{i}(1:end-4) '_VS_plot.png']);
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

lessthanvalid = (densitometry_fit_amplitude<0.05) & valid & valid_densitometry;

for i=1:length(single_cone_mat_files)
    figure; hold on;    
    plot(densitometry_fit_amplitude(valid,i),single_cone_response(valid,i),'.');
    plot(densitometry_fit_amplitude(lessthanvalid(:,i),i),single_cone_response(lessthanvalid(:,i),i),'r.');
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
        saveas(gcf, [single_cone_mat_files{i}(1:end-4) '_VS_plot.png']);
    end
    
    generate_spatial_map(single_cone_response(:,i), allcoords, lessthanvalid(:,i), single_cone_mat_files(i), '_VS', saveplots);
end
