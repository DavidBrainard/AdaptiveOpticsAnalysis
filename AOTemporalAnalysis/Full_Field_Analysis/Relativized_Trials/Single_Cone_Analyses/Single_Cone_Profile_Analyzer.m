% 2018-08-29 Robert F Cooper
%
% This script analyzes the output from Rel_FF_Single_Cone_Analyses.
%

clear;
close all;

saveplots = false;


%% Load our data, and calculate some basic stats.
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
    
    single_cone_response(:,i) = log10(AmpResp+abs(MedianResp)+1);
    single_cone_control_response(:,i) = log10(ControlAmpResp+abs(ControlMedianResp)+1);
end

valid = all(~isnan(single_cone_response),2);


%% Individual Spatal maps

generate_spatial_map(single_cone_response, allcoords, valid, single_cone_mat_files, saveplots);

generate_spatial_map(single_cone_control_response, allcoords, valid, single_cone_mat_files, saveplots);

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
donelist=[];
for i=1:length(single_cone_mat_files)
    figure; hold on;
    plot(single_cone_control_response(:,i),single_cone_response(:,i),'.');
    plot([-10 10],[-10 10],'k');
    axis square; axis([-0.5 1.5 -0.5 1.5]); 
end

%% Display Cells under the 1:1 line
close all

lessthanvalid= (single_cone_response<single_cone_control_response) & valid;

donelist=[];
for i=1:length(single_cone_mat_files)
    figure; hold on;    
    plot(single_cone_control_response(valid,i),single_cone_response(valid,i),'.');
    plot(single_cone_control_response(lessthanvalid(:,i),i),single_cone_response(lessthanvalid(:,i),i),'r.');
    plot([-10 10],[-10 10],'k');
    axis square; axis([-0.5 1.5 -0.5 1.5]); 
    thelessthan{i} = find(lessthanvalid(:,i)==1);
    xlabel('Log Control Response (mean control subtracted)')
    ylabel('Log Stimulus Response (mean control subtracted)');
end

generate_spatial_map(single_cone_response, allcoords, lessthanvalid, single_cone_mat_files, saveplots);
