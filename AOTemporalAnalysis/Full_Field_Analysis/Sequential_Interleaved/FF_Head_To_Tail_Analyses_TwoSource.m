function [fitCharacteristics]=FF_Aggregate_Analyses_TwoSource(stimRootDir, controlRootDir)
% Robert F Cooper
% 06-20-2017
% This script calculates pooled variance across a set of given signals.

if ~exist('stimRootDir','var')
    close all force;
    stimRootDir = uigetdir(pwd, 'Select the directory containing the stimulus profiles');
    controlRootDir = uigetdir(pwd, 'Select the directory containing the control profiles');
end

profileSDataNames = read_folder_contents(stimRootDir,'mat');
profileCDataNames = read_folder_contents(controlRootDir,'mat');

thatstimmax=0;
thatcontrolmax=0;
%% Code for determining variance across all signals at given timepoint

allmax=500;
num_control_cones = 0;
num_stim_cones =0;
max_cones = 0;
min_cones = 10000000000000;
mean_control_reflectance = zeros(500,1);
num_contributors_control_reflectance = zeros(500,1);

for j=1:length(profileCDataNames)

    profileCDataNames{j}
    load(fullfile(controlRootDir,profileCDataNames{j}));
    
    % Remove empty cells
    norm_control_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance)  );
    control_cell_times            = cell_times( ~cellfun(@isempty,cell_times) );

    
    num_control_cones = num_control_cones+length(control_cell_times);
    
    thatcontrolmax(j) = max( cellfun(@max,control_cell_times) );

    % Pooled variance of all cells before first stimulus
    [ ref_variance_control{j},ref_control_times{j}, ref_control_count{j} ] = reflectance_pooled_variance( control_cell_times, norm_control_cell_reflectance, thatcontrolmax(j) );    


    i=1;
    while i<= length( ref_control_times{j} )

        % Remove times from both stim and control that are NaN, or 0
        if  isnan(ref_control_times{j}(i)) || ref_control_times{j}(i) ==0 || isnan(ref_variance_control{j}(i)) || ref_variance_control{j}(i) ==0

            ref_control_count{j}(i) = [];
            ref_control_times{j}(i) = [];
            ref_variance_control{j}(i) = [];        
        else
            i = i+1;
        end

    end
    precontrol = ref_variance_control{j}(ref_control_times{j}<=66)./ref_control_count{j}(ref_control_times{j}<=66);
    figure(8); plot(ref_control_times{j}/16.6, sqrt(ref_variance_control{j}./ref_control_count{j})-sqrt(mean(precontrol)) ); hold on; drawnow;      
end
hold off;axis([0 max(thatcontrolmax)/16.6 0 4])

pooled_variance_control = zeros(allmax,1);
pooled_variance_control_count = zeros(allmax,1);
pooled_stddev_control_times = 0:allmax;

% Create the pooled variance for the control data
for j=1:length(profileCDataNames) 
    for i=1:length(ref_control_times{j})
    
        % Create the upper and lower halves of our pooled variance
        pooled_variance_control( ref_control_times{j}(i) ) = pooled_variance_control( ref_control_times{j}(i) ) + ref_variance_control{j}(i);
        pooled_variance_control_count( ref_control_times{j}(i) ) = pooled_variance_control_count( ref_control_times{j}(i) ) + (ref_control_count{j}(i)-1);

    end
end

for i=1:length(pooled_variance_control)    
    pooled_variance_control(i) = pooled_variance_control(i)/pooled_variance_control_count(i);
end

pooled_stddev_control_times( isnan(pooled_variance_control(i)) ) = NaN;


for j=1:length(profileSDataNames)

    profileSDataNames{j}
    load(fullfile(stimRootDir,profileSDataNames{j}));
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance) );
    stim_cell_times            = cell_times(  ~cellfun(@isempty,cell_times) );

    num_stim_cones = num_stim_cones+length(stim_cell_times);
    
    thatstimmax(j) = max( cellfun(@max,stim_cell_times) );    
    
    % Pooled variance of all cells before first stimulus
    [ ref_stddev_stim{j}, ref_stim_times{j} ] = reflectance_std_dev( stim_cell_times( ~cellfun(@isempty,stim_cell_times) ), ...
                                                                       norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) ), thatstimmax(j) );

    i=1;
    while i<= length( ref_stim_times{j} )

        % Remove times from both stim and control that are NaN, or 0
        if isnan(ref_stim_times{j}(i)) || ref_stim_times{j}(i) == 0 || isnan(ref_stddev_stim{j}(i)) || ref_stddev_stim{j}(i) == 0

            ref_stim_times{j}(i) = [];

            ref_stddev_stim{j}(i) = [];
      
        else
            i = i+1;
        end

    end
    
    figure(9); plot(ref_stim_times{j}/16.6, ref_stddev_stim{j} ); hold on; drawnow;
end
hold off;axis([0 max(thatstimmax)/16.6 0 4])

if (length(stim_cell_times) + length(control_cell_times)) < min_cones
    min_cones = (length(stim_cell_times) + length(control_cell_times));
end
if (length(stim_cell_times) + length(control_cell_times)) > max_cones
    max_cones = (length(stim_cell_times) + length(control_cell_times));
end

pooled_stddev_control = sqrt(pooled_variance_control);
stddev_train_times = nan(1, sum(thatstimmax));
stddev_train = nan(1, sum(thatstimmax));
stdiz_breaks = zeros(length(ref_stddev_stim)*2,2 );

% Subtract our average control signal from each of the stimulus trials
offset =0;
for j=1:length(ref_stddev_stim)
    for k=1:length(ref_stim_times{j})       
        if ~isnan( pooled_stddev_control_times(ref_stim_times{j}(k)) )
            stddev_train_times(k+offset) = pooled_stddev_control_times(ref_stim_times{j}(k))+offset;
            stddev_train(k+offset) = ref_stddev_stim{j}(k)-pooled_stddev_control(k);
        end        
    end
    
    stdiz_breaks( (2*j-1):(2*j), 1 ) = [offset; offset];
    if mod(j,2) == 0
        stdiz_breaks( (2*j-1):(2*j), 2 ) = [-4; 4];
    else
        stdiz_breaks( (2*j-1):(2*j), 2 ) = [4; -4];
    end
    
    offset = offset+max(ref_stim_times{j});
end

figure(10); plot(stddev_train_times/16.6,stddev_train,'b'); hold on;
plot(stdiz_breaks(:,1)/16.6,stdiz_breaks(:,2),'k');

plot((stdiz_breaks(1:2:end,1)+66)/16.6, repmat(3,[1 length(stdiz_breaks(1:2:end,1))]) ,'r*');
axis([0 max(stddev_train_times/16.6) -0.5 3.5]) 


% For structure: /stuff/id/intensity/time/data
[remain kid] = getparent(stimRootDir);

% [remain region] = getparent(remain);
[remain stim_time] = getparent(remain);
[remain stim_intensity] = getparent(remain);
[remain stimwave] = getparent(remain);
[~, id] = getparent(remain);

title(['Head to Tail analysis: ' id ',' stim_time ',' stim_intensity ''])

outFname = [id '_' stimwave '_' stim_intensity '_' stim_time '_head_to_tail_' num2str(length(profileSDataNames)) '_signals_twosource'];

figure(8);  title('All control signals'); xlabel('Frame #'); ylabel('Standard deviations'); %axis([0 249 -20 75]);
saveas(gcf, fullfile(pwd, [outFname '_allcontrol.png']), 'png' );
figure(9);  title('All stimulus signals'); xlabel('Frame #'); ylabel('Standard deviations'); %axis([0 249 -20 75]);
saveas(gcf, fullfile(pwd, [outFname '_allstim.png']), 'png' );

hz=16.66666666;
timeBase = ((1:allmax)/hz)';

% dlmwrite(fullfile(pwd, [outFname '.csv']), [timeBase sqrt(pooled_variance_stim) sqrt(pooled_variance_control)], ',' );


figure(10);
saveas(gcf, fullfile(pwd, [outFname '.png']), 'png' );
saveas(gcf, fullfile(pwd, [outFname '.svg']), 'svg' );

% save( fullfile(pwd,['pooled_var_aggregate_' num2str(length(profileDataNames)) '_signals.mat' ] ), 'pooled_std_stim', 'timeBase' );

% save thisshit.mat
% [fitCharacteristics, residuals] = modelFit(timeBase, pooled_std_stim);
% figure(2); hold on;
% plot(trainlocs, (.2+max(pooled_std_stim))*ones(size(trainlocs)),'y*'); hold off;
% 
% % saveas(gcf, fullfile(pwd, [outFname '_wfit.png']) );
% saveas(gcf, fullfile(pwd, [outFname '_wfit.fig']) );
% saveas(gcf, fullfile(pwd, [outFname '_wfit.svg']), 'svg' );
% % figure(1);
% % saveas(gcf, fullfile(pwd, [outFname '_meanratio.svg']), 'svg' );
% % close(8);
% 
% fitCharacteristics.min_cones = min_cones;
% fitCharacteristics.max_cones = max_cones;
% fitCharacteristics.avg_num_cones = num_control_cones/length(profileSDataNames);
% fitCharacteristics.num_pooled = length(profileSDataNames);
% fitCharacteristics.subject = id;
% fitCharacteristics.stim_intensity = stim_intensity;
% fitCharacteristics.stim_length = stimlen;
% fitCharacteristics.stim_wavelength = stimwave;
% fitCharacteristics

