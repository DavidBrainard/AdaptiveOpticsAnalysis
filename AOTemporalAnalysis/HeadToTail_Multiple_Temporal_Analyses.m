
% Robert F Cooper
% 12-31-2015
% This script calculates pooled variance across a set of given signals.

clear;
close all force;



profileDataNames = read_folder_contents(pwd,'mat');

thatstimmax=0;
thatcontrolmax=0;
hz = 16.66666666;
%% Code for determining variance across all signals at given timepoint

allmax=0;

mean_control_reflectance = zeros(500,1);

for j=1:length(profileDataNames)
    profileDataNames{j}
    load(profileDataNames{j});
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
    control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );
    
    thatstimmax = max( cellfun(@max,stim_cell_times) );    
    
    thatcontrolmax = max( cellfun(@max,control_cell_times) );

    if thatstimmax ~= thatcontrolmax
        error('The control and stimulus number of frames do not match! Cannot perform analysis...');
    end
    if thatstimmax > allmax
        allmax=thatstimmax;
    end

    
    % Pooled variance of all cells before first stimulus
[ ref_stddev_stim{j}, ref_stim_times{j} ]    = reflectance_std_dev( stim_cell_times( ~cellfun(@isempty,stim_cell_times) ), ...
                                                              norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) ), thatstimmax );
[ ref_stddev_control{j} ,ref_control_times{j} ] = reflectance_std_dev( control_cell_times( ~cellfun(@isempty,control_cell_times) ), ...
                                                                norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance) ), thatcontrolmax );


    i=1;
    while i<= length( ref_control_times{j} )

        % Remove times from both stim and control that are NaN
        if isnan(ref_stim_times{j}(i)) || isnan(ref_control_times{j}(i))

            ref_stim_times{j}(i) = [];
            ref_control_times{j}(i) = [];

            ref_stddev_stim{j}(i) = [];
            ref_stddev_control{j}(i) = [];        
        else
            i = i+1;
        end

    end
    
    ref_stddev_stim_sub{j} = (ref_stddev_stim{j}-ref_stddev_control{j})';
    ref_stim_times{j}=ref_stim_times{j}' + 249*(j-1);% Add an offset so they aren't all on top of each other...
    trainlocs{j} = ((66:1:(65+33)) + 249*(j-1))/hz;

end

ref_stddev_stim_sub = cell2mat(ref_stddev_stim_sub);
ref_stim_times = cell2mat(ref_stim_times);
trainlocs = cell2mat(trainlocs);

figure(10); hold on;
plot( ref_stim_times/hz,ref_stddev_stim_sub,'r');
legend('Stimulus signal');

% Stim train
plot(trainlocs, max(ref_stddev_stim_sub)*ones(size(trainlocs)),'r*'); hold off;

% plot(stim_locs, max([ref_variance_stim; ref_variance_control])*ones(size(stim_locs)),'r*'); hold off;
ylabel('Pooled Standard deviation'); xlabel('Time (s)'); title( ['Pooled standard deviation of ' num2str(length(profileDataNames)) ' signals.'] );
saveas(gcf, fullfile(pwd,['pooled_var_HeadToTail_' num2str(length(profileDataNames)) '_signals.png' ] ) );
% save( fullfile(pwd,['pooled_var_aggregate_' num2str(length(profileDataNames)) '_signals.mat' ] ), 'pooled_std_stim', 'timeBase' );


% brainardbasedModelFit(timeBase, pooled_std_stim)
