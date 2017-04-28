function [characteristics]=Aggregate_Multiple_Temporal_Analyses_bootstrap(rootDir)
% Robert F Cooper
% 12-31-2015
% This script calculates pooled variance across a set of given signals.

if ~exist('rootDir','var')
    close all force;
    rootDir = uigetdir(pwd);
end

profileDataNames = read_folder_contents(rootDir,'mat');

thatstimmax=0;
thatcontrolmax=0;
%% Code for determining variance across all signals at given timepoint

allmax=0;
num_control_cones = 0;
num_stim_cones = 0;
max_cones = 0;
min_cones = 10000000000000;
mean_control_reflectance = zeros(500,1);

for j=1:length(profileDataNames)
    clear mean_ratio;
    profileDataNames{j}
    load(fullfile(rootDir,profileDataNames{j}));
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
    control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );

    if exist('mean_ratio','var')
        
        mean_ratio_times = unique(cell2mat(control_cell_times'));
        all_ratio_times{j} = mean_ratio_times;    
        all_mean_ratio{j} = mean_ratio;
    end
    
    num_control_cones = num_control_cones+length(control_cell_times);
    num_stim_cones = num_stim_cones+length(stim_cell_times);
    
    if (length(stim_cell_times) + length(control_cell_times)) < min_cones
        min_cones = (length(stim_cell_times) + length(control_cell_times));
    end
    if (length(stim_cell_times) + length(control_cell_times)) > max_cones
        max_cones = (length(stim_cell_times) + length(control_cell_times));
    end
    
    thatstimmax = max( cellfun(@max,stim_cell_times) );    
    
    thatcontrolmax = max( cellfun(@max,control_cell_times) );

    if thatstimmax ~= thatcontrolmax
        error('The control and stimulus number of frames do not match! Cannot perform analysis...');
    end
    if thatstimmax > allmax
        allmax=thatstimmax;
    end

    
    % Pooled variance of all cells before first stimulus
    [ ref_variance_stim{j}, ref_stim_times{j}, ref_stim_count{j} ]    = reflectance_pooled_variance( stim_cell_times, norm_stim_cell_reflectance, allmax );
    [ ref_variance_control{j},ref_control_times{j}, ref_control_count{j} ] = reflectance_pooled_variance( control_cell_times, norm_control_cell_reflectance, allmax );    

    i=1;
    while i<= length( ref_stim_times{j} )

        % Remove times from both stim and control that are NaN, or 0
        if isnan(ref_stim_times{j}(i)) || ref_stim_times{j}(i) == 0 || ...
           isnan(ref_variance_stim{j}(i)) || ref_variance_stim{j}(i) == 0

            ref_stim_count{j}(i) = [];
            
            ref_stim_times{j}(i) = [];

            ref_variance_stim{j}(i) = [];       
        else
%             ref_times = [ref_times; ref_stim_times(i)];
            i = i+1;
        end

    end
    
    i=1;
    while i<= length( ref_control_times{j} )

        % Remove times from both stim and control that are NaN, or 0
        if isnan(ref_control_times{j}(i)) || ref_control_times{j}(i) ==0 || ...
           isnan(ref_variance_control{j}(i)) || ref_variance_control{j}(i) ==0

            ref_control_count{j}(i) = [];
            
            ref_control_times{j}(i) = [];

            ref_variance_control{j}(i) = [];        
        else
%             ref_times = [ref_times; ref_stim_times(i)];
            i = i+1;
        end

    end
    precontrol = ref_variance_control{j}(ref_control_times{j}<=66)./ref_control_count{j}(ref_control_times{j}<=66);
    prestim = ref_variance_stim{j}(ref_stim_times{j}<=66)./ref_stim_count{j}(ref_stim_times{j}<=66);
    
    figure(8); plot(ref_stim_times{j}, sqrt( ref_variance_stim{j}./ref_stim_count{j})-sqrt(mean(prestim)) ); hold on; drawnow;
    figure(9); plot(ref_control_times{j}, sqrt(ref_variance_control{j}./ref_control_count{j})-sqrt(mean(precontrol)) ); hold on; drawnow;
      
    for i=1 : length(norm_control_cell_reflectance)
        for k=1 : length( norm_control_cell_reflectance{i} )

            if ~isnan( norm_control_cell_reflectance{i}(k) )
                if mean_control_reflectance(k) == 0
                    mean_control_reflectance(k) = norm_control_cell_reflectance{i}(k);
                else
                    mean_control_reflectance(k) = (mean_control_reflectance(k) + norm_control_cell_reflectance{i}( (k) ) )/2;
                end
            end
        end
    end
      

end
figure(8); hold off;
figure(9); hold off;

hz=16.66666666;
timeBase = ((1:allmax)/hz)';

[remain kid] = getparent(rootDir);
[remain region] = getparent(remain);
[remain stim_time] = getparent(remain);
[remain stim_intensity] = getparent(remain);
[remain stimwave] = getparent(remain);
[~, id] = getparent(remain);

stimlen = str2double( strrep(stim_time(1:3),'p','.') );

%% Bootstrap using the above signals
parfor b=1:1500

    rng('shuffle'); % Reshuffle the RNG after each loop to make sure we're getting a good mix.
    
    num_signals_stim = cell(allmax, 1);
    num_signals_control = cell(allmax, 1);
    
    % Get an array of which signals have data for each time point    
    for j=1:length(profileDataNames)
        for i=1:length(ref_stim_times{j})
            num_signals_stim{ ref_stim_times{j}(i) } = [num_signals_stim{ ref_stim_times{j}(i) }; j ];
        end
    end
    
    for j=1:length(profileDataNames) 
        for i=1:length(ref_control_times{j})
            num_signals_control { ref_control_times{j}(i) } = [num_signals_control{ ref_control_times{j}(i) } j];
        end
    end
    
    stim_signal_inds = cell(1,length(num_signals_stim));
    control_signal_inds = cell(1,length(num_signals_control));
    
    % Randomly pick the datapoints we'll average from the signals we have
    % at each time point
    for i=1:length(num_signals_stim)
        if ~isempty(num_signals_stim{i})
            signal_picks = randi( length(num_signals_stim{i}), length(num_signals_stim{i}),1);
            stim_signal_inds{i} = num_signals_stim{i}( signal_picks );
        end
    end
    
    for i=1:length(num_signals_control)
        if ~isempty(num_signals_control{i})
            signal_picks = randi( length(num_signals_control{i}), length(num_signals_control{i}),1);
            control_signal_inds{i} = num_signals_control{i}( signal_picks );
        end
    end       

    
    pooled_variance_stim = zeros(allmax, 1);
    pooled_variance_stim_count = zeros(allmax, 1);
    pooled_variance_control = zeros(allmax, 1);
    pooled_variance_control_count = zeros(allmax, 1);

    %% Create the pooled variance for each of these
    for i=1:length(pooled_variance_stim)
        for j=1:length(stim_signal_inds{i})
            which_ind = stim_signal_inds{i}(j);
            which_time = ref_stim_times{ which_ind };
            which_signal = ref_variance_stim{ which_ind };
            which_count = ref_stim_count{ which_ind };
            
            % Create the upper and lower halves of our pooled variance
            pooled_variance_stim( i ) = pooled_variance_stim( i ) + which_signal( which_time==i );
            pooled_variance_stim_count( i ) = pooled_variance_stim_count( i ) + which_count( which_time==i );
        end
    end

    for i=1:length(pooled_variance_control)
        for j=1:length(control_signal_inds{i})
            which_ind = control_signal_inds{i}(j);
            which_time = ref_control_times{ which_ind };
            which_signal = ref_variance_control{ which_ind };
            which_count = ref_control_count{ which_ind };
            
            % Create the upper and lower halves of our pooled variance
            pooled_variance_control( i ) = pooled_variance_control( i ) + which_signal( which_time==i );
            pooled_variance_control_count( i ) = pooled_variance_control_count( i ) + which_count( which_time==i );
        end
    end

    for i=1:length(pooled_variance_stim)    
        pooled_variance_stim(i) = pooled_variance_stim(i)/pooled_variance_stim_count(i);
    end
    for i=1:length(pooled_variance_control)    
        pooled_variance_control(i) = pooled_variance_control(i)/pooled_variance_control_count(i);
    end

%     figure(10); 
%     plot( timeBase,sqrt(pooled_variance_stim),'r'); hold on;
%     plot( timeBase,sqrt(pooled_variance_control),'b');

    pooled_std_stim    = sqrt(pooled_variance_stim)-sqrt(pooled_variance_control);

%     plot( timeBase(~isnan(pooled_std_stim)), pooled_std_stim(~isnan(pooled_std_stim)),'k'); hold on;
%     legend('Stimulus cones','Control cones','Subtraction');


    
    [fitCharacteristics, residuals] = modelFit(timeBase, pooled_std_stim);    

    if (residuals < 10) && (fitCharacteristics.amplitude < 3)
       all_amps(b) = fitCharacteristics.amplitude;
       all_res(b) = residuals; 
    else
       all_amps(b) = -1;
       all_res(b) = -1;
    end
end

all_amps = all_amps(all_amps>=0);
all_res = all_res(all_res>=0);

all_amps = all_amps(1:1000);
all_res = all_res(1:1000);

characteristics.avg_num_control_cones = num_control_cones/length(profileDataNames);
characteristics.avg_num_stim_cones = num_stim_cones/length(profileDataNames);
characteristics.num_control_pooled = length(profileDataNames);
characteristics.num_stim_pooled = length(profileDataNames);
characteristics.subject = id;
characteristics.stim_intensity = stim_intensity;
characteristics.stim_length = stimlen;
characteristics.stim_wavelength = stimwave;

characteristics.all_amps = all_amps;
characteristics.mean_amp = mean(all_amps);
characteristics.std_amp = std(all_amps);

characteristics.mean_mse = mean(all_res);
characteristics.std_mse =std(all_res);
characteristics.min_mse = min(all_res);
characteristics.max_mse = max(all_res);
characteristics
