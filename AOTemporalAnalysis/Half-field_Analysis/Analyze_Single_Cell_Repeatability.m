
% Change normalization, then analyze


clear;
close all force;


profileDataNames = read_folder_contents(pwd,'mat');

if length(profileDataNames) < 2
    error('Requires more than one dataset to do conservation analysis.');
end


load(profileDataNames{1});

control_inds = cell(length(profileDataNames),1);
stimulus_inds = cell(length(profileDataNames),1);

control_times = cell(length(profileDataNames),1);
control_reflectance = cell(length(profileDataNames),1);

stimulus_cones_times = cell(length(profileDataNames),1);
stimulus_cones_reflectance = cell(length(profileDataNames),1);

for j=1:length(profileDataNames)

    load(profileDataNames{j});
    
    control_inds{j} = contcellinds;       
    control_times{j} = control_cell_times;
    control_reflectance{j} = norm_control_cell_reflectance;
    
    stimulus_inds{j} = stimcellinds;
    stimulus_cones_times{j} = stim_cell_times;
    stimulus_cones_reflectance{j} = norm_stim_cell_reflectance;
        stimulus_inds{j}
end


[cons_control_inds]= intersect(control_inds{1}, control_inds{2});
[cons_stimulus_inds]= intersect(stimulus_inds{1}, stimulus_inds{2});

for i=3:1:length(profileDataNames)
   
    [cons_control_inds] = intersect(cons_control_inds, control_inds{i});
    [cons_stimulus_inds] = intersect(cons_stimulus_inds, stimulus_inds{i});
    
end

% stim_cone_mean_fft = fft(control_reflectance{1}(cons_control_inds));

for j=1:length(profileDataNames)

    control_cones_times{j} = control_times{j}(cons_control_inds);
    control_cones_reflectance{j} = control_reflectance{j}(cons_control_inds);
    
    thatcontrolmax = max( cellfun(@max, control_cones_times{j}( ~cellfun(@isempty, control_cones_times{j} )) ) );  
    [ ref_stddev_control{j}, ref_control_times{j} ] = reflectance_std_dev( control_cones_times{j}( ~cellfun(@isempty,control_cones_times{j}) ), ...
                                                                           control_cones_reflectance{j}( ~cellfun(@isempty,control_cones_reflectance{j}) ), thatcontrolmax );
%     stim_cone_mean_fft = mean([stim_cone_mean_fft fft(stimulus_cones_reflectance{j})]);
%     plot(isnan(ref_stddev_control{j}))
    stimulus_cones_times{j} = stimulus_cones_times{j}(cons_stimulus_inds);
    stimulus_cones_reflectance{j} = stimulus_cones_reflectance{j}(cons_stimulus_inds);
    
end

hz = 16.666666;
stim_locs = 68:(68+33);

avg_order = 1;
filt_coeff = ones(avg_order,1)/avg_order;

%% Determine % response using the control standard deviation as the
%% threshold.

% Determine the control pooled variance to use as a threshold.
% Aggregate_Multiple_Temporal_Analyses;


    

control_stim_response_detected = zeros( length(cons_stimulus_inds), length(profileDataNames) );
percent_resp_control = zeros( length(cons_stimulus_inds),1 );

for i=1: length(cons_stimulus_inds)
    
    maxresp = 0;
    maxtime = 0;
    power_spect=[];
    

    for j=1:length(profileDataNames)

        reflectance = stimulus_cones_reflectance{j}{i}; %conv(stimulus_cones_reflectance{j}{i}, filt_coeff,'same');
        stim_frames = stimulus_cones_times{j}{i}>=stim_locs(1) & stimulus_cones_times{j}{i}<=stim_locs(end) & ~isnan( reflectance );
        stim_times  = stimulus_cones_times{j}{i}(stim_frames);
        % Determine if the cell responded in this trial.
        stim_reflectance  = reflectance(stim_frames)';
        

        figure(2); plot( stimulus_cones_times{j}{i},  reflectance );  hold on;
        

%         plot( control_times{j}{i}, control_reflectance{j}{i}); hold on;
%         maxresp = max(maxresp, max(control_reflectance{j}{i}) );
%         maxtime = max(maxtime, length(control_times{j}{i}) );
    end
    hold off;
    saveas(gcf, fullfile(pwd, ['repeated_cone_profile' num2str(i) '.png']), 'png' );
end
