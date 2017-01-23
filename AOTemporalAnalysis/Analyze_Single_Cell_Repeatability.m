
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
        
        num_stim_pts = length(stim_times);
        pos_interval_bound =  tinv(0.95, num_stim_pts-1)*ref_stddev_control{j}( stim_frames )*sqrt(1+(1/num_stim_pts));
        neg_interval_bound = -tinv(0.95, num_stim_pts-1)*ref_stddev_control{j}( stim_frames )*sqrt(1+(1/num_stim_pts));

        stim_over  = sum( stim_reflectance > pos_interval_bound )/num_stim_pts;
        stim_under = sum( stim_reflectance < neg_interval_bound )/num_stim_pts;
        
%         interval_bound = prestim_mean + 1.96*prestim_std/sqrt(num_prestim_pts);                               % Confidence interval

        allstim = stim_under + stim_under;

        if allstim > 0
            control_stim_response_detected(i,j) = 1;
%             title('Response detected');

        else
%             title('Response NOT detected');
        end

        figure(2); plot( stimulus_cones_times{j}{i},  reflectance );  hold on;
        plot( stim_times,  pos_interval_bound,'b' );
        plot( stim_times,  neg_interval_bound,'b' );
        
        maxresp = max(maxresp, max(reflectance) );
        maxtime = max(maxtime, length(stimulus_cones_times{j}{i}) );
        
%         plot( control_times{j}{i}, control_reflectance{j}{i}); hold on;
%         maxresp = max(maxresp, max(control_reflectance{j}{i}) );
%         maxtime = max(maxtime, length(control_times{j}{i}) );
    end
    
    percent_resp_control(i)  = 100*( sum(control_stim_response_detected(i,:))/length(profileDataNames) );
    
     
    plot( stim_times, maxresp*ones(num_stim_pts,1 ),'r*'); title([num2str(percent_resp_control(i)) '% over threshold.']);  hold off;
%     frame = getframe(gcf);
%     writeVideo(profile_vid,frame);

end
figure(4); hist(percent_resp_control,6); title('Histogram of cell response repeatability- control standard devs'); ylabel('# of cells'); xlabel('% times reponded to stimulus');


%% Determine % response using the prediction interval of each cell as the threshold.

% profile_vid = VideoWriter('NC_11049_0ND_control_cell_profiles.avi','Uncompressed AVI');
% open(profile_vid);

% pred_int_stim_response_detected = zeros( length(cons_stimulus_inds), length(profileDataNames) );
% percent_resp = zeros( size(pred_int_stim_response_detected,1),1 );
% 
% for i=1: length(cons_stimulus_inds)
%     
%     maxresp = 0;
%     maxtime = 0;
%     power_spect=[];
%     
%     for j=1:length(profileDataNames)
% %         conv(stim_reflectance{j}{i}, filt_coeff,'same')
% 
%         reflectance = (conv(stimulus_cones_reflectance{j}{i}, filt_coeff,'same') );
% 
%         % Determine if the cell responded in this trial.
%         prestim_frames = stimulus_cones_times{j}{i}( stimulus_cones_times{j}{i}<stim_locs(1) & ~isnan( reflectance ) )';
%         prestim_reflectance = reflectance( stimulus_cones_times{j}{i}<stim_locs(1) & ~isnan( reflectance ) )';
%         num_prestim_pts = length( prestim_reflectance );
%         
%         prestim_reg_coef  = [ones(length(prestim_frames),1) prestim_frames]\prestim_reflectance;
%         prestim_regression= prestim_reg_coef(1) + prestim_frames*prestim_reg_coef(2);
%         
%         prestim_mean =   mean( prestim_reflectance );
% %         prestim_std  =   std( prestim_reflectance-prestim_regression );
%         prestim_std  =   std( prestim_reflectance );
%         
%         SS = (length(prestim_reflectance)-1) * var(prestim_reflectance);
%         resid = norm(prestim_reflectance-prestim_regression);
%         R_squared = 1-resid.^2 / SS;
%         
%                 
%         stim_frames = stimulus_cones_times{j}{i}( stimulus_cones_times{j}{i}>=stim_locs(1) & stimulus_cones_times{j}{i}<=stim_locs(end) & ~isnan( reflectance ) );
%         stim_reflectance = reflectance( stimulus_cones_times{j}{i}>=stim_locs(1) & stimulus_cones_times{j}{i}<=stim_locs(end) & ~isnan( reflectance ) );
%         num_stim_pts = length(stim_reflectance);
%         stim_regression   = prestim_reg_coef(1) + stim_frames*prestim_reg_coef(2);
% 
%         pos_interval_bound = stim_regression + tinv(0.99, num_prestim_pts-1)*prestim_std*sqrt(1+(1/num_prestim_pts));  % 99% Prediction interval
%         neg_interval_bound = stim_regression - tinv(0.99, num_prestim_pts-1)*prestim_std*sqrt(1+(1/num_prestim_pts));
%         
%         stim_over  = sum( stim_reflectance > pos_interval_bound )/num_stim_pts;
%         stim_under = sum( stim_reflectance < neg_interval_bound )/num_stim_pts;
%         
% %         interval_bound = prestim_mean + 1.96*prestim_std/sqrt(num_prestim_pts);                               % Confidence interval
% 
%         figure(2); hold on;
%         
%         allstim = stim_under + stim_under;
% 
%         if allstim >= .05
%             pred_int_stim_response_detected(i,j) = 1;
% %             figure(2); plot( stimulus_cones_times{j}{i},  reflectance ); hold on;
% %             plot( stim_frames,  stim_reflectance );
% %             plot( [prestim_frames' stim_frames], [prestim_regression' stim_regression] );
% %             plot(stim_frames, pos_interval_bound);
% %             plot(stim_frames, neg_interval_bound); hold off;
% %             R_squared
%             title('Response detected');
% %         elseif stim_under >= .05
% %             pred_int_stim_response_detected(i,j) = 1; 
% %             title('Response detected');
%         else
%             title('Response NOT detected');
%         end
% 
%         plot( stimulus_cones_times{j}{i},  reflectance ); 
%         plot( stim_frames,  stim_reflectance );
%         plot( [prestim_frames' stim_frames], [prestim_regression' stim_regression] );
%         plot(stim_frames, pos_interval_bound);
%         plot(stim_frames, neg_interval_bound); hold off;
% 
% 
% 
%         
%         maxresp = max(maxresp, max(reflectance) );
%         maxtime = max(maxtime, length(stimulus_cones_times{j}{i}) );
%         
% %         plot( stimulus_cones_times{j}{i}, reflectance); hold on;
% %         maxresp = max(maxresp, max(control_reflectance{j}{i}) );
% %         maxtime = max(maxtime, length(control_times{j}{i}) );
%     end
%     
%     percent_resp(i)  = 100*( sum(pred_int_stim_response_detected(i,:))/length(profileDataNames) );
%     
%     plot(stim_locs, maxresp*ones(34),'r*'); title([num2str(percent_resp(i)) '% over threshold.']);  hold off;
% %     frame = getframe(gcf);
% %     writeVideo(profile_vid,frame);
% percent_resp;
% end
% % close(profile_vid);
% 
% figure(3); hist(percent_resp,6); title('Histogram of cell response repeatability- prediction interval'); ylabel('# of cells'); xlabel('% times reponded to stimulus');
