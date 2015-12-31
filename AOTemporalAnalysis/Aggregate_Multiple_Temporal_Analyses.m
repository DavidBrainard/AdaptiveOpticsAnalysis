
clear;
close all force;



profileDataNames = read_folder_contents(pwd,'mat');

thatstimmax=0;
thatcontrolmax=0;
%% Code for determining variance across all signals at given timepoint
% for j=1:length(profileDataNames)
% 
%     load(profileDataNames{j});
%     
%     % Remove the empty cells
%     norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
%     stim_cell_times = stim_cell_times( ~cellfun(@isempty,stim_cell_times) );
%     norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance) );
%     control_cell_times = control_cell_times( ~cellfun(@isempty,control_cell_times) );
%     
%     
%     thismax = max( cellfun(@max,stim_cell_times) );
%     if thismax > thatstimmax
%        thatstimmax = thismax; 
%     end
%     
%     thismax = max( cellfun(@max,control_cell_times) );
%     if thismax > thatcontrolmax
%        thatcontrolmax = thismax; 
%     end
% end
% 
% all_stim_ref = cell(thatstimmax,1);
% 
% all_control_ref = cell(thatcontrolmax,1);
% 
% % Consolidate all of the data
% for j=1:length(profileDataNames)
% 
%     load(profileDataNames{j});
%     
%     % Remove the empty cells
%     norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
%     stim_cell_times = stim_cell_times( ~cellfun(@isempty,stim_cell_times) );
%     
%     norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance) );
%     control_cell_times = control_cell_times( ~cellfun(@isempty,control_cell_times) );
%     
%     for i=1:length(stim_cell_times)
%         for k=1:length(stim_cell_times{i})
%             all_stim_ref{ stim_cell_times{i}(k) } = [all_stim_ref{ stim_cell_times{i}(k) } norm_stim_cell_reflectance{i}(k)];
%         end
%     end
%     
%     for i=1:length(control_cell_times)
%         for k=1:length(control_cell_times{i})
%             all_control_ref{ control_cell_times{i}(k) } = [all_stim_ref{ control_cell_times{i}(k) } norm_control_cell_reflectance{i}(k)];
%         end
%     end
%     
% end

allmax=0;


for j=1:length(profileDataNames)

    load(profileDataNames{j});
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times = stim_cell_times( ~cellfun(@isempty,stim_cell_times) );
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance) );
    control_cell_times = control_cell_times( ~cellfun(@isempty,control_cell_times) );
    
    
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
    while i<= length( ref_control_times{j} )

        % Remove times from both stim and control that are NaN
        if isnan(ref_stim_times{j}(i)) || isnan(ref_control_times{j}(i))

            ref_stim_count{j}(i) = [];
            ref_control_count{j}(i) = [];
            
            ref_stim_times{j}(i) = [];
            ref_control_times{j}(i) = [];

            ref_variance_stim{j}(i) = [];
            ref_variance_control{j}(i) = [];        
        else
%             ref_times = [ref_times; ref_stim_times(i)];
            i = i+1;
        end

    end
    

end

pooled_variance_stim = zeros(allmax,1);
pooled_variance_stim_count = zeros(allmax,1);

pooled_variance_control = zeros(allmax,1);
pooled_variance_control_count = zeros(allmax,1);

% Create the pooled variance for each of these

for j=1:length(profileDataNames)
    
    for i=1:length(ref_stim_times{j})
    
        % Create the upper and lower halves of our pooled variance
        pooled_variance_stim( ref_stim_times{j}(i) ) = pooled_variance_stim( ref_stim_times{j}(i) ) + ref_variance_stim{j}(i);
        pooled_variance_stim_count( ref_stim_times{j}(i) ) = pooled_variance_stim_count( ref_stim_times{j}(i) ) + (ref_stim_count{j}(i)-1);

    end
    
    for i=1:length(ref_control_times{j})
    
        % Create the upper and lower halves of our pooled variance
        pooled_variance_control( ref_control_times{j}(i) ) = pooled_variance_control( ref_control_times{j}(i) ) + ref_variance_control{j}(i);
        pooled_variance_control_count( ref_control_times{j}(i) ) = pooled_variance_control_count( ref_control_times{j}(i) ) + (ref_control_count{j}(i)-1);

    end
end

for i=1:length(pooled_variance_stim)    
    pooled_variance_stim(i) = pooled_variance_stim(i)/pooled_variance_stim_count(i);
end
for i=1:length(pooled_variance_control)    
    pooled_variance_control(i) = pooled_variance_control(i)/pooled_variance_control_count(i);
end

% If its in the normalization, subtract the control value from the stimulus
% value
% if ~isempty( strfind(norm_type, 'sub') )
        
        pooled_std_stim    = sqrt(pooled_variance_stim)-sqrt(pooled_variance_control);
        pooled_std_control = sqrt(pooled_variance_control)-sqrt(pooled_variance_control);
    
% end
hz=16.6;
figure(10); hold off;

plot( (1:allmax)/hz,pooled_std_stim,'r'); hold on;
plot( (1:allmax)/hz,pooled_std_control,'b'); hold on;
legend('Stimulus cones','Control cones');

% plot(stim_locs, max([ref_variance_stim; ref_variance_control])*ones(size(stim_locs)),'r*'); hold off;
ylabel('Pooled Standard deviation'); xlabel('Time (s)'); title( ['Pooled variance of ' num2str(length(profileDataNames)) ' signals.'] );
saveas(gcf, fullfile(pwd,['pooled_var_aggregate_' num2str(length(profileDataNames)) '_signals.png' ] ) );