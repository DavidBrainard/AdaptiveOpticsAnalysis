function [fitCharacteristics residuals]=Aggregate_Multiple_Temporal_Analyses(rootDir)
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
num_cones = 0;
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
    
    num_cones = num_cones+length(stim_cell_times) + length(control_cell_times);
    
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
    while i<= length( ref_control_times{j} )

        % Remove times from both stim and control that are NaN, or 0
        if isnan(ref_stim_times{j}(i)) || isnan(ref_control_times{j}(i)) || ref_stim_times{j}(i) == 0 || ref_control_times{j}(i) ==0 || ...
           isnan(ref_variance_stim{j}(i)) || isnan(ref_variance_control{j}(i)) || ref_variance_stim{j}(i) == 0 || ref_variance_control{j}(i) ==0

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
    
    figure(8); plot(ref_stim_times{j}, sqrt(ref_variance_stim{j})-sqrt(ref_variance_stim{j}(1)) ); hold on; drawnow;
    figure(9); plot(ref_control_times{j}, sqrt(ref_variance_control{j})-sqrt(ref_variance_control{j}(1)) ); hold on; drawnow;
      
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

% For structure: /stuff/id/intensity/time/region_cropped/data
[remain kid] = getparent(rootDir);
% [remain kid] = getparent(remain);
[remain region] = getparent(remain);
[remain stim_time] = getparent(remain);
[remain stim_intensity] = getparent(remain);
[remain stimwave] = getparent(remain);
[~, id] = getparent(remain);

outFname = [id '_' stimwave '_' stim_intensity '_' stim_time '_aggregate_' num2str(length(profileDataNames)) '_signals'];

hz=16.66666666;
timeBase = ((1:allmax)/hz)';

%% Fitting an exp to the curve
% options = optimset('fmincon');
% options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','interior-point');
% 
% std_control = sqrt(pooled_variance_control);
% 
% timeBase = timeBase(~isnan(std_control) );
% std_control = std_control(~isnan(std_control));
% 
% % Exponential
% % parameters(1) = std_control(1);
% % parameters(2) = 1/(max(std_control) - min(std_control));
% % parameters(3) = 2;
% % parameters(4) = 1;
% 
% %Logistic
% parameters(1) = mean(std_control(1:66));
% parameters(2) = max(std_control) - min(std_control);
% parameters(3) = 1;
% parameters(4) = timeBase(124);
% 
% vlb = [-1 0               0 0];
% vub = [ 1 max(std_control) 2 15];
% 
% cutoutTime = timeBase([find(timeBase<66) find(timeBase>150)]);
% cutoutControl = std_control([find(timeBase<66) find(timeBase>150)]);
% 
% out_params = fmincon(@(parameters)control_fit(parameters, cutoutControl, cutoutTime), parameters,[],[],[],[],vlb,vub,[],options);
% out_params
% % prediction = out_params(1) + out_params(2)*out_params(3).^( out_params(4) * timeBase );
% prediction = out_params(1) + ( out_params(2) ./ ( 1+exp(-out_params(3).*(timeBase-out_params(4)) ) ) );
% 
% plot(cutoutTime, cutoutControl, timeBase, prediction );


%%
dlmwrite(fullfile(pwd, [outFname '.csv']), [timeBase sqrt(pooled_variance_stim) sqrt(pooled_variance_control)], ',' );


pooled_std_stim    = sqrt(pooled_variance_stim)-sqrt(pooled_variance_control);

% Stim train
stimlen = str2double( strrep(stim_time(1:3),'p','.') );

%% Manual marking
% figure(7); plot(timeBase, sqrt(pooled_variance_stim),'r'); hold on;
trainlocs = 68/hz:1/hz:(68/hz+stimlen);
% plot(timeBase, sqrt(pooled_variance_control),'b');
% plot(trainlocs, max(pooled_std_stim)*ones(size(trainlocs)),'r*'); hold off;
% [timepts, heightpts]=ginput(2);


figure(10); 
plot( timeBase,pooled_std_stim,'r'); hold on;
% plot( timeBase,pooled_std_control,'b');
% legend('Stimulus cones','Control cones');

trainlocs = 66/hz:1/hz:(66/hz+stimlen);
plot(trainlocs, max(pooled_std_stim)*ones(size(trainlocs)),'r*'); hold off;

% plot(stim_locs, max([ref_variance_stim; ref_variance_control])*ones(size(stim_locs)),'r*'); hold off;
ylabel('Pooled Standard deviation'); xlabel('Time (s)'); title( [stim_intensity ' ' stim_time 'pooled standard deviation of ' num2str(length(profileDataNames)) ' signals.'] );
axis([0 16 -1 3])
hold off;
saveas(gcf, fullfile(pwd, [outFname '.png']), 'png' );
% saveas(gcf, fullfile(pwd, [outFname '.svg']), 'svg' );

% save( fullfile(pwd,['pooled_var_aggregate_' num2str(length(profileDataNames)) '_signals.mat' ] ), 'pooled_std_stim', 'timeBase' );

dlmwrite(fullfile(pwd, [date '_all_plots.csv']), [ [str2double(id(4:end)), str2double(stim_intensity(1:3)), stimlen] ;[ timeBase sqrt(pooled_variance_stim) sqrt(pooled_variance_control) ] ]',...
         '-append', 'delimiter', ',', 'roffset',1);

% save thisshit.mat
[fitCharacteristics, residuals] = modelFit(timeBase, pooled_std_stim);
figure(2); hold on;
plot(trainlocs, (.2+max(pooled_std_stim))*ones(size(trainlocs)),'y*'); hold off;

saveas(gcf, fullfile(pwd, [outFname '_wfit.png']) );
saveas(gcf, fullfile(pwd, [outFname '_wfit.fig']) );
saveas(gcf, fullfile(pwd, [outFname '_wfit.svg']), 'svg' );
% figure(1);
% saveas(gcf, fullfile(pwd, [outFname '_meanratio.svg']), 'svg' );
% close(8);
if exist('heightpts','var')
    fitCharacteristics.absolute_height = abs(heightpts(2)-heightpts(1));
end

fitCharacteristics.min_cones = min_cones;
fitCharacteristics.max_cones = max_cones;
fitCharacteristics.avg_num_cones = num_cones/length(profileDataNames);
fitCharacteristics.num_pooled = length(profileDataNames);
fitCharacteristics.subject = id;
fitCharacteristics.stim_intensity = stim_intensity;
fitCharacteristics.stim_length = stimlen;
fitCharacteristics.stim_wavelength = stimwave;

%% Mean ratio analyses
% if exist('mean_ratio','var') && length(profileDataNames) == length(all_ratio_times)
%     maxpts = cellfun(@max, all_ratio_times);
%     maxlen = max(maxpts);
%     mean_mean_ratios = nan(1,maxlen);
% 
%     for i=1:length(all_ratio_times)
%         for j=1:length(all_ratio_times{i})
% 
%             if isnan( mean_mean_ratios( all_ratio_times{i}(j) ) )            
%                 mean_mean_ratios( all_ratio_times{i}(j) ) = all_mean_ratio{i}(j);
%             else
%                 mean_mean_ratios( all_ratio_times{i}(j) ) = ( mean_mean_ratios( all_ratio_times{i}(j) ) + all_mean_ratio{i}(j) )/2;
%             end
%         end
%         before(i) = mean(all_mean_ratio{i}(1:66));
%         during(i) = mean(all_mean_ratio{i}( uint8(trainlocs*hz) ));
%         after(i) = mean(all_mean_ratio{i}( uint8(trainlocs(end)*hz):end));
%     end
% %     figure(1); plot(mean_mean_ratios);
% 
%     minvals = min([before during after]);
%     maxvals = max([before during after]);
%     figure(1); plot(before,during,'r*',before,after,'b*' );
%     dlmwrite( [stim_intensity '_' stim_time '_beforeduringafter.csv'],[before' during' after'],'-append','delimiter',',');
% end


