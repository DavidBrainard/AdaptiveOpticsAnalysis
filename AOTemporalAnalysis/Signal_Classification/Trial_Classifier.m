% Robert Cooper 2017-08-01
%
% This script uses a previously trained classifier (called either
% randforestclassifier or svmclassifier) and previously processed data
% to create a classification "forest" from an image
clear;
cls_path = 'trained_classifier_11-Aug-2017.mat';

load(cls_path);


% Stimulus

load(['/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/Shutter-based-Data/NC_11049/TRC/550b10nm/660nW/2p0/Profile_Data/NC_11049_20170802_OD_confocal_0013_ref_84_crop_affine_n181_box_global_norm_prestim_stdiz_stimulus_profiledata.mat']);
% load('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/Pathology_Experiments/CHM/CH_13131_20161209/550b10nm/312nW/stimulus/region_cropped/Profile_Data/CH_13131_20161208_OD_confocal_0155_ref_75_crop_affine_box_global_norm_prestim_stdiz_stimulus_profiledata.mat')

% Legacy support
if ~exist('norm_cell_reflectance','var')
    norm_cell_reflectance = norm_stim_cell_reflectance;
    cell_times            = stim_cell_times;
end
norm_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance) );
cell_times            = cell_times(  ~cellfun(@isempty,cell_times) );

norm_cell_reflectance = norm_cell_reflectance(coords_used(:,1)>40);
cell_times = cell_times(coords_used(:,1)>40);
coords_used = coords_used(coords_used(:,1)>40,:);

[coeffs,labels] = extract_features( cell_times, norm_cell_reflectance, [66 99], vid_type );

varNames = {'Affinity','std_derv_pkpk','Max_Deriv_Ampl', 'Max_Ampl', ...
            'Mean_Abs_1D','Mean_Abs_2D','Mean_Abs_3D','Mean_Abs_4D','Mean_Abs_A4',...
            'Mean_Power_1D','Mean_Power_2D','Mean_Power_3D','Mean_Power_4D','Mean_Power_A4',...
            'Std_Dev_1D','Std_Dev_2D','Std_Dev_3D','Std_Dev_4D','Std_Dev_A4',...
            'Coeff_Ratio_1D2D','Coeff_Ratio_2D3D','Coeff_Ratio_3D4D','Coeff_Ratio_4D4A', 'Training_Labels'};
testTable = [array2table(coeffs) cell2table(labels)];
% testTable.Properties.VariableNames = varNames;

[predictions, score] = randforest_classifier.predictFcn(testTable);

confusionmat(labels,predictions)

score(:,2) = sign(score(:,2)).*log( abs(score(:,2))+1);
score(:,2) = (score(:,2)-min(score(:,2)));
score(:,2) = (score(:,2)./max(score(:,2)))*100+1;


cmapped = parula(101);
figure(1); imagesc(ref_image); colormap gray; hold on;
for i=1:length(coords_used)    
    plot( coords_used(i,1), coords_used(i,2),'*','Color',cmapped( ceil( score(i,2) ),:) );
    
%     if i==239
%         strcmp(predictions{i},labels{i})
%         figure(10);plot(cell_times{i},norm_cell_reflectance{i});
%         figure(1);
%     end
end
%     plot( coords_used(strcmp(predictions,'stimulus'),1), coords_used(strcmp(predictions,'stimulus'),2),'r*');
%     plot( coords_used(strcmp(predictions,'control'),1), coords_used(strcmp(predictions,'control'),2),'b*'); hold off;

title('Stimulated cones');

% Control
clear norm_cell_reflectance;
load(['/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/Shutter-based-Data/NC_11049/TRC/550b10nm/control/Profile_Data/NC_11049_20170802_OD_confocal_0005_ref_84_crop_affine_n231_box_global_norm_prestim_stdiz_control_profiledata.mat']);
% load('/local_data/Projects/AOSLO-Intrinsic-Reflectivity-2015_10_01-/data/Pathology_Experiments/CHM/CH_13131_20161209/550b10nm/312nW/control/region_cropped/Profile_Data/CH_13131_20161208_OD_confocal_0146_ref_130_crop_affine_box_global_norm_prestim_stdiz_control_profiledata.mat')

% Legacy support
if ~exist('norm_cell_reflectance','var')
    norm_cell_reflectance = norm_control_cell_reflectance;
    cell_times            = control_cell_times;
end
norm_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance) );
cell_times            = cell_times(  ~cellfun(@isempty,cell_times) );

[coeffs,labels]=extract_features( cell_times, norm_cell_reflectance, [66 99], vid_type );

% [predictions, score] = randforest_classifier.predict(coeffs);

testTable = [array2table(coeffs) cell2table(labels)];
testTable.Properties.VariableNames = varNames;

[predictions, score] = randforest_classifier.predictFcn(testTable);

confusionmat(labels,predictions)


score(:,2) = sign(score(:,2)).*log( abs(score(:,2))+1);
score(:,2) = (score(:,2)-min(score(:,2)));
score(:,2) = (score(:,2)./max(score(:,2)))*100+1;

figure(2); imagesc(ref_image); colormap gray; hold on;

for i=1:length(coords_used)    
    plot( coords_used(i,1), coords_used(i,2),'*','Color',cmapped( ceil( score(i,2) ),:) );
end

% plot( coords_used(strcmp(predictions,'stimulus'),1), coords_used(strcmp(predictions,'stimulus'),2),'r*');
% plot( coords_used(strcmp(predictions,'control'),1), coords_used(strcmp(predictions,'control'),2),'b*'); hold off;
title('Unstimulated cones');

