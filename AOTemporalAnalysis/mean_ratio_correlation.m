
clear;
close all force;
clc;

rootDir = uigetdir(pwd);

fPaths = read_folder_contents(rootDir,'mat');


wbh = waitbar(0,['Dataset 0 of ' num2str(length(fPaths)) '.']);
ref=[];
ratio=[];
for i=1:size(fPaths,1)
       
    waitbar(i/length(fPaths), wbh, ['Dataset (' num2str(i) ' of ' num2str(length(fPaths)) ').']);

    load(fullfile(rootDir,fPaths{i}));
    
    ref   = [ref; prestim_ref];
    ratio = [ratio; prestim_ratio'];
    
end

close(wbh);


mdl= fitlm(ratio,ref,'linear');

coeff = mdl.Coefficients.Estimate;

plot(ratio,ref,'b.',ratio, ratio*coeff(2)+coeff(1),'r');
xlabel('Mean ratio (stim/control)-1');
ylabel('Stimulus Reflectance Standard Deviation (Subtracted)');
title(['Prestimulus reflectance using CONTROL normalization vs mean ratio w/regression: (R^2 ' num2str(mdl.Rsquared.Adjusted) ')']);