%% TypeAOSLOCones
%
% Use MDS to try to type cones based on microspectrophometry data provided
% by Ram Sebasen.  
%
% The basic data values provided are intensity reflectance values and not
% absorptance. So if (post - pre) is a large number, this would mean that
% cone underwent a large change in its reflectivity when it was bleached
% and and then imaged with the 543 nm light. Large post/pre value would correspond
% to bright pixels and vice-versa.
%
% Here is the data format for file absorptance_cones_per_bleaching_cycle.mat, as described by Ram.
% These are the data for the individual runs.
%     It should be pretty self explanatory. When you load this into the workspace, you
%     will see 3 structures, one for each bleach condition : L-bleach, M-bleach,
%     No-bleach. 
% 
%     Each structure is 11 elements long, corresponding to the bleaching cycles.
%     For every bleach cycle, there is absorptance, pre-bleach intensity
%     (Intensity immediately after a selective L/M/No bleach) and post-intensity
%     (Intensity after a full bleach with 550nm imaging beam, following the selective bleaches).
% 
%     Just a word about how I got this. 
%     Cones were identified across all bleaching cycles, their intensities were
%     fit with an exponential, and the intensity difference at the start
%     and end were computed. These are the numbers in the mat file.
%
%     The absorptance field is the absolute value of the difference between
%     post and pre, which we don't actually need given that we have the pre
%     and post.  And Ram subsequently indicated that taking the absolute
%     value doesn't really make sense.  So I don't examine that field here.
%
% We also have mean data, in file absorbtance_cones_mean.mat.  This has three
% variables, absorptance_L_bleach, absorbtance_M_bleach, and absorbtance_No_bleach.
% Ram writes:
%     Please find attached a mat file with the mean absorptance for each cone.
%     Loading the file onto the workspace will put 3 variables into the
%     workspace, each with 1x924 dimensions. The names of the variables are
%     self explanatory. These are calculated by clumping together the
%     intensities from all the bleaching runs and fitting them together with
%     one curve. Similar to the individual bleaching runs, each absorptance
%     value is the difference of post and pre intensity calculated from the
%     curve . I have observed that this tends to be more robust.
% 
%     For parsing out S-cones from this mean data, try plotting the histogram
%     of the 3-d euclidean distance from these absorptance values. This would
%     look like :  sqrt( absorptance_L_bleach(ii).^2 +
%     absorptance_M_bleach(ii).^2 + absorptance_No_bleach(ii).^2)
%
% Dependencies:
%     This depends on having the Matlab version of libsvm
%     (http://www.csie.ntu.edu.tw/~cjlin/libsvm/) on your path, for the svm
%     classification.  I have found this easier to work with tha Matlab's
%     version of svm.  This code is included as part of theisetbio toolbox available
%     on gitHub.
%       https://github.com/isetbio/isetbio
% 
%     There are also likely some dependencies on Brainard Lab standard routines.  The
%     BrainardLabToolbox is available from the DavidBrainard repositories on gitHub.
%       https://github.com/DavidBrainard/BrainardLabToolbox
%     It is also possible that there are calls to routines in the Psychophysics Toolbox.
%     You can also get that from gitHub.
%       https://github.com/Psychtoolbox-3/Psychtoolbox-3
%     Or you can go all in and use our installer to produce our standard Brainard Lab
%     setup, which will get you both the BrainardLab Toolbox and the Psychophysics Toolbox.
%       https://github.com/DavidBrainard/BrainardLabInstaller
%     (But you'll still need to add isetbio by hand.)
%
%     The one BrainardLab dependency I am sure of is FigureSave().  Calls to this could just be
%     commented out and then the figures wouldn't save to disk but nothing else bad would happen.
%     Could replace calls to FigureSave with calls to Matlab's saveas(), but the arg
%     order is different for saveas (reverse first two arguements).
%
%     Some plotting things break under pre-release 2014b, presumably because the
%     plotting stuff in Matlab was extensively redone.  I think the time to worry
%     about this is with 2015a, when Matlab has had a chance to beat out the bugs
%     that are surely in 2014b.
% 
% Ideas to try:
%   1) Explicitly fit mixtures of Gaussians to univariate diagnostic
%   variables rather than just set a threshold by eye.  This is what Ram
%   currently does and it would be nice to have it working here.  The stats
%   toolbox has an implementation of em for this.  Could also try mixture
%   of multivariate Gaussians in higher dimensions.
%
%   2) Try kmeans clustering and perhaps other unsupervised clustering
%   methods, both on the full multivariate data and on the MDS solution.
%
%   3) Try svm classification based on high-dimensional data spaces and see
%   if this leads to anything good that is too hard to visualize otherwise.
%
%   4) Try some sort of combination of clustering based on separate
%   classifications from data and MDS solutions.
%
%   5) Try to look systematically at whether data varies slowly over space, which
%   I think is the thing that would need to be true for the correlation method to
%   work well.
%
%   6) For the data set we're working with here, the S cones seem as hard to get, or 
%   harder, than the L and M cones.  Ram says:
%     This is what I would expect from the data. The bimodal S-cone plot I showed you
%     in Philadelphia in July is from a little further eccentricity (1.5 -2 deg), so hopefully
%     better spatial resolution. Given that there are no bleach conditions for
%     this particular plot, the overlap of the modes can mostly be explained by optical
%     quality with AO. In the L-bleach and  M-bleach cases, the overlap of the histogram also
%     has to do with insufficient/inoptimal bleach. In this subject, (Heidi's
%     husband) we had to stick to a closer eccentricity (1 deg) because we were comparing
%     from old data. I was able to fit a sum of Gaussians to this histogram and parse out S-cones,
%     which matched well with Heidi's old AO fundus camera densitometry.
%   It is probably worth spending a little time optimizing the S cone segregation and then 
%   seeing how the L/M cones segregation works when we first parse out the S cones well.
% 
% 7/22/14  dhb  Wrote it.
% 8/29/14  dhb  Added new mean data file from Ram, cleaning it up.
% 8/30/14  dhb  Played around with Gaussian mixtures.  Some promise, more work needed.
%               Also played wround with k-means classification.  That didn't work so well.
% 9/1/14   dhb  Fit S-cone diagnostic variable with mixture of 2 Gaussians.  

%% Clear
clear; close all;

%% Some basic parameters
figureDir = 'xFigures';
figureType = 'png';
scatterMarkerSize = 4;
doMeanCheckFigure = false;
doCorrelationDistanceFigure = false;
doMDSQualityCheckFigure = false;

%% Setup
if (~exist(figureDir,'dir'))
    mkdir(figureDir);
end

%% Load the individual run data and put it into standard variable names.
fileName = 'absorptance_cones_per_bleaching_cycle';
loadedData = load(fileName);
nCycles = length(loadedData.L_bleach_cycles_absorptance_pre_post);
for i = 1:nCycles
    LConeBleachPre(:,i) = loadedData.L_bleach_cycles_absorptance_pre_post(i).pre';
    LConeBleachPost(:,i) = loadedData.L_bleach_cycles_absorptance_pre_post(i).post';
    
    MConeBleachPre(:,i) = loadedData.M_bleach_cycles_absorptance_pre_post(i).pre';
    MConeBleachPost(:,i) = loadedData.M_bleach_cycles_absorptance_pre_post(i).post';
    
    NoBleachPre(:,i) = loadedData.No_bleach_cycles_absorptance_pre_post(i).pre';
    NoBleachPost(:,i) = loadedData.No_bleach_cycles_absorptance_pre_post(i).post';
end
LConeBleachDiff = LConeBleachPost-LConeBleachPre;
MConeBleachDiff = MConeBleachPost-MConeBleachPre;
NoBleachDiff = NoBleachPost-NoBleachPre; 
clear loadedData

%% Get the mean data from the individual run data, for each cone/bleach condition
LConeBleachMeanPre = mean(LConeBleachPre,2);
LConeBleachMeanPost = mean(LConeBleachPost,2);
LConeBleachMeanDiff = mean(LConeBleachDiff,2);
MConeBleachMeanPre = mean(MConeBleachPre,2);
MConeBleachMeanPost = mean(MConeBleachPost,2);
MConeBleachMeanDiff = mean(MConeBleachDiff,2);
NoBleachMeanPre = mean(NoBleachPre,2);
NoBleachMeanPost = mean(NoBleachPost,2);
NoBleachMeanDiff = mean(NoBleachDiff,2);

%% Load the mean data as we got it from Ram.
%
% This takes the average of the data and fits it, rather
% than being the average of the fits to the individual
% runs.  So it won't be exactly the same as what we compute
% just above, but it would be surprising if it were vastly
% different. We will use variable name segment "RamMean" for
% these versions.
fileName = 'absorptance_cones_mean';
loadedData = load(fileName);
LConeBleachRamMeanDiff = loadedData.absorptance_L_bleach';
MConeBleachRamMeanDiff = loadedData.absorptance_M_bleach';
NoBleachRamMeanDiff = loadedData.absorptance_No_bleach';
clear loadedData

% Plot to see if the two methods of averaging bear sensible relation
% to each other.  I would say yes they are reasonably similar, so
% this minute part of the multiverse is relatively sane.
if (doMeanCheckFigure)
    limLow = -0.2; limHigh = 0.8;
    meanCheckFig = figure; clf;
    set(gcf,'Position',[100,200,1500,500]);
    subplot(1,3,1); hold on
    set(gca,'FontName','Hevetica','FontSize',14);
    plot(LConeBleachMeanDiff,LConeBleachRamMeanDiff,'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
    plot([limLow limHigh],[limLow limHigh],'r','LineWidth',2);
    xlabel('Mean of individual runs fit separately','FontSize',18);
    ylabel('Mean of individual run data fit once','FontSize',18);
    title('L cone bleach difference ','FontSize',18);
    axis('square');
    axis([limLow limHigh limLow limHigh]);
    subplot(1,3,2); hold on
    set(gca,'FontName','Hevetica','FontSize',14);
    plot(MConeBleachMeanDiff,MConeBleachRamMeanDiff,'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
    plot([limLow limHigh],[limLow limHigh],'g','LineWidth',2);
    xlabel('Mean of individual runs fit separately','FontSize',18);
    ylabel('Mean of individual run data fit once','FontSize',18);
    title('M cone bleach difference','FontSize',18);
    axis('square');
    axis([limLow limHigh limLow limHigh]);
    subplot(1,3,3); hold on
    set(gca,'FontName','Hevetica','FontSize',14);
    plot(NoBleachMeanDiff,NoBleachRamMeanDiff,'ko','MarkerFaceColor','k','MarkerSize',scatterMarkerSize);
    plot([limLow limHigh],[limLow limHigh],'k','LineWidth',2);
    xlabel('Mean of individual runs fit separately','FontSize',18);
    ylabel('Mean of individual run data fit once','FontSize',18);
    title('No bleach difference','FontSize',18);
    axis('square');
    axis([limLow limHigh limLow limHigh]);
    FigureSave(fullfile(figureDir,'meanCheckFig'),meanCheckFig,figureType);
end

%% Try to identify S cones in a quick and drty univariate classification.
%
% The separation I can get so far is not spectacular.
%
% Select which type of mean data to use
%       'meanofindfits' - Use mean of fits to individual runs
%       'fittomeandata' - Use fit to mean of individual run data
whichSDiagnoseMeanDataType = 'meanofindfits';
switch (whichSDiagnoseMeanDataType)
    case 'meanofindfits'
        SDiagnose3D = [LConeBleachMeanDiff MConeBleachMeanDiff NoBleachMeanDiff];
        SDiagnoseDataTypeStr = 'mean of ind fits';
    case 'fittomeandata'
        SDiagnose3D = [LConeBleachRamMeanDiff MConeBleachRamMeanDiff NoBleachRamMeanDiff];
        SDiagnoseDataTypeStr = 'fit to mean data';
end

% Set a univariate variable to threshold on to find S cones, and find them
%       'allbleachdiffdistance' - Euclidean distance of three bleach differences
%       'nobleachpre' - Mean of pre measurements in no bleach (dark adapted) condition
%       'nobleachdiff' - Mean difference from no bleach condition
whichSDiagnoseVariable = 'allbleachdiffdistance';
switch (whichSDiagnoseVariable)
    % This is the method Ram recommended.  The logic is that S cones
    % would change less over the course of an imaging run, because they
    % hardly absorb light at all at the imaging wavelengths.
    case 'allbleachdiffdistance'
        SDiagnoseVariable = sqrt(diag(SDiagnose3D*SDiagnose3D'));
        SDiagnoseThresh = 0.31;
        SDiagnoseSign = -1;
        SDiagnoseLow = 0;
        SDiagnoseHigh = 1;
        xLabelStr = sprintf('Euclidean Distance of All Bleach Diffs; %s',SDiagnoseDataTypeStr);
    case 'allprebleachdistance'
        SDDiagnoseRaw = [NoBleachMeanPre(:) LConeBleachMeanPre(:) MConeBleachMeanPre(:)];
        SDiagnoseVariable = sqrt(diag(SDDiagnoseRaw*SDDiagnoseRaw'));
    
    % This would be the same logic as above, but would just look at the no bleach (pre-dark adapted)
    % condition, on the grounds that t
    case 'nobleachdiff'
        SDiagnoseVariable = SDiagnose3D(:,3);
        SDiagnoseThresh = 0.27;
        SDiagnoseSign = -1;
        SDiagnoseLow = 0;
        SDiagnoseHigh = 1;
        xLabelStr = sprintf('No Bleach Diffs; %s',SDiagnoseDataTypeStr);
        
        % The logic here is that S cones would have high reflectance at the imaging
        % wavelength, and right after dark adaptation L and M cnoes would have low
        % reflectance.
    case 'nobleachpre'
        SDiagnoseVariable = NoBleachPre;
        SDiagnoseThresh = 0.52;
        SSign = 1;
        SDiagnoseLow = -1.5;
        SDiagnoseHigh = 1.5;
        xLabelStr = sprintf('No Bleach Pre');     
end

% Classify S cones by fitting a mixture of two Gaussians.
%
% Determine which index goes with which class by looking at means
% of the fit.
gaussMixtureOptions = statset('MaxIter',1000);
gaussMixtureOnSDiagnoseObj = gmdistribution.fit(SDiagnoseVariable,2,'options',gaussMixtureOptions );
gaussMixtureOnSDiagnoseIDX = gaussMixtureOnSDiagnoseObj.cluster(SDiagnoseVariable);
if (gaussMixtureOnSDiagnoseObj.mu(1) > gaussMixtureOnSDiagnoseObj.mu(2))
    SIDXVal = 2;
    LMIDXVal = 1;
else
    SIDXVal = 1;
    LMIDXVal = 2;
end
SFastIndex = find(gaussMixtureOnSDiagnoseIDX == SIDXVal);
LMLogical = (gaussMixtureOnSDiagnoseIDX == LMIDXVal);
LMFastIndex = find(LMLogical);

% Report on S-cone findings
fprintf('%%S cones identified: %0.0f\n',round(100*length(SFastIndex)/length(SDiagnoseVariable)));
titleStr = sprintf('Threshold %0.2f; %%S cones identified: %d\n',SDiagnoseThresh,round(100*length(SFastIndex)/length(SDiagnoseVariable)));

% Histogram of classified S diagnostc variable with Guassian
% mixture fit.
%
% Note that the patches seem to come back in the opposite
% order from that in which the calls to hist were made, so
% that h(1) corresponds to the last call to hist().
SDiagnoseFig = figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
[~,x] = hist(SDiagnoseVariable,50);
nS = hist(SDiagnoseVariable(SFastIndex),x)';
nLM = hist(SDiagnoseVariable(LMFastIndex),x)';
bar(x,[nS nLM],'stacked');
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','y'); set(h(1),'EdgeColor','k');
set(h(2),'FaceColor','b'); set(h(2),'EdgeColor','k');
gaussMixtureOnSDiagnosePred = gaussMixtureOnSDiagnoseObj.pdf(x')*(x(2)-x(1))*length(SDiagnoseVariable);
plot(x,gaussMixtureOnSDiagnosePred,'k','LineWidth',2);
xlim([SDiagnoseLow SDiagnoseHigh]);
xlabel(xLabelStr,'FontSize',18);
ylabel('Count','FontSize',18);
title(titleStr,'FontSize',18);
FigureSave(fullfile(figureDir,['SConeDiagnose_' whichSDiagnoseVariable '_' whichSDiagnoseMeanDataType '_Fig']),SDiagnoseFig ,figureType);

%% Classify L versus M using threshold set by eye, based on theta histogram
%
% The purpose here is to get some idea of what the right answer
% is, so that we're not completely exploring in the dark.  This
% uses a threshold on the 'theta' variable that I believe the Roorda
% lab uses.
%
% Here we have not implemented different choices of averaging and are just
% using the mean of the fits to the individual runs.
%
% The thetaThresh variable sets the boundary between L and M in this space.
% the thetaBoundary variable sets a 'neutral zone' around the boundary that
% we ignore when setting up the trainging set for SVM classification.
%
% Not 100% sure I have the logic of the sign correct for which are L and
% which are M.
thetaThresh = 0.76;
thetaBoundary = 0.05;
LMTheta = atan2(MConeBleachMeanDiff,LConeBleachMeanDiff);
LFastIndex = find(LMTheta > thetaThresh & LMLogical);
MFastIndex = find(LMTheta <= thetaThresh & LMLogical);
LFastBoundaryIndex = find(LMTheta > (thetaThresh+thetaBoundary) & LMLogical);
MFastBoundaryIndex = find(LMTheta <= (thetaThresh-thetaBoundary) & LMLogical);
fprintf('LM ratio identified fast: %0.2f\n',length(LFastIndex)/length(MFastIndex));
fprintf('LM training set ratio identified fast: %0.2f\n',length(LFastBoundaryIndex)/length(MFastBoundaryIndex));
fprintf('Training set is %d of %d LM cones\n',length(LFastBoundaryIndex)+length(MFastBoundaryIndex),length(LFastIndex)+length(MFastIndex));

% Histogram of theta
%
% Note that the patches seem to come back in the opposite
% order from that in which the calls to hist were made, so
% that h(1) corresponds to the last call to hist().
LMThetaFig = figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
[~,x] = hist(LMTheta,50);
hist(LMTheta(LFastIndex),x);
hist(LMTheta(MFastIndex),x);
plot([thetaThresh,thetaThresh],[0 10],'k','LineWidth',3);
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','g'); set(h(1),'EdgeColor','k');
set(h(2),'FaceColor','r'); set(h(2),'EdgeColor','k');
xlim([0 pi]);
xlabel('LM Bleach Diff Theta','FontSize',18);
ylabel('Count','FontSize',18);
title('Histogram of LM theta','FontSize',18);
FigureSave(fullfile(figureDir,['LMThetaHistogramFig']),LMThetaFig,figureType);

% Also plot the L/M bleach differences in 2D
LMBleachDiffFig = figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot(LConeBleachMeanDiff(LFastIndex),MConeBleachMeanDiff(LFastIndex),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(MFastIndex),MConeBleachMeanDiff(MFastIndex),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(SFastIndex),MConeBleachMeanDiff(SFastIndex),'bo','MarkerFaceColor','b','MarkerSize',scatterMarkerSize);
axis('square'); xlim([-0.1 1]); ylim([-0.1 1.0]);
xlabel('L Bleach Difference','FontSize',18); 
ylabel('M Bleach Difference','FontSize',18);
title('Fast L/M Classification Based on Theta','FontSize',18);
FigureSave(fullfile(figureDir,['LMBleachDiff2DFig']),LMThetaFig,figureType);

%% Look at the difference (post-pre) data in 3D 
% And try to classify in various ways.

% Get the mean data
prePostMean3DData = [LConeBleachMeanDiff MConeBleachMeanDiff NoBleachMeanDiff];

% Make a plot
az = 2; el = -22;
figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot3(prePostMean3DData(LFastIndex,1),prePostMean3DData(LFastIndex,2),prePostMean3DData(LFastIndex,3),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot3(prePostMean3DData(MFastIndex,1),prePostMean3DData(MFastIndex,2),prePostMean3DData(MFastIndex,3),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot3(prePostMean3DData(SFastIndex,1),prePostMean3DData(SFastIndex,2),prePostMean3DData(SFastIndex,3),'bo','MarkerFaceColor','b','MarkerSize',scatterMarkerSize);
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
zlabel('Z','FontSize',18);
title('Pre-Post Mean Data, Fast Classifications','FontSize',18);
view([az el]);

% Try k-means on mean data.
%
% This produces a sensible result at a gross level, but is not so hot
% in detail.  I'm not sure k-means does so well when the distribution
% within class is not circular in the data space.  
%
% Which returned index corresponds to which cone class is arbitrary,
% and these were just figured out here by looking.
[kMeansOnMeanDataIDX,kmeansOnDataCentroid,~,kmeansOnDataDistances] = kmeans(prePostMean3DData,3);
LKMeansOnMeanDataIndex = find(kMeansOnMeanDataIDX == 2);
MKMeansOnMeanDataIndex = find(kMeansOnMeanDataIDX == 3);
SKMeansOnMeanDataIndex = find(kMeansOnMeanDataIDX == 1);

kMeansOnMeanData2DFig = figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot(LConeBleachMeanDiff(LKMeansOnMeanDataIndex ),MConeBleachMeanDiff(LKMeansOnMeanDataIndex ),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(MKMeansOnMeanDataIndex ),MConeBleachMeanDiff(MKMeansOnMeanDataIndex ),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(SKMeansOnMeanDataIndex),MConeBleachMeanDiff(SKMeansOnMeanDataIndex),'bo','MarkerFaceColor','b','MarkerSize',scatterMarkerSize);
axis('square'); xlim([-0.1 1]); ylim([-0.1 1.0]);
xlabel('L Bleach Difference','FontSize',18); 
ylabel('M Bleach Difference','FontSize',18);
title('K-Means on 3D Mean Data Classification','FontSize',18);
FigureSave(fullfile(figureDir,['KMeansOnMeanData2DFig']),kMeansOnMeanData2DFig,figureType);

% Try mixture of Gaussians on 3D Mean data
%
% This works better than k-means, to my eye.  Since we don't actually
% know ground truth, it is a bit hard to say but I think there are too
% many S-cones identified.  This might be improved by filtering those
% out before trying to find the L and M cones this way.
%
% Which returned index corresponds to which cone class is arbitrary,
% and these were just figured out here by looking.  Unfortunately,
% it is stochastic and comes out differently in each run.  Can fix
% by looking at relation between mixture means, I suppose, but didn't
% do that yet.
gaussMixtureOptions = statset('MaxIter',1000);
gaussMixtureOnMeanDataObj = gmdistribution.fit(prePostMean3DData,3,'options',gaussMixtureOptions );
gaussMixtureOnMeanDataIDX = gaussMixtureOnMeanDataObj.cluster(prePostMean3DData);
LGaussMixtureOnMeanDataIndex = find(gaussMixtureOnMeanDataIDX == 2);
MGaussMixtureOnMeanDataIndex = find(gaussMixtureOnMeanDataIDX == 1);
SGaussMixtureOnMeanDataIndex = find(gaussMixtureOnMeanDataIDX == 3);

gaussMixtureOnMeanData2DFig = figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot(LConeBleachMeanDiff(LGaussMixtureOnMeanDataIndex ),MConeBleachMeanDiff(LGaussMixtureOnMeanDataIndex ),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(MGaussMixtureOnMeanDataIndex ),MConeBleachMeanDiff(MGaussMixtureOnMeanDataIndex ),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(SGaussMixtureOnMeanDataIndex),MConeBleachMeanDiff(SGaussMixtureOnMeanDataIndex),'bo','MarkerFaceColor','b','MarkerSize',scatterMarkerSize);
axis('square'); xlim([-0.1 1]); ylim([-0.1 1.0]);
xlabel('L Bleach Difference','FontSize',18); 
ylabel('M Bleach Difference','FontSize',18);
title('Gaussian Mixture on 3D Mean Data Classification','FontSize',18);
FigureSave(fullfile(figureDir,['GaussMixtureOnMeanData2DFig']),gaussMixtureOnMeanData2DFig,figureType);

%% Build the correlation matrix. 
%
% There are a number of ways to do this.

% How to aggregate over trials
%   'useMean'
%   'useRaw'
correlationMeanMethod = 'useRaw';
switch (correlationMeanMethod)
    case 'useMean'
        theLConeBleachPre = LConeBleachMeanPre;
        theLConeBleachPost = LConeBleachMeanPost;
        theLConeBleachDiff = LConeBleachMeanDiff;
        theMConeBleachPre = MConeBleachMeanPre;
        theMConeBleachPost = MConeBleachMeanPost;
        theMConeBleachDiff = MConeBleachMeanDiff;
        theNoBleachPre = NoBleachMeanPre;
        theNoBleachPost = NoBleachMeanPost;
        theNoBleachDiff = NoBleachMeanDiff;
    case 'useRaw'          
        theLConeBleachPre = LConeBleachPre;
        theLConeBleachPost = LConeBleachPost;
        theLConeBleachDiff = LConeBleachDiff;
        theMConeBleachPre = MConeBleachPre;
        theMConeBleachPost = MConeBleachPost;
        theMConeBleachDiff = MConeBleachDiff;
        theNoBleachPre = NoBleachPre;
        theNoBleachPost = NoBleachPost;
        theNoBleachDiff = NoBleachDiff;
    otherwise
        error('Unknown meanMethod specified');
end

% What to stick into the correlation matrix
%  'LandMBleachDiff'
%  'LandMBleachPrePost'
%  'LandMBleachPre' 
%  'LandMBleachPost'
%  'AllPrePost'
correlationWhatToInclude = 'LandMBleachPrePost';
switch (correlationWhatToInclude)
    case 'LandMBleachDiff'
        theCorrData = [theLConeBleachDiff theMConeBleachDiff];
    case 'LandMBleachPrePost'
        theCorrData = [theLConeBleachPre theLConeBleachPost theMConeBleachPre theMConeBleachPost];
    case 'LandMBleachPre'
        theCorrData = [theLConeBleachPre theMConeBleachPre];
     case 'LandMBleachPost'
        theCorrData = [theLConeBleachPost theMConeBleachPost];
    case 'AllPrePost'
        theCorrData = [theLConeBleachPre theLConeBleachPost theMConeBleachPre theMConeBleachPost theNoBleachPre theNoBleachPost];
    otherwise
        error('Unknown whatToInclude specified');
end
theCorr = corr(theCorrData');

%% Build 'distance' matrix
%
% Have to condition the raw distance matrix because of numerical issues.
theDistance = -log10(theCorr+2);
theDistance = theDistance-min(theDistance(:));
for i = 1:size(theDistance,1)
    if (theDistance(i,i) > 1e-10)
        %error('Non-zero diag entry of raw distance matrix');
    end
    theDistance(i,i) = 0;
end
for i = 1:size(theDistance,1)
    for j = 1:i
        if (abs((theDistance(i,j) - theDistance(j,i)) > 1e-10))
            error('Asymmetric raw distance matrix');
        end
        theDistance(i,j) = theDistance(j,i);
    end
end

% Check that the correlation versus distance plot looks sensible.
if (doCorrelationDistanceFigure)
    figure; clf; hold on
    set(gca,'FontName','Hevetica','FontSize',14);
    plot(theCorr(:),theDistance(:),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
    xlabel('Correlation','FontSize',18);
    ylabel('Distance','FontSize',18);
    title('Check of Correlation->Distance Conversion','FontSize',18);
end

%% Get the MDS solution
fullSolution = mdscale(theDistance,3);

% Plot quality of solution' account of the distances it was given.
% This should generally be an upward heading scattergram.  There
% will be considerable scatter because the individual runs are pretty
% noisy.
if (doMDSQualityCheckFigure)
    predDistance = squareform(pdist(fullSolution));
    figure; clf; hold on
    set(gca,'FontName','Hevetica','FontSize',14);
    plot(theDistance(:),predDistance(:),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
    xlabel('Correlation-Based Distance','FontSize',18);
    ylabel('MDS Predicted Distance','FontSize',18);
    title('Quality of MDS Solution','FontSize',18);
    axis('square');
end

% Plot solution with LMS from fast classifiation in red, green, and blue.
%
% Whether it is useful to show purported S cones here is a matter to think about.
az = 2; el = -22;
figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot3(fullSolution(LFastIndex,1),fullSolution(LFastIndex,2),fullSolution(LFastIndex,3),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot3(fullSolution(MFastIndex,1),fullSolution(MFastIndex,2),fullSolution(MFastIndex,3),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
%plot3(fullSolution(SFastIndex,1),fullSolution(SFastIndex,2),fullSolution(SFastIndex,3),'bo','MarkerFaceColor','b','MarkerSize',scatterMarkerSize);
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
zlabel('Z','FontSize',18);
title('MDS Solution 3D Plot','FontSize',18);
view(az,el);

%% Get the training set from the solution
LLabel = 1;
MLabel = -1;
LTrainingSet = fullSolution(LFastBoundaryIndex,:);
MTrainingSet = fullSolution(MFastBoundaryIndex,:);
LLabels = LLabel*ones(size(LFastBoundaryIndex));
MLabels = MLabel*ones(size(MFastBoundaryIndex));

% Plot training set in 3D
az = 5; el = 16;
figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot3(LTrainingSet(:,1),LTrainingSet(:,2),LTrainingSet(:,3),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot3(MTrainingSet(:,1),MTrainingSet(:,2),MTrainingSet(:,3),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
zlabel('Z','FontSize',18)
title('MDS Solution 3D Plot (Training Set)');
view(az,el);

% Histogram of first MDS solution dimension.
%
% This is a rough and ready plot.
%   1) here is no reason why the first dimension of the solution needs to be the best
%   one to separate the L and M cones.
%   2) The way I'm doing the histogram the second set of data plotted overwrites
%   the first.  Need to figure out how to do a stacked histogram, which probably involves
%   the bar() function.
% Those points noted, the separation on the first dimension is pretty good.
figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
[~,x] = hist([LTrainingSet(:,1) ; MTrainingSet(:,1)],50);
hist(LTrainingSet(:,1),x);
hist(MTrainingSet(:,1),x);
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','g'); set(h(1),'EdgeColor','k');
set(h(2),'FaceColor','r'); set(h(2),'EdgeColor','k');
xlabel('MDS Dimension 1','FontSize',18);
ylabel('Count','FontSize',18);
title('MDS Dimension 1 (Training Set)','FontSize',18);

%% Use theta-based rough classifcation with boundary to try build a classifier
% on the solution space.
%
% Get path to svmlib.  Ugh.  This is to call it 
% which full path and avoid convlicts with 
% matlab functions of same name.
libsvmPath = fileparts(which('svmpredict'));
curDir = pwd;

% Train classifier and make predictions
LIBSVM_QUIET = false;
trainOpts = '-s 1 -t 0';
if (LIBSVM_QUIET)
    trainOpts = [trainOpts ' -q'];
end
predictOpts = '';
if (LIBSVM_QUIET)
    predictOpts = [predictOpts ' -q'];   
end
cd(libsvmPath);
classifyInfo = svmtrain([LLabels ; MLabels], [LTrainingSet ; MTrainingSet], trainOpts);
testClassify = svmpredict([LLabels ; MLabels], [LTrainingSet ; MTrainingSet],classifyInfo, predictOpts);
cd(curDir);

% Predict using the classifier
cd(libsvmPath);
[predict,~,svmDecisionVals] = svmpredict(rand(size(fullSolution,1),1), fullSolution, classifyInfo, predictOpts);
cd(curDir);
LSVMIndex = find(predict == LLabel & LMLogical);
MSVMIndex = find(predict == MLabel & LMLogical);
fastlabels = ones(size(predict));
fastlabels(LFastIndex) = LLabel;
fastlabels(MFastIndex) = MLabel;
disagreeIndex = find(predict ~= fastlabels & LMLogical);
LDisagreeIndex = find(predict == LLabel & fastlabels == MLabel & LMLogical);
MDisagreeIndex = find(predict == MLabel & fastlabels == LLabel & LMLogical);
fprintf('Fast and SVM based answers differ for %d of %d LM cones\n',length(disagreeIndex),length(LSVMIndex)+length(MSVMIndex));

% Histogram of svm decision values on training set
figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
[~,x] = hist([svmDecisionVals(LFastBoundaryIndex) ; svmDecisionVals(MFastBoundaryIndex)],50);
nL = hist(svmDecisionVals(LFastBoundaryIndex),x)';
nM = hist(svmDecisionVals(MFastBoundaryIndex),x)';
bar(x,[nL nM],'stacked');
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','g'); set(h(1),'EdgeColor','k');
set(h(2),'FaceColor','r'); set(h(2),'EdgeColor','k');
xlabel('SVM Decision Value','FontSize',18);
ylabel('Count','FontSize',18);
title('SVM Decision Value (Training Set)','FontSize',18);

% Histogram of svm decision values
svmDecisionValFig = figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
[~,x] = hist([svmDecisionVals(LSVMIndex) ; svmDecisionVals(MSVMIndex)],50);
nL = hist(svmDecisionVals(LSVMIndex),x)';
nM = hist(svmDecisionVals(MSVMIndex),x)';
bar(x,[nL nM],'stacked');
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','g'); set(h(1),'EdgeColor','k');
set(h(2),'FaceColor','r'); set(h(2),'EdgeColor','k');
xlabel('SVM Decision Value','FontSize',18);
ylabel('Count','FontSize',18);
title('SVM Decision Values','FontSize',18);

% Plot in 2D
figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot(LConeBleachMeanDiff(LSVMIndex),MConeBleachMeanDiff(LSVMIndex),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(MSVMIndex),MConeBleachMeanDiff(MSVMIndex),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(LDisagreeIndex),MConeBleachMeanDiff(LDisagreeIndex),'ko','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(MDisagreeIndex),MConeBleachMeanDiff(MDisagreeIndex),'ko','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(SFastIndex),MConeBleachMeanDiff(SFastIndex),'bo','MarkerFaceColor','b','MarkerSize',scatterMarkerSize);
axis('square');
xlim([-0.1 1]); ylim([-0.1 1.0]);
xlabel('L bleach difference','FontSize',18);
ylabel('M bleach difference','FontSize',18);
title('L/M Classification Based on SVM of MDS');

% SVM classification 3D plot
az = 5; el = 16;
figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot3(fullSolution(LSVMIndex,1),fullSolution(LSVMIndex,2),fullSolution(LSVMIndex,3),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot3(fullSolution(MSVMIndex,1),fullSolution(MSVMIndex,2),fullSolution(MSVMIndex,3),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot3(fullSolution(LDisagreeIndex,1),fullSolution(LDisagreeIndex,2),fullSolution(LDisagreeIndex,3),'ko','MarkerFaceColor','r','MarkerSize',scatterMarkerSize+2);
plot3(fullSolution(MDisagreeIndex,1),fullSolution(MDisagreeIndex,2),fullSolution(MDisagreeIndex,3),'ko','MarkerFaceColor','g','MarkerSize',scatterMarkerSize+2);
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
zlabel('Z','FontSize',18);
title('SVM 3D Plot','FontSize',18);
view(az,el);

%% Try unsupervised classification on the solution matrix
%
% This works better than k-means, to my eye.  Since we don't actually
% know ground truth, it is a bit hard to say but I think there are too
% many S-cones identified.  This might be improved by filtering those
% out before trying to find the L and M cones this way.
%
% Which returned index corresponds to which cone class is arbitrary,
% and these were just figured out here by looking.  Unfortunately,
% it is stochastic and comes out differently in each run.  Can fix
% by looking at relation between mixture means, I suppose, but didn't
% do that yet.
gaussMixtureOptions = statset('MaxIter',1000);
gaussMixtureOnMDSDataObj = gmdistribution.fit(fullSolution,3,'options',gaussMixtureOptions );
gaussMixtureOnMDSDataIDX = gaussMixtureOnMDSDataObj.cluster(fullSolution);
LGaussMixtureOnMDSDataIndex = find(gaussMixtureOnMDSDataIDX == 2);
MGaussMixtureOnMDSDataIndex = find(gaussMixtureOnMDSDataIDX == 1);
SGaussMixtureOnMDSDataIndex = find(gaussMixtureOnMDSDataIDX == 3);

gaussMixtureOnMDSData2DFig = figure; clf; hold on
set(gca,'FontName','Hevetica','FontSize',14);
plot(LConeBleachMeanDiff(LGaussMixtureOnMDSDataIndex ),MConeBleachMeanDiff(LGaussMixtureOnMDSDataIndex ),'ro','MarkerFaceColor','r','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(MGaussMixtureOnMDSDataIndex ),MConeBleachMeanDiff(MGaussMixtureOnMDSDataIndex ),'go','MarkerFaceColor','g','MarkerSize',scatterMarkerSize);
plot(LConeBleachMeanDiff(SGaussMixtureOnMDSDataIndex),MConeBleachMeanDiff(SGaussMixtureOnMDSDataIndex),'bo','MarkerFaceColor','b','MarkerSize',scatterMarkerSize);
axis('square'); xlim([-0.1 1]); ylim([-0.1 1.0]);
xlabel('L Bleach Difference','FontSize',18); 
ylabel('M Bleach Difference','FontSize',18);
title('Gaussian Mixture on MDS Classification','FontSize',18);
FigureSave(fullfile(figureDir,['GaussMixtureOnMDSData2DFig']),gaussMixtureOnMDSData2DFig,figureType);


        