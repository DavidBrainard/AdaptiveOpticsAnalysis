%%AOAnova  Anova for Ray's optoretinogram data
%
% Description:
%    Based on BL anovan tutorial.  This is what we said in the
%    pre-registration document.
%
%    After the above analysis, we will calculate descriptive statistics for
%    thepeak amplitude for each stimulus condition, (e.g.mean and median),
%    and compare and contrast the descriptive statistics between stimulus
%    conditions, sessions, and subjects. We will determine if reciprocity is
%    established by computing for each subject a two-wayANOVA (stimulus
%    condition x session) and determining if there is a significant main
%    effect on the average cone reflectance response of stimulus conditions,
%    of session, and whether there is any interaction. If necessary, a
%    post-hoc t-test will be used to determine where the differences in
%    cone reflectances lie. We will also determine if there is
%    intrasession variation (i.e.with the 3.53 μW stimulus) by computing a
%    one-way ANOVA to analyze significant differences between response
%    means within a session, using a random split of the within-session
%    3.53 uW trials.

% 12/21/20  dhb  Wrote it


%% Clear
clear; close all;

%% Dependent variable, peak reflectance response
sessionOneReflectance = [4.541
2.273
3.306
1.959
4.539
2.75
3.234
2.274
4.753
2.843
2.962
1.836
4.721
2.564
3.071
2.277];
sessionTwoReflectance = [4.341
3.249
3.095
2.484
4.476
3.362
3.271
2.827
4.228
3.587
2.917
2.893
4.12
3.349
3.144
2.595];

peakReflectance = [sessionOneReflectance ; sessionTwoReflectance];
    
%% Independent variables
%
% These are subject, illuminant, and target.  For
% our demonstration purposes, it doesn't matter what
% these mean, except that it is natural to consider
% subject a random variable (each subject is a draw
% from a population) and illuminant and target to
% be independent variables that are under experimental
% controlled (and thus fixed variables).
subjectOrder = [11102
11057
11108
11112
11102
11057
11108
11112
11102
11057
11108
11112
11102
11057
11108
11112];
subject = [subjectOrder ; subjectOrder];

stimulusOrder = {'153nW'
'153nW'
'153nW'
'153nW'
'306nW'
'306nW'
'306nW'
'306nW'
'917nW'
'917nW'
'917nW'
'917nW'
'2.75uW'
'2.75uW'
'2.75uW'
'2.75uW'
};
stimulus = {stimulusOrder{:} stimulusOrder{:}};

sessionOrder = ones(size(subjectOrder));
session = [sessionOrder ; 2*sessionOrder];

%% Specs for anovan
%
% random is a vector containing the indices of variables that should
% be considered random.  One's not listed are considered fixed.
%
% Set to [] to have all variables considered fixed.  Making the
% change does effect the significance of the main effect of some
% of the fixed variables, in some cases.  So we ought to understand
% clearly which we want when we design the anova.  Searching on
% 'fixed versus random effect anoval' on the web returns a number
% of reasonable clear descriptions.
random = [1 3];

% Names vfor the variables, to make the anova table more readable
varNames = strvcat('Subject', 'Stimulus', 'Session');

%% Run the full 3-way model with all interactions
%
% The full model is the default
[pFull, tabFull] = anovan(peakReflectance',{subject stimulus session}, 'model','full','varnames', varNames, 'random', random);

% Specify a 2-way model that focusses on stimulus and session.
% 
% It looks like that in this case, subject variability is treated as measurement noise,
% and the significance values of the various effects change. At least I think that
% is why they change.  I am not completely sure that this is the correct interpretation.
model = [0 1 0; 0 0 1; 0 1 1];
[pTwoWay, tabTwoWay] = anovan(peakReflectance',{subject stimulus session}, 'model',model,'varnames', varNames, 'random', random);

%% Plot of data averaged over stimuli
theAvgData = [4.6385 4.29125
    2.6075 3.38675
    3.14325 3.10675
    2.0865 2.69975];
theStdErr = 2*[0.057244359 0.076340656
    0.125681144 0.07136336
    0.077886001 0.073319364
    0.111972095 0.096202369];

figure; clf; hold on; theLim = 5;
plot(theAvgData(:,1),theAvgData(:,2),'ro','MarkerSize',12,'MarkerFaceColor','r');
errorbarX(theAvgData(:,1),theAvgData(:,2),theStdErr(:,1),'ro');
errorbarY(theAvgData(:,1),theAvgData(:,2),theStdErr(:,2),'ro');
plot([0 theLim],[0 theLim],'k:')
xlabel('Session 1');
ylabel('Session 2');
axis('square');
xlim([1 theLim]); ylim([1 theLim]);

