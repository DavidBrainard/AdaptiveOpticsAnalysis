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
clear;

%% Dependent variable, peak reflectance response
sessionOneReflectance = [4.395
2.501
3.062
2.183
2.785
4.412
2.862
2.913
2.508
2.965
4.733
2.907
2.885
1.924
2.530
4.659
2.623
2.983
2.332
2.407];
sessionTwoReflectance = [4.503
3.358
2.939
2.599
3.367
4.525
3.275
3.059
2.847
3.302
4.291
3.533
2.868
2.838
3.212
4.167
3.377
2.966
2.675
3.332];

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
11115
11102
11057
11108
11112
11115
11102
11057
11108
11112
11115
11102
11057
11108
11112
11115];
subject = [subjectOrder ; subjectOrder];

stimulusOrder = {'153nW'
'153nW'
'153nW'
'153nW'
'153nW'
'306nW'
'306nW'
'306nW'
'306nW'
'306nW'
'917nW'
'917nW'
'917nW'
'917nW'
'917nW'
'2.75uW'
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
random = [1];

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

% Specify a 3-way model but without the three way interaction
% 
% It looks like that in this case, subject variability is treated as measurement noise,
% and the significance values of the various effects change. At least I think that
% is why they change.  I am not completely sure that this is the correct interpretation.
model = [1 0 0 ; 0 1 0; 0 0 1; 0 1 1; 1 1 0; 1 0 1];
[pTwoWay, tabTwoWay] = anovan(peakReflectance',{subject stimulus session}, 'model',model,'varnames', varNames, 'random', random);

%% Plot of data averaged over stimuli
theAvgData = [4.549592 4.371892
    2.723052 3.385613
    2.96045 2.957948
    2.236925 2.747275
    2.671575 3.30355];
theStdErr = 2*[0.086085 0.086191
    0.096844 0.053773
    0.039609 0.039425
    0.123565 0.065832
    0.125514 0.033115];

figure; clf; hold on; theLim = 5;
plot(theAvgData(:,1),theAvgData(:,2),'ro','MarkerSize',6,'MarkerFaceColor','k');
errorbarx = errorbar(theAvgData(:,1),theAvgData(:,2),theStdErr(:,1),'ko', 'horizontal', 'LineWidth', 0.5);
errorbary = errorbar(theAvgData(:,1),theAvgData(:,2),theStdErr(:,2),'ko', 'vertical', 'LineWidth', 0.5);
plot([0 theLim],[0 theLim],'k--', 'LineWidth', 1)
set(gca,'FontName','Times New Roman','FontSize',20, 'FontWeight', 'bold');
xlabel('Session 1');
ylabel('Session 2');
axis('square');
xlim([1 theLim]); ylim([1 theLim]);

plot(theAvgData(:,1),theAvgData(:,2),'ro','MarkerSize',6,'MarkerFaceColor','r');
