%%AOAnovaWithin  Anova for Ray's optoretinogram data
%
% Description:
%    Based on BL anovan tutorial.  This is what we said in the
%    pre-registration document.  This file does the within session ANOVAs
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
%
% Columns here are subject, trial, group, response
sessionOneData = [
11002	1	1	4.5551
11002	2	1	4.521
11002	3	1	4.8205
11002	4	1	4.9163
11002	5	1	4.8445
11057	1	1	2.1589
11057	2	1	2.6341
11057	3	1	2.542
11057	4	1	2.998
11057	5	1	2.812
11108	1	1	2.8097
11108	2	1	3.1785
11108	3	1	3.4561
11108	4	1	3.7456
11108	5	1	2.6009
11112	1	1	1.672
11112	2	1	2.2177
11112	3	1	2.4321
11112	4	1	2.7351
11112	5	1	2.1021
11002	1	2	3.984
11002	2	2	4.2623
11002	3	2	4.6293
11002	4	2	4.9283
11002	5	2	5.0428
11057	1	2	2.972
11057	2	2	2.322
11057	3	2	2.7309
11057	4	2	2.8187
11057	5	2	2.9215
11108	1	2	2.7021
11108	2	2	3.2
11108	3	2	3.2255
11108	4	2	2.935
11108	5	2	3.1174
11112	1	2	2.3041
11112	2	2	2.083
11112	3	2	2.3457
11112	4	2	2.6531
11112	5	2	2.9795];
sessionTwoData = [
11002	1	1	4.1835
11002	2	1	4.1134
11002	3	1	4.1486
11002	4	1	3.1766
11002	5	1	4.1554
11057	1	1	3.3216
11057	2	1	3.3796
11057	3	1	3.302
11057	4	1	3.4218
11057	5	1	3.605
11108	1	1	3.2883
11108	2	1	3.2995
11108	3	1	3.253
11108	4	1	3.115
11108	5	1	3.1489
11112	1	1	2.6928
11112	2	1	3.0409
11112	3	1	2.5383
11112	4	1	2.5951
11112	5	1	2.6928
11002	1	2	4.5785
11002	2	2	4.5646
11002	3	2	4.1681
11002	4	2	4.0477
11002	5	2	4.0164
11057	1	2	3.5662
11057	2	2	3.59
11057	3	2	3.1557
11057	4	2	3.1134
11057	5	2	3.58
11108	1	2	2.911
11108	2	2	2.8379
11108	3	2	3.1657
11108	4	2	3.3003
11108	5	2	3.7035
11112	1	2	2.7676
11112	2	2	2.2237
11112	3	2	2.6068
11112	4	2	3.139
11112	5	2	2.4181];

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

%% Run the one way for session 1
varNames = strvcat('Subject_1', 'Trial_1', 'Group_1');
subject1 = sessionOneData(:,1);
trial1 = sessionOneData(:,2);
group1 = sessionOneData(:,3);
peakReflectance1 = sessionOneData(:,4);

% Use model to specify one-way.  Not really
% sure if trial and group should be random, since
% they are ordered.  But it doesn't change the ANOVA
% whether they are specified as random or not.  Leaving
% only subject as random.
random = [1];
model = [0 0 1];
[pFull, tabFull] = anovan(peakReflectance1',{subject1 trial1 group1}, 'model',model,'varnames', varNames, 'random', random);

%% Run the one way for session 2
varNames = strvcat('Subject_2', 'Trial_2', 'Group_2');
subject2 = sessionTwoData(:,1);
trial2 = sessionTwoData(:,2);
group2 = sessionTwoData(:,3);
peakReflectance2 = sessionTwoData(:,4);

% Use model to specify one-way
random = [1];
model = [0 0 1];
[pFull, tabFull] = anovan(peakReflectance2',{subject2 trial2 group2}, 'model',model,'varnames', varNames, 'random', random);

%% Specify a 2-way model that focusses on stimulus and session.
%
% Considering session random as well as subject.
random = [1 4];
varNames = strvcat('Subject', 'Trial', 'Group', 'Session');
subject = [subject1 ; subject2];
trial = [trial1 ; trial2];
group = [group1 ; group2];
session1 = ones(size(group1));
session2 = 2*ones(size(group2));
session = [session1 ; session2];
peakReflectance = [peakReflectance1 ; peakReflectance2];
model = [0 0 1 0; 0 0 0 1; 0 0 1 1];
[pTwoWay, tabTwoWay] = anovan(peakReflectance',{subject trial group session}, 'model',model,'varnames', varNames, 'random', random);

