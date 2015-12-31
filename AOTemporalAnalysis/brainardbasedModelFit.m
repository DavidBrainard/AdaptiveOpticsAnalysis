function brainardbasedModelFit(timeBase, pooled_std_stim)
% gammaFitTutorial
%
% Illustrates using fmincon to fit a gamma function to data.
%
% 12/31/15 dhb       Wrote from parameter search tutorial.
% 12/31/15 rfc      Added temporal analysis specific parameters.

%% Initialize
close all;

%% Generate some simulated data for fitting
trueParams.type = 'gammapdf';
trueParams.preStimValue = 0;
trueParams.stimOnsetTime = 2;
trueParams.responseDelay = 0.4;
trueParams.scale = 3;
trueParams.gammaA = 4;
trueParams.gammaB = 0.24;

startTime = 0;
endTime = max(timeBase);
% nTrueData = 200;
% timeBase = linspace(startTime,endTime,nTrueData);
% trueParams.noiseSd = 0;
% theResponse = ComputeModelPreds(trueParams,timeBase);
% trueParams.noiseSd = 0.3;
% pooled_std_stim = ComputeModelPreds(trueParams,timeBase);

%% Start plot
thePlot = figure; clf; hold on
set(gca,'FontName','Helvetica','FontSize',14);
plot(timeBase,pooled_std_stim,'ro','MarkerFaceColor','r','MarkerSize',6);
% figure(thePlot); plot(timeBase,theResponse,'r','LineWidth',4);
xlim([startTime endTime]);
ylim([-1 2]);
xlabel('Time (secs)','FontSize',18);
ylabel('Pooled Standard deviation','FontSize',18);
title('Pooled standard deviation data and fit');
drawnow;

%% Set up initial guess for fit parameters

% These are known
fitParams0.type = 'gammapdf';
fitParams0.preStimValue = -0.1;
fitParams0.stimOnsetTime = 5;

% These we make up based on our excellent judgement
fitParams0.responseDelay = 0;
fitParams0.gammaA = 2;

% Then rate parameter so that mode of distribution happens at
% the time of maximum response;
[maxResp,maxRespIndex] = max(pooled_std_stim);
maxTime = timeBase(maxRespIndex(1));
fitParams0.gammaB = (fitParams0.gammaA-1)/(maxTime - fitParams0.stimOnsetTime);
if (fitParams0.gammaB <= 0)
   fitParams0.gammaB = 0.5;
end

% And set scale so that initial max predition matches max data
fitParams0.scale = 1;
tempPreds = ComputeModelPreds(fitParams0,timeBase);
tempMax = max(tempPreds);
fitParams0.scale = maxResp/tempMax;

% Add initial guess to the plot
predictions0 = ComputeModelPreds(fitParams0,timeBase);
figure(thePlot); plot(timeBase,predictions0,'k:','LineWidth',2);

%% Fit

% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');

% Initial guess
x0 = ParamsToX(fitParams0);

% First seach on gammaA and scale only, and add to plot
vlb = [x0(1) 0.01 x0(3) 0.01];
vub = [x0(1) 100 x0(3) 100];
x1 = fmincon(@(x)FitModelErrorFunction(x,timeBase,pooled_std_stim,fitParams0),x0,[],[],[],[],vlb,vub,[],options);
predictions1 = ComputeModelPreds(XToParams(x1,fitParams0),timeBase);
figure(thePlot); plot(timeBase,predictions1,'c:','LineWidth',2);

% Then full search
vlb = [0 0.01 0.01 0.1];
vub = [2 100 100 1000];
x = fmincon(@(x)FitModelErrorFunction(x,timeBase,pooled_std_stim,fitParams0),x1,[],[],[],[],vlb,vub,[],options);

% Extract fit parameters
fitParams = XToParams(x,fitParams0);

% Add final fit to plot
predictions = ComputeModelPreds(fitParams,timeBase);
figure(thePlot); plot(timeBase,predictions,'g','LineWidth',2);
legend({' Data', 'Underlying Fcn' ' Initial Guess', ' Intermediate Fit' ' Final Fit'},'FontSize',14,'Location','NorthEast');

end

% f = FitModelErrorFunction(x,timeBase,theResponse,fitParams)
%
% Search error function
function f = FitModelErrorFunction(x,timeBase,theResponse,fitParams)

% Extract parameters into meaningful structure
fitParams = XToParams(x,fitParams);

% Make predictions
preds = ComputeModelPreds(fitParams,timeBase)';

% Compute fit error as RMSE
nPoints = length(theResponse);
theDiff2 = (theResponse-preds).^2;
f = 100*sqrt(sum(theDiff2)/nPoints);

end

% x = ParamsToX(params)
%
% Convert parameter structure to vector of parameters to search over
function x = ParamsToX(params)
switch (params.type)
    case 'gammapdf'
        x = [params.responseDelay params.gammaA params.gammaB params.scale];
    otherwise
        error('Unknown model type');
end
end


% fitParams = XToParams(x,params)
%
% Convert search params and base structure to filled in structure.
function params = XToParams(x,params)
switch (params.type)
    case 'gammapdf'
        params.delay = x(1);
        params.gammaA = x(2);
        params.gammaB = x(3);
        params.scale = x(4);
    otherwise
        error('Unknown model type');
end
end

% preds =  ComputeModelPreds(params,t)
%
% Compute the predictions of the model
function preds = ComputeModelPreds(params,timeBase)

% Allow multiple model types
switch (params.type)
    case 'gammapdf'
        preds = params.preStimValue*ones(size(timeBase));
        stimZeroedTime = timeBase-params.stimOnsetTime;
        delayZeroedTime = stimZeroedTime-params.responseDelay;
        index = find(delayZeroedTime >= 0);
        preds(index) = params.scale*gampdf(delayZeroedTime(index),params.gammaA,params.gammaB);
    otherwise
        error('Unknown model type');
end

% Add simulated noise, optionally
% if (isfield(params,'noiseSd') & params.noiseSd > 0)
%     preds = preds + normrnd(0,params.noiseSd,size(preds));
% end

end
