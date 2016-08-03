function [max_resp_ampl] = modelFit(timeBase, pooled_std_stim)
% gammaFitTutorial
%
% Illustrates using fmincon to fit a gamma function to data.
%
% 12/31/15 dhb       Wrote from parameter search tutorial.
% 12/31/15 rfc      Added temporal analysis specific parameters.

%% Initialize
% close all;

%% Generate some simulated data for fitting
trueParams.type = '2xgammapdf';
trueParams.preStimValue = 0;
trueParams.stimOnsetTime = 2;
trueParams.responseDelay1 = 0.4;
trueParams.scale1 = 3;
trueParams.gammaA1 = 4;
trueParams.gammaB1 = 0.24;
trueParams.responseDelay2 = 1.5;
trueParams.scale2 = 3;
trueParams.gammaA2 = 4;
trueParams.gammaB2 = 0.24;

startTime = 0;
endTime = max(timeBase);
% nTrueData = 200;
% timeBase = linspace(startTime,endTime,nTrueData);
% trueParams.noiseSd = 0;
% theResponse = ComputeModelPreds(trueParams,timeBase);
% trueParams.noiseSd = 0.3;
% pooled_std_stim = ComputeModelPreds(trueParams,timeBase);

%Remove values below zero


%% Start plot
thePlot = figure(2); clf; hold on
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
fitParams0.type = '2xgammapdf';

fitParams0.stimOnsetTime = 3.96;

% Remove any nan.
maskout = ~isnan(pooled_std_stim);

pooled_std_stim = pooled_std_stim(maskout);
timeBase = timeBase(maskout);

fitParams0.preStimValue = mean( pooled_std_stim( timeBase < fitParams0.stimOnsetTime ) );

% These we make up based on our excellent judgement
fitParams0.responseDelay1 = 0;
fitParams0.gammaA1 = 2;
fitParams0.responseDelay2 = 2;
fitParams0.gammaA2 = 2;

% Then rate parameter so that mode of distribution happens at
% the time of maximum response;
[maxResp,maxRespIndex] = max(pooled_std_stim);
maxTime = timeBase(maxRespIndex(1));

fitParams0.gammaB1 = (fitParams0.gammaA1-1)/(maxTime - fitParams0.stimOnsetTime);
if (fitParams0.gammaB1 <= 0)
   fitParams0.gammaB1 = 0.5;
end
fitParams0.gammaB2 = (fitParams0.gammaA2-1)/(fitParams0.responseDelay2 - fitParams0.stimOnsetTime);
if (fitParams0.gammaB2 <= 0)
   fitParams0.gammaB2 = 0.5;
end

% And set scale so that initial max predition matches max data
fitParams0.scale1 = 1;
fitParams0.scale2 = 1;
tempPreds = ComputeModelPreds(fitParams0,timeBase);
tempMax = max(tempPreds);
fitParams0.scale1 = maxResp/tempMax;
fitParams0.scale2 = maxResp/(tempMax*2);

% Add initial guess to the plot
predictions0 = ComputeModelPreds(fitParams0,timeBase);
figure(thePlot); plot(timeBase,predictions0,'k:','LineWidth',2);

%% Fit

% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');

% Initial guess
x0 = ParamsToX(fitParams0);
% fieldnames(fitParams0);
% First seach on gammaA and scale only, and add to plot
vlb = [x0(1) 0.01 x0(3) 0.01 x0(5) 0.01 x0(7) -100];
vub = [x0(1) 100 x0(3) 100 x0(5) 100 x0(7) 100];
x1 = fmincon(@(x)FitModelErrorFunction(x,timeBase,pooled_std_stim,fitParams0),x0,[],[],[],[],vlb,vub,[],options);
predictions1 = ComputeModelPreds(XToParams(x1,fitParams0),timeBase);
figure(thePlot); plot(timeBase,predictions1,'c:','LineWidth',2);

% Then full search
vlb = [0 0.01 0.01 0.1 .5 0.01 0.01 0.1];
vub = [100 100 100 1000 100 100 100 1000];
x = fmincon(@(x)FitModelErrorFunction(x,timeBase,pooled_std_stim,fitParams0),x1,[],[],[],[],vlb,vub,[],options);

% Extract fit parameters
fitParams = XToParams(x,fitParams0)

% Add final fit to plot
predictions = ComputeModelPreds(fitParams,timeBase);
figure(thePlot); plot(timeBase,predictions,'g','LineWidth',2);
legend({' Data', 'Underlying Fcn' ' Initial Guess', ' Intermediate Fit' ' Final Fit'},'FontSize',14,'Location','NorthEast');

figure(3); plot(timeBase,pooled_std_stim-predictions');
title('Residuals');

[max_ampl, max_ind ] = max(predictions);

max_resp_ampl = max_ampl-fitParams0.preStimValue

max_prestim_val = max(pooled_std_stim( timeBase < fitParams0.stimOnsetTime ) )

resp_start_time = timeBase( min( find( predictions > max_prestim_val ) ) )

time_to_peak  = timeBase(max_ind) - fitParams0.stimOnsetTime


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
    case '2xgammapdf'
        x = [params.responseDelay1 params.gammaA1 params.gammaB1 params.scale1 params.responseDelay2 params.gammaA2 params.gammaB2 params.scale2];
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
        params.responseDelay = x(1);
        params.gammaA = x(2);
        params.gammaB = x(3);
        params.scale = x(4);
    case '2xgammapdf'
        params.responseDelay1 = x(1);
        params.gammaA1 = x(2);
        params.gammaB1 = x(3);
        params.scale1 = x(4);
        params.responseDelay2 = x(5);
        params.gammaA2 = x(6);
        params.gammaB2 = x(7);
        params.scale2 = x(8);
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
        firstDelayZeroedTime = stimZeroedTime-params.responseDelay;
        index = find(firstDelayZeroedTime >= 0);
        
        preds(index) = preds(index)+params.scale*gampdf(firstDelayZeroedTime(index),params.gammaA,params.gammaB);
    case '2xgammapdf'
        preds = params.preStimValue*ones(size(timeBase));
        
        stimZeroedTime = timeBase-params.stimOnsetTime;
        firstDelayZeroedTime = stimZeroedTime-params.responseDelay1;
        secondDelayZeroedTime = stimZeroedTime-params.responseDelay2;
        
        index1 = find(firstDelayZeroedTime >= 0); % Don't use find?
        index2 = find(secondDelayZeroedTime >= 0);
        
        preds(index1) = preds(index1) + params.scale1*gampdf(firstDelayZeroedTime(index1),params.gammaA1,params.gammaB1);
        preds(index2) = preds(index2) + params.scale2*gampdf(secondDelayZeroedTime(index2),params.gammaA2,params.gammaB2);

    otherwise
        error('Unknown model type');
end

% Add simulated noise, optionally
% if (isfield(params,'noiseSd') & params.noiseSd > 0)
%     preds = preds + normrnd(0,params.noiseSd,size(preds));
% end

end
