function [fitCharacteristics, residual] = modelFit_beta(timeBase, pooled_std_stim, plotstuff)
% [fitCharacteristics] = modelFit_beta(timeBase, pooled_std_stim)
%
%
% 12/31/15 dhb       Wrote from parameter search tutorial.
% 12/31/15 rfc      Added temporal analysis specific parameters.
%
% @params:
%    timeBase: A 1xN array of time stamps corresponding pooled_std_stim
%
%    pooled_std_stim: The subtracted pooled standard deviation of the 
%                     stimulus and control videos
%
% @outputs:
%    fitCharacteristics: A struct containing information about the fit
%                       the contents of this struct are subject to change.



%% Generate some simulated data for fitting
if ~exist('plotstuff','var')
    plotstuff = false;
end

%Keep values before cutoff time
cutofftime = 7;

pretime = timeBase <= cutofftime;
timeBase = timeBase(pretime);
pooled_std_stim = pooled_std_stim(pretime);


% pooled_std_stim = wavelet_denoise(pooled_std_stim);


%% Start plot
if plotstuff
    thePlot = figure(2); hold on; clf; 
    % set(gca,'FontName','Helvetica','FontSize',14);
    plot(timeBase,pooled_std_stim,'r','MarkerFaceColor','r','MarkerSize',6);
    % figure(thePlot); plot(timeBase,theResponse,'r','LineWidth',4);
    xlim([0 16]);
    ylim([-1 3]);
    xlabel('Time (secs)','FontSize',18);
    ylabel('Pooled Standard deviation','FontSize',18);
    title('Pooled standard deviation data and fit');
    drawnow;
end


%% Set up initial guess for fit parameters

% These are known
fitParams0.type = 'gammapdf';

fitParams0.stimOnsetTime = 3.96;

% Remove any nan.
maskout = ~isnan(pooled_std_stim);

pooled_std_stim = pooled_std_stim(maskout);
timeBase = timeBase(maskout);

fitParams0.preStimValue = mean( pooled_std_stim( timeBase < fitParams0.stimOnsetTime ) );

% Then rate parameter so that mode of distribution happens at
% the time of maximum response;
[maxResp,maxRespIndex] = max(pooled_std_stim);
maxTime = timeBase(maxRespIndex(1));
preResp = pooled_std_stim(find(timeBase <= fitParams0.stimOnsetTime+.01 & timeBase >= fitParams0.stimOnsetTime-.01));

max_prestim_val = max(pooled_std_stim( timeBase < fitParams0.stimOnsetTime ) );

prestim_stddev = std(pooled_std_stim( timeBase < fitParams0.stimOnsetTime ) );

prestim_PI = fitParams0.preStimValue + prestim_stddev;

% Kludgy setup, but works for a simple analysis...
global peaked;

peaked = true;

% And set scale so that initial max predition matches max data
if ~isempty( strfind( fitParams0.type, 'gammapdf') )
    
    overPreStim = find( pooled_std_stim(timeBase > fitParams0.stimOnsetTime) > prestim_PI );
    
    % Supress the gamma function if there is no obvious response.
    if length(overPreStim)>15
        fitParams0.scale1 = 1;
        fitParams0.gammaA1 = 1.25;
        fitParams0.gammaB1 = (maxTime - fitParams0.stimOnsetTime)/(fitParams0.gammaA1);
        if (fitParams0.gammaB1 <= 1)
           fitParams0.gammaB1 = 1;
        end
        peaked = true;
    else
        peaked = false;
        
        fitParams0.gammaA1 = 1;
        fitParams0.gammaB1 = 1;
        fitParams0.scale1 = 0.01;
    end
        
        fitParams0.responseDelay1 = .15;
    
    if strcmp( fitParams0.type, '2xgammapdf')
        fitParams0.responseDelay2 = 2;
        fitParams0.gammaA2 = 2;
        
        fitParams0.gammaB2 = (fitParams0.gammaA2-1)/(fitParams0.responseDelay2 - fitParams0.stimOnsetTime);
        if (fitParams0.gammaB2 <= 2)
           fitParams0.gammaB2 = 2;
        end
        
        fitParams0.scale2 = 1;
        fitParams0.offset = 0; %min(pooled_std_stim);
    end
    
    if strcmp( fitParams0.type, 'gammapdfexp')

        fitParams0.responseDelay2 = 0;
        fitParams0.decay   = .25;
        fitParams0.offset = max(pooled_std_stim)-min(pooled_std_stim); 
    end
    
    if peaked
        tempPreds = ComputeModelPreds(fitParams0,timeBase);
        tempMax = max(tempPreds);
        fitParams0.scale1 = maxResp/tempMax;
    end
    
    if strcmp( fitParams0.type, '2xgammapdf')
%         poststimind = find(timeBase>fitParams0.stimOnsetTime);
%         [minval, minind] = min(pooled_std_stim(poststimind(10:end)));
%         minind = minind+ poststimind(10)-1;
%         if minval > 0            
            fitParams0.scale2 = maxResp/(tempMax*2);
%         else
%             fitParams0.responseDelay2 = timeBase(minind)-fitParams0.stimOnsetTime;
%             fitParams0.scale2 = -maxResp/(tempMax*2);
%         end
    end

end
    
% fitParams0;
% Add initial guess to the plot
predictions0 = ComputeModelPreds(fitParams0,timeBase);
if plotstuff
    figure(thePlot); hold on; plot(timeBase,predictions0,'k:','LineWidth',2); hold off;
end

%% Fit

% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','sqp');

% Initial guess
x0 = ParamsToX(fitParams0);

% Then full search
switch fitParams0.type
    case 'gammapdf'
        vlb = [-2 0 0 -10];
        vub = [2  3 3 10];        
    case '2xgammapdf'
        vlb = [-2 0.1 0.1 0.001  0  1    1  -6   0];
        vub = [ 1 10    10    10   2  10   10  6   0];
    case 'gammapdfexp'
        if peaked
            vlb = [-3 0 0 -20  0.01 0];
            vub = [ 3 3 3  20  3    3];
        else
            vlb = [-3 0 0 -20  0.01 0];
            vub = [ 3 3 3  0.5 3    3];
        end
end

x = fmincon(@(x)FitModelErrorFunction(x,timeBase,pooled_std_stim,fitParams0),x0,[],[],[],[],vlb,vub,[],options);

% Extract fit parameters
fitParams = XToParams(x,fitParams0);


delta=min(diff(timeBase));
fitTimeBase = 0:delta:max(timeBase);

% Add final fit to plot
predictions = ComputeModelPreds(fitParams,fitTimeBase);

residual = mean((pooled_std_stim-predictions(1:length(pooled_std_stim))').^2);

if plotstuff
    figure(thePlot); hold on; plot(fitTimeBase,predictions,'g','LineWidth',2);
end

if peaked
    [max_ampl, max_ind ] = max((predictions));
else
    [max_ampl, max_ind ] = max(abs(predictions)); 
end

threeQind = min( find( predictions > (max_ampl/2) ) );
threeQval = predictions(threeQind);

% Interpolate to find the exact spot it becomes greater than the 2sd+mean of the prestim
% value- should help with quantization issues
prestim_PI = fitParams0.preStimValue + 2*prestim_stddev;

fitCharacteristics.amplitude = predictions(max_ind)-fitParams0.preStimValue;

afterPIval = min( find( abs(predictions) > prestim_PI ) );

fitCharacteristics.resp_start = 0;
fitCharacteristics.time_to_peak = fitTimeBase(max_ind) - fitParams0.stimOnsetTime;

if fitCharacteristics.time_to_peak < 0
    fitCharacteristics.time_to_peak = Inf;
end
%     figure; plot(fitTimeBase,predictions); hold on;
%     plot(timeBase, pooled_std_stim,'k')
%     plot(fitTimeBase(max_ind), predictions(max_ind), '*r');
%     xlabel('Time'); ylabel('Response');
%     hold off;



if afterPIval<=length(timeBase)            
    beforePIval = afterPIval-1;

    interpslope = (predictions(afterPIval)-predictions(beforePIval))/(timeBase(afterPIval)-timeBase(beforePIval));

    fitCharacteristics.resp_start = timeBase(beforePIval) + ((prestim_PI-predictions(beforePIval))/interpslope);
end


if strcmp( fitParams0.type, 'gammapdfexp') 
    fitCharacteristics.decay_initval = fitParams.offset;
end

if exist('fitParams.decay','var')
    fitCharacteristics.decay_constant = fitParams.decay;
end

if strcmp( fitParams0.type, '2xgammapdf') 
    
    tmpparam1.responseDelay1 = fitParams.responseDelay1;
    tmpparam1.gammaA1 = fitParams.gammaA1;
    tmpparam1.gammaB1 = fitParams.gammaB1;
    tmpparam1.scale1 = fitParams.scale1;
    tmpparam1.type = 'gammapdf';
    tmpparam1.preStimValue = fitParams0.preStimValue;
    tmpparam1.stimOnsetTime = fitParams0.stimOnsetTime;

    tmpparam2.responseDelay1 = fitParams.responseDelay2;
    tmpparam2.gammaA1 = fitParams.gammaA2;
    tmpparam2.gammaB1 = fitParams.gammaB2;
    tmpparam2.scale1 = fitParams.scale2;
    tmpparam2.type = 'gammapdf';
    tmpparam2.preStimValue = fitParams0.preStimValue;
    tmpparam2.stimOnsetTime = fitParams0.stimOnsetTime;
    
    gamma1preds = ComputeModelPreds(tmpparam1, timeBase);
    gamma2preds = fitParams.offset + ComputeModelPreds(tmpparam2, timeBase);
    
    [pks loc1]=findpeaks(gamma1preds,timeBase);
    [~, loc2]=findpeaks(gamma2preds,timeBase);
    
%     figure(thePlot); hold on; plot(timeBase,gamma1preds,'k', timeBase,gamma2preds,'b'); hold off;
    
    fitCharacteristics.gamma_separation = abs(loc1-loc2);
end

fitParams;
% afterPIval
% threeQind

% if threeQind ~= afterPIval
%     fitCharacteristics.max_slope = (threeQval-predictions(afterPIval))/(timeBase(threeQind)-timeBase(afterPIval))
% else
    fitCharacteristics.max_slope = max(diff(predictions))/(timeBase(2)-timeBase(1));
% end

% figure(thePlot); hold on; 
%     plot(fitCharacteristics.resp_start,1.5,'g*');
%     plot(timeBase(max_ind),1.5,'b*');
% hold off;


end

% f = FitModelErrorFunction(x,timeBase,theResponse,fitParams)
%
% Search error function
function f = FitModelErrorFunction(x,timeBase,theResponse,fitParams)

% Extract parameters into meaningful structure
fitParams = XToParams(x,fitParams);

% Make predictions
preds = ComputeModelPreds(fitParams,timeBase);

% Compute fit error as RMSE
nPoints = length(theResponse);
% theDiff2 = (theResponse-preds).^2;
% f = mean(sqrt(theDiff2));
theDiff2 = abs(theResponse-preds);
f = mean(theDiff2);

end

% x = ParamsToX(params)
%
% Convert parameter structure to vector of parameters to search over
function x = ParamsToX(params)
switch (params.type)
    case 'gammapdf'
        x = [params.responseDelay1 params.gammaA1 params.gammaB1 params.scale1];
    case '2xgammapdf'
        x = [params.responseDelay1 params.gammaA1 params.gammaB1 params.scale1 params.responseDelay2 params.gammaA2 params.gammaB2 params.scale2 params.offset];
    case 'gammapdfexp'
%         x = [params.responseDelay1 params.gammaA1 params.gammaB1 params.scale1 params.responseDelay2 params.decay params.offset];
        x = [params.responseDelay1 params.gammaA1 params.gammaB1 params.scale1 params.decay params.offset];
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
        params.responseDelay1 = x(1);
        params.gammaA1 = x(2);
        params.gammaB1 = x(3);
        params.scale1 = x(4);
    case '2xgammapdf'
        params.responseDelay1 = x(1);
        params.gammaA1 = x(2);
        params.gammaB1 = x(3);
        params.scale1 = x(4);
        params.responseDelay2 = x(5);
        params.gammaA2 = x(6);
        params.gammaB2 = x(7);
        params.scale2 = x(8);
        params.offset = x(9);
    case 'gammapdfexp'
        params.responseDelay1 = x(1);
        params.gammaA1 = x(2);
        params.gammaB1 = x(3);
        params.scale1  = x(4);
%         params.responseDelay2 = x(5);
        params.decay   = x(5);        
        params.offset  = x(6);        
    otherwise
        error('Unknown model type');
end
end

% preds =  ComputeModelPreds(params,t)
%
% Compute the predictions of the model
function preds = ComputeModelPreds(params,timeBase)

global peaked;

% Allow multiple model types
switch (params.type)
    case 'gammapdf'
        preds = params.preStimValue*ones(size(timeBase));
        stimZeroedTime = timeBase-params.stimOnsetTime;
        firstDelayZeroedTime = stimZeroedTime-params.responseDelay1;
        index = find(firstDelayZeroedTime >= 0);
        
        preds(index) = preds(index)+params.scale1*gampdf(firstDelayZeroedTime(index),params.gammaA1,params.gammaB1);

    case '2xgammapdf'
        preds = params.preStimValue*ones(size(timeBase));
        
        stimZeroedTime = timeBase-params.stimOnsetTime;
        firstDelayZeroedTime = stimZeroedTime-params.responseDelay1;
        secondDelayZeroedTime = stimZeroedTime-params.responseDelay2;
        
        index1 = find(firstDelayZeroedTime >= 0); % Don't use find?
        index2 = find(secondDelayZeroedTime >= 0);
        
        preds(index1) = preds(index1) + params.scale1*gampdf(firstDelayZeroedTime(index1),params.gammaA1,params.gammaB1);
        
        preds(index2) = preds(index2)+ params.offset + params.scale2*gampdf(secondDelayZeroedTime(index2),params.gammaA2,params.gammaB2);
        
        % Clamp the pre-stim region to the prestim value.
%         preds(find(firstDelayZeroedTime <= 0)) = params.preStimValue;
    case 'gammapdfexp'
        preds = params.preStimValue*ones(size(timeBase));
        stimZeroedTime = timeBase-params.stimOnsetTime;

        firstDelayZeroedTime  = stimZeroedTime-params.responseDelay1;
        
        index1 = find(firstDelayZeroedTime >= 0); % Don't use find?
%         if peaked        
            preds(index1) = preds(index1) + params.scale1*gampdf(firstDelayZeroedTime(index1),params.gammaA1,params.gammaB1);

        
            % Using maxval as the anchor
            [maxval, ind] = max(preds);
            secondDelayZeroedTime = timeBase-timeBase(ind(1));          
            index2 = find( secondDelayZeroedTime >= 0 );            
%         else
%             % Using maxval as the anchor
%             maxval= preds(index1(1));
%             secondDelayZeroedTime = firstDelayZeroedTime;
%             index2 = index1;
%         end
        
        % The exponential must always line up with the gamma function's
        % value!

        maxmatch = maxval-max(params.offset*exp( -params.decay * secondDelayZeroedTime(index2) ));
        
%         maxmatch = preds(index2(1))-params.offset;
        
        preds(index2) = maxmatch + params.offset*exp( -params.decay * secondDelayZeroedTime(index2) );
        
%         figure(400);plot(maxmatch, preds(index2(1)), 'r*'); hold on; drawnow;
    otherwise
        error('Unknown model type');
end

% Add simulated noise, optionally
% if (isfield(params,'noiseSd') & params.noiseSd > 0)
%     preds = preds + normrnd(0,params.noiseSd,size(preds));
% end

end
