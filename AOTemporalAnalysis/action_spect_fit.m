function [ rel_shifts ] = action_spect_fit( wavelengths, irradiances, datavals, type )

normalize = false;

% Set up the structure
fitParams.max_or_intercept = [];
fitParams.irr_shift=[];
fitParams.slope=0;

max_irr = 10^5;
log_irr_range = 0:.25:11;


% Rescale amplitudes so that we can fit them all the same- each column is
% a different wavelength
handle = figure;

for w=1:length(wavelengths)
    data = datavals{w};
    irr = irradiances{w};
    
    irr = irr(~isinf(data));
    data = data(~isinf(data));
    datavals{w} = data;
    irradiances{w} = irr;
    
    if wavelengths(w) ~= 480 && normalize
        datavals{w} = data./max(data);
    end
    switch(type)
        case 'sigmoid'
            fitParams.max_or_intercept = [fitParams.max_or_intercept data(end)];       
            irrval(w) = {log(irr)};
            fitParams.slope = fitParams.slope+1;
        case 'linear'
            irr=log(irr);
            m = (data(end)-data(1))/(irr(end)-irr(1));
            fitParams.slope = fitParams.slope+m;
            
            c = data(1)-(m*irr(1));
            
            fitParams.max_or_intercept = [fitParams.max_or_intercept c];
            irrval(w) = {irr};
        case 'exponential'
            fitParams.max_or_intercept =[fitParams.max_or_intercept data(end)];
            m = (data(end)-data(1))/(irr(end)-irr(1));
            fitParams.slope = fitParams.slope+m;
            
            irrval(w) = {log(irr)};
    end
    
    fitParams.irr_shift = [fitParams.irr_shift mean(irrval{w})];
    
    % The midpoint of the fit should never move beyond the edges of the
    % actual curve
    vub(w) = max(irrval{w}); 
    vlb(w) = min(irrval{w});
    
    plot(irrval{w}, datavals{w}); hold on;
    legends{w} = num2str(wavelengths(w));
end
legend(legends);

% Determine the average slope we want to use.
fitParams.slope = fitParams.slope/length(wavelengths);

% The starting maximum value is the median of our maximum amplitudes
fitParams.max_or_intercept = median(fitParams.max_or_intercept);


for i=1:size(datavals,2)
   
    thisParam.max_or_intercept = fitParams.max_or_intercept;
    thisParam.slope = fitParams.slope;
    thisParam.midpoint = fitParams.irr_shift(i);
        
%     plot(log_irr_range, ComputeModel(thisParam, log_irr_range),'k.');
end
hold off;


options = optimset('fmincon');
options = optimset(options,'LargeScale','off','Algorithm','interior-point'); % 'Diagnostics','off','Display','off'

switch(type)
    case 'sigmoid'
        vub = [3 3 vub];
        vlb = [.5 0 vlb];
    case 'linear'
        vub = [3 2 vub+2];
        vlb = [-3 -2 vlb-2];
    case 'exponential'
        vub = [3 3 vub];
        vlb = [.2 -3 vlb];
end

v0 = ParamsToVector(fitParams);

finalParams = fmincon(@(v)FitModelErrorFunction(v,irrval,datavals,type),v0,[],[],[],[],vlb,vub,[],options);

%% Lock the parameters, then shift the curves up and down to fit the CIE

rel_shifts = 1./( exp(finalParams(3:end))./exp(finalParams(5)) );

% Load CIE data
ciefunc = dlmread('/local_data/Projects/AdaptiveOpticsAnalysis/AOTemporalAnalysis/linCIE2008v2e_5.csv');

cierow = [];
for w=1:length(wavelengths)    
    cierow = [cierow find(ciefunc(:,1) == wavelengths(w))];
end

cievals = ciefunc(cierow,2)';

vub = 0.6;
vlb = min(rel_shifts)-0.02; % Can't shift below what is possible, nor can it be 0!

v0 = cievals(3)-1;
    
finalShifts = fmincon(@(v)FitCIEErrorFunction(v,cievals,rel_shifts),v0,[],[],[],[],vlb,vub,[],options);


rel_shifts= finalShifts+rel_shifts;

hold on;
for i=1:size(datavals,2)
    params = VectorToParams(finalParams, i)
    plot(log_irr_range, ComputeModel(params, log_irr_range, type),'r.-');
    
end
hold off;
axis([0 max(cellfun(@max,irrval)) 0 max(cellfun(@max,datavals))]);


end

function f = FitCIEErrorFunction(v, cie, shifts)

    preds = ComputeCIEModel(v,shifts);
    
    diff = (preds-cie).^2;

    f = 100*sqrt(sum(diff)/length(diff));
end

function values = ComputeCIEModel(v, shifts)

values = shifts+v;

end

function f = FitModelErrorFunction(v, irradiances, values,type)

    diff=[];
    for i=1:size(values,2)

        params = VectorToParams(v,i);

        preds{i} = ComputeModel(params, irradiances{i},type);

        diff = [diff; (values{i}-preds{i}).^2];
    end
    
    f = 100*sqrt(sum(diff)/length(diff));
end


function v = ParamsToVector(params)
    v = [params.max_or_intercept, params.slope, params.irr_shift];
end

function params = VectorToParams(v, shiftnum)
    params.max_or_intercept = v(1);
    params.slope = v(2);
    params.midpoint = v(2+shiftnum);
end

function values = ComputeModel(params, irrads, type)
    switch(type)
        case 'sigmoid'
            values = params.max_or_intercept./( 1+exp(-params.slope*(irrads-params.midpoint)) );
        case 'linear'
            values = params.slope*(irrads-params.midpoint)+params.max_or_intercept;
        case 'exponential'
            values = params.max_or_intercept*exp( -params.slope * (irrads-params.midpoint) );        
            
    end
end