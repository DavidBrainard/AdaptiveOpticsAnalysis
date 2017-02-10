function [ rel_shifts ] = action_spect_fit( wavelengths, irradiances, amplitudes )


normalize = false;

% Set up the structure
fitParams.max_val = [];
fitParams.midpoint_shifts=[];

max_irr = 10^5;
log_irr_range = 0:.25:11;


% Rescale amplitudes so that we can fit them all the same- each column is
% a different wavelength
handle = figure;

for w=1:length(wavelengths)
    amps = amplitudes{w};
    
    if wavelengths(w) ~= 480 && normalize
        amplitudes{w} = amps./max(amps);
    end
    
    fitParams.max_val = [fitParams.max_val amps(end)];    
    
    log_irr(w) = {log(irradiances{w})};
    
    fitParams.midpoint_shifts = [fitParams.midpoint_shifts mean(log_irr{w})];
    
    % The midpoint of the sigmoid should never move beyond the edges of the
    % actual curve
    vub(w) = max(log_irr{w}); 
    vlb(w) = min(log_irr{w});
    
    plot(log_irr{w}, amplitudes{w}); hold on;
    legends{w} = num2str(wavelengths(w));
end
 legend(legends);

% The starting maximum value is the median of our maximum amplitudes
fitParams.max_val = median(fitParams.max_val);
fitParams.steepness = 1;


for i=1:size(amplitudes,2)
   
    thisParam.max_val = fitParams.max_val;
    thisParam.steepness = fitParams.steepness;
    thisParam.midpoint = fitParams.midpoint_shifts(i);
        
%     plot(log_irr_range, ComputeModel(thisParam, log_irr_range),'k.');
end
hold off;


options = optimset('fmincon');
options = optimset(options,'LargeScale','off','Algorithm','interior-point'); % 'Diagnostics','off','Display','off'

vub = [3 3 vub];
vlb = [.5 0 vlb];

v0 = ParamsToVector(fitParams);

finalParams = fmincon(@(v)FitModelErrorFunction(v,log_irr,amplitudes),v0,[],[],[],[],vlb,vub,[],options);

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
vlb = min(rel_shifts)-0.02 % Can't shift below what is possible, nor can it be 0!

v0 = cievals(3)-1;
    
finalShifts = fmincon(@(v)FitCIEErrorFunction(v,cievals,rel_shifts),v0,[],[],[],[],vlb,vub,[],options);

finalShifts
rel_shifts= finalShifts+rel_shifts;

hold on;
for i=1:size(amplitudes,2)
   
    plot(log_irr_range, ComputeModel(VectorToParams(finalParams, i), log_irr_range),'r.-');
end
hold off;


end

function f = FitCIEErrorFunction(v, cie, shifts)

    preds = ComputeCIEModel(v,shifts);
    
    diff = (preds-cie).^2;

    f = 100*sqrt(sum(diff)/length(diff));
end

function values = ComputeCIEModel(v, shifts)

values = shifts+v;

end

function f = FitModelErrorFunction(v, irradiances, amplitudes)

    diff=[];
    for i=1:size(amplitudes,2)

        params = VectorToParams(v,i);

        preds{i} = ComputeModel(params, irradiances{i});

        diff = [diff; (amplitudes{i}-preds{i}).^2];
    end
    
    f = 100*sqrt(sum(diff)/length(diff));
end


function v = ParamsToVector(params)
    v = [params.max_val, params.steepness, params.midpoint_shifts];
end

function params = VectorToParams(v, shiftnum)
    params.max_val = v(1);
    params.steepness = v(2);
    params.midpoint = v(2+shiftnum);
end

function values = ComputeModel(params, irrads)

values = params.max_val./( 1+exp(-params.steepness*(irrads-params.midpoint)) );

end