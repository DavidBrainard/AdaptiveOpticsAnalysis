function [ rel_shifts, fits ] = action_spect_fit( wavelengths, irradiances, datavals, ciefunc, type )
% [ rel_shifts, fits ] = action_spect_fit( wavelengths, irradiances, datavals, ciefunc, type )
%
% This function is responsible for fitting the reflectance response vs irradiance functions for each
% wavelength from a single subject, and determining their relative offset. The fit to all wavelengths
% is performed simultaneously. The only fit parameter allowed to vary on a per-wavelength basis is
% the relative X-axis (irradiance) shift of each curve. Everything else (amplitude, shape, etc)
% is allowed to vary across all wavelengths' curves simultaneously.
%
% @params:
%    wavelengths: A 1xN array of wavelengths (in nm) that we'll be fitting. Each wavelength in this
%                 array corresponds to a single cell in parameters "irradiances" and "datavals".
%
%    irradiances: A 1xN cell array where each cell corresponds to an index/wavelength
%                 in the "wavelengths" parameter. Each cell the contains the irradiances used
%                 for the corresponding wavelength array value.
%
%    datavals: A 1xN cell array where each cell corresponds to an index/wavelength
%              in the "wavelengths" parameter. Each cell the contains the reflectance response
%              of the corresponding wavelength array value and irradiance cell.
%
%    ciefunc: The CIE function we're attempting to fit to.
%
%    type: The type of fit to use for the reflectance response vs irradiance. Options are:
%         "Sigmoid", "Linear", or "Exponential"
%
% @outputs:
%
%    rel_shifts: A 1xN array of relative shifts of each wavelength.
%
%    fits: A 1xN cell array with the fits of each wavelength's reflectance reponse vs irradiance
%         function
%
% Created by Robert F Cooper 01-24-2017
%
% The analyses performed in this script are from:
% Cooper RF, Tuten WS, Dubra A, Brainard BH, Morgan JIW. 
% "Non-invasive assessment of human cone photoreceptor function." 
% Biomed Opt Express 8(11): 5098-5112 and are
% encompassed in Figures 6/7, Equation 4.

normalize = false;

% Set up the structure
fitParams.max_or_intercept = [];
fitParams.irr_shift=[];
fitParams.slope=0;

max_irr = 10^5;
log_irr_range = -10:.1:10;


% Rescale amplitudes so that we can fit them all the same- each column is
% a different wavelength
handle = figure(99); clf;

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
            irrval(w) = {log10(irr)};
            fitParams.slope = fitParams.slope+1;
        case 'linear'
            irr=log10(irr);
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
    vub(w) = 3*max(irrval{w}); 
    vlb(w) = min(irrval{w});
    

    hold on;
    plot(irrval{w}, datavals{w},'.','MarkerSize',15);
    hold off;
end


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
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','interior-point'); % 'Diagnostics','off','Display','off'

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

%% Lock the parameters, then scale the curves up and down to fit the CIE

% Put the shifts back into a linear scale
linearshifts = 10.^finalParams(3:end);

rel_shifts = 1./( linearshifts ./ linearshifts(2) );

figure(100);
plot(ciefunc(:,1),ciefunc(:,2)); hold on;
plot(wavelengths,rel_shifts); title('Preshift');set(gca,'yscale','log');axis([450 700 10^-3 10^1]); hold off;

cierow = []; % Find the rows (wavelengths) in the CIE function that we care about
for w=1:length(wavelengths)    
    cierow = [cierow find(ciefunc(:,1) == wavelengths(w))];
end

cievals = ciefunc(cierow,2)';

vub = 4;
vlb = -3; 

v0 = 0;
% Shift our relative actions around to get the best fit    
finalShifts = fmincon(@(v)FitCIEErrorFunction(v, log10(cievals), log10(rel_shifts)),v0,[],[],[],[],vlb,vub,[],options);
% finalShifts = rel_shifts'\cievals';

rel_shifts= 10.^(log10(rel_shifts)+finalShifts);

figure(101);
plot(ciefunc(:,1),ciefunc(:,2)); hold on;
plot(wavelengths,rel_shifts); title('Postshift'); set(gca,'yscale','log');
axis([450 700 10^-3 10^1]);hold off;


wavecolors=['b' 'c', 'g', 'y','r'];
figure(handle);hold on;
for i=1:size(datavals,2)
    params = VectorToParams(finalParams, i);
    fits{i} = ComputeModel(params, log_irr_range, type);
    plot(log_irr_range, ComputeModel(params, log_irr_range, type),wavecolors(i));
    
end
legend(num2str(wavelengths));
hold off;
axis([0 max(cellfun(@max,irrval)) 0 max(cellfun(@max,datavals))]);


end

function f = FitCIEErrorFunction(v, cie, rel_action)

    preds = ComputeCIEModel(v,rel_action);
    
    diff = (cie-preds).^2; 

    f = sqrt(sum(diff)/length(diff));
end

function values = ComputeCIEModel(v, rel_action)

values = rel_action+v;

end

function f = FitModelErrorFunction(v, irradiances, values,type)

    diff=[];
    for i=1:size(values,2)

        params = VectorToParams(v,i);

        preds{i} = ComputeModel(params, irradiances{i},type);

        diff = [diff; (values{i}-preds{i}).^2];
        
%         plot(0:.1:11,ComputeModel(params,0:.1:11,type),'k-.', irradiances{i}, values{i}); hold on;
%         axis([0 5 0 3]);
    end        
%     drawnow;  hold off;
    
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
