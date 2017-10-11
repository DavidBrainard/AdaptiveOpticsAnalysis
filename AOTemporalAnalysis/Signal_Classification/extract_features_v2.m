function [ feature_coeffs, labels, varNames ] = extract_features_v2( cell_times, norm_cell_reflectance, critical_period, vid_type )

    

    % Remove empty cells
%     norm_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance)  );
%     cell_times    = cell_times( ~cellfun(@isempty,cell_times) );
    
    naners = ~cellfun(@any, cellfun(@isnan, norm_cell_reflectance, 'UniformOutput',false));
    norm_cell_reflectance = norm_cell_reflectance( naners );
    cell_times            = cell_times( naners );

    stimulus_length = (critical_period(2)-critical_period(1));
    wavelet = 'gaus4';
    scales = [10; 25; 40];
    sampling_period = 1/17;
    
    
    varNames = {'Affinity_Scale_10','During_Pwr_Scale_10','After_Pwr_Scale_10',...
                'Affinity_Scale_25','During_Pwr_Scale_25','After_Pwr_Scale_25',...
                'Affinity_Scale_40','During_Pwr_Scale_40','After_Pwr_Scale_40','Training_Labels'};
    
    feature_coeffs=zeros(length(norm_cell_reflectance),length(varNames)-1);
    labels=cell(length(norm_cell_reflectance),1);
        
    parfor i=1:length(norm_cell_reflectance)
        
        times = cell_times{i};
        signal = norm_cell_reflectance{i};
        
        if length(signal) >= 160 % The signal has to have enough data that we can reliably analyze it.
            

            interptimes = min(times):max(times);
%             datadiff = setdiff(times,interptimes);
            interpsignal = interp1(times,signal,interptimes,'linear');

            interptimes = interptimes(~isnan(interpsignal));
            interpsignal = interpsignal(~isnan(interpsignal));
            
            [coefs, freq] = cwt(interpsignal,scales,wavelet,sampling_period);            
%             conofinf(wavelet,1:128,length(coefs), 128,'plot');
            
            normalization = repmat(scales, 1, size(coefs,2));
            
            normcoefs = (coefs.^2)./normalization; % Per Liu et al 2007, "Rectification of the Bias in the Wavelet Power Spectrum"
            
            figure(1); title([vid_type ' cones']);
            plot(times,signal);  plot(interptimes, interpsignal);
            axis([0 max(times) -4 4]);
           
            figure(2); plot(normcoefs'); axis([0 size(normcoefs,2) -5 5]);
            legend( num2str(scales) )
            
            before = times>10 & times<=critical_period(1); % Avoid edge effects
            during = times>critical_period(1) & times<=critical_period(2);
            after = times>critical_period(2) & times<=max(times)-10;
            
            affinitycurve = normpdf(0:max(times)-1, critical_period(1)+16, 33);
            affinitycurve = affinitycurve./max(affinitycurve);
            
            avg_pwr=zeros(1,length(scales)*3);
            for s=1:length(scales)
                
                [maxrespval,maxrespind]=max( normcoefs(s, :) );
                stim_affinity = affinitycurve(maxrespind)*maxrespval;
                
                avg_pwr( (3*(s-1))+1:(3*(s-1))+3 ) = [stim_affinity mean(normcoefs(s, during)) mean(normcoefs(s, after))]; 
            end

            feature_coeffs(i, :) = avg_pwr;
            labels{i} = vid_type;
        else
            feature_coeffs(i,:) = nan(1,length(varNames)-1);
            labels{i} = NaN;
        end
    end
    


end

