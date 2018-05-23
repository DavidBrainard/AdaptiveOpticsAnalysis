function [ feature_coeffs, labels ] = extract_features( cell_times, norm_cell_reflectance, critical_period, vid_type )

    numfeatures  =23;

    % Remove empty cells
    norm_cell_reflectance = norm_cell_reflectance( ~cellfun(@isempty,norm_cell_reflectance)  );
    cell_times    = cell_times( ~cellfun(@isempty,cell_times) );
    
    naners = ~cellfun(@any, cellfun(@isnan, norm_cell_reflectance, 'UniformOutput',false));
    norm_cell_reflectance = norm_cell_reflectance( naners );
    cell_times            = cell_times( naners );
    
    feature_coeffs=zeros(length(norm_cell_reflectance),numfeatures);
    labels=cell(length(norm_cell_reflectance),1);
    
    stimulus_length = (critical_period(2)-critical_period(1));
    wavelet = 'db4';
    level=5;
    sorh = 's';    % Specified soft or hard thresholding
    maxind = 220;
    
        
    parfor i=1:length(norm_cell_reflectance)
        
        times = cell_times{i};
        signal = norm_cell_reflectance{i};
        
        if length(signal) >= 160 % The signal has to have enough data that we can reliably denoise it.
%         poststim = find(times>critical_period(1));

            nextpowpad = (2^nextpow2(length(signal)))-length(signal);
            padsignal = padarray(signal, [0 nextpowpad],'symmetric','pre');

            % *** Generated using wavemenu input. ***
            % Denoise.
            %---------        
            [SWC] = swt(padsignal,level, wavelet);
            valTHR = wthrmngr('sw1ddenoLVL','heursure',SWC,'one')';
            maxTHR = max(abs(SWC(1:level,:)),[],2);        
            valTHR = min(valTHR,maxTHR);

            len = length(padsignal);

            thrSettings = repmat([1 len],[level 1]);
            thrSettings = [thrSettings valTHR];
            SWC_thresh = SWC;

            for k = 1:level
                thr_par = thrSettings(k,:);
                if ~isempty(thr_par)
                    NB_int = size(thr_par,1);
                    x      = [thr_par(:,1) ; thr_par(NB_int,2)];
                    x      = round(x);
                    x(x<1) = 1;
                    x(x>len) = len;
                    thr = thr_par(:,3);
                    for j = 1:NB_int
                        if j==1 , d_beg = 0; else d_beg = 1; end
                        j_beg = x(j)+d_beg;
                        j_end = x(j+1);
                        j_ind = (j_beg:j_end);
                        SWC_thresh(k,j_ind) = wthresh(SWC_thresh(k,j_ind),sorh,thr(j));
                    end
                end
            end
            % Reconstruct the denoised signal using ISWT.        
            denoised_signal = iswt(SWC_thresh, wavelet);
            %-------------------------------------------


            denoised_signal = denoised_signal((nextpowpad+1):end);
            
            interptimes = 0:maxind;
            poststim = find( interptimes>critical_period(1) );
            interpsignal = interp1(times,denoised_signal,interptimes,'pchip');

            interptimes = interptimes(~isnan(interpsignal));
            interpsignal = interpsignal(~isnan(interpsignal));
            
            % Remove the signal before the stimulus delivery
            interptimes = interptimes(poststim(1):end);
            interpsignal = interpsignal(poststim(1):end);
%             denoised_signal = denoised_signal(poststim(1)+nextpowpad:end);
%             denoised_signal = denoised_signal-denoised_signal(1); % zero-center the signal.
    %         times = times(poststim(1):end);


            during = denoised_signal( 1:(stimulus_length*2) );
    %         after = denoised_signal( stimulus_length+1:end );

            % Determine distance from max response to stim
            affinitycurve = normpdf(0:length(denoised_signal)-1, 16, 33);
            affinitycurve = affinitycurve./max(affinitycurve);
    %         % If there are things happening during the stimulus, then that's really high affinity
    %         affinitycurve(1:stimulus_length) = 1; 

            [maxrespval,maxrespind]=max( abs( during ));

            stim_affinity = affinitycurve(maxrespind)*maxrespval;

            denoised_derv_during = diff(during);
            denoised_derv = diff(denoised_signal);

            [peakvals, peaks] = findpeaks(denoised_derv);
            [troughvals, troughs] = findpeaks(-denoised_derv);
            numzeros=0;

            auc = sum(abs(denoised_derv_during));

            if isempty(peaks) 
               pk_pk = troughs; 
            elseif isempty(troughs)
               pk_pk = peaks;
            else        
                if peaks(1) < troughs(1)
                   listlen = length(peaks);
                else
                   listlen = length(troughs);
                end

                for p=1:listlen
                   if p > length(peaks) || p > length(troughs)
                      break; 
                   end

                    pk_pk(p) = peakvals(p)+troughvals(p);

                end
            end



    %         [c, l] = wavedec(interpsignal, level, wavelet);
    %         
    %         DWC = detcoef(c,l,1:level);
            DWC = modwt(signal,wavelet,4);

    %         meanabscoeff=zeros(1,length(DWC));
    %         meanpowercoeff=zeros(1,length(DWC));
    %         stddevcoeff=zeros(1,length(DWC));

    %         for p=1:size(DWC,1)
                % Features from Subasi et al 2007
    %             meanabscoeff(p) = mean(abs( DWC{p} ));
    %             meanpowercoeff(p) = sum( DWC{p}.^2);
    %             stddevcoeff(p) = std( DWC{p} );
    %         end

            meanabscoeff = mean(abs( DWC' ));
            meanpowercoeff = sum( DWC'.^2);
            stddevcoeff = std( DWC' ); 

            coeffratio = zeros(1,length(stddevcoeff)-1);

            for p=1:size(DWC,1)-1
                coeffratio(p) = meanabscoeff(p)/meanabscoeff(p+1);
            end


           figure(1); title([vid_type ' cones']); hold off;
           plot(times,signal); hold on; plot(times,denoised_signal); plot(interptimes, interpsignal);


    %         stim_affinity
    %         auc
    % length([stim_affinity auc max(denoised_derv_during)-min(denoised_derv_during)...
    %                                max(during)-min(during) meanabscoeff meanpowercoeff stddevcoeff coeffratio])
            % Put together the feature lists
            feature_coeffs(i,:) = [stim_affinity auc std(pk_pk)...
                                   max(during)-min(during) meanabscoeff meanpowercoeff stddevcoeff coeffratio];

            labels{i} = vid_type;
        else
            feature_coeffs(i,:) = nan(1,numfeatures);
            labels{i} = NaN;
        end
    end
end

