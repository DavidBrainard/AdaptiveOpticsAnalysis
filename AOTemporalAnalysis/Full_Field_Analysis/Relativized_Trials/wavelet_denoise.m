function [ denoised_signal ] = wavelet_denoise( signal )

wavelet = 'haar';
level=5;
sorh = 's';

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

end

