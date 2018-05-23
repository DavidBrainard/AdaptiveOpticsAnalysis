

clear;
close all;

rng('shuffle');

siggen = dsp.ColoredNoise('SamplesPerFrame',165,'NumChannels',10000);

% Make a bunch of fake null signals
std_dev_sub = step(siggen);

mean_sub = std_dev_sub(:,5001:end)./10;
std_dev_sub = std_dev_sub(:,1:5000)./10;
timeBase = ((1:165)/16.6)';

fitAmp = nan(size(std_dev_sub,1),1);
fitMean = nan(size(std_dev_sub,1),1);
fitAngle = nan(size(std_dev_sub,1),1);


parfor j=1:size(std_dev_sub,2)

    j
    std_dev_sig = std_dev_sub(:,j);
%     padding_amt = ceil((2^(nextpow2(length(std_dev_sig)))-length(std_dev_sig)) /2);
%     padded_stddev_sig = padarray(std_dev_sig, [0  padding_amt],'symmetric', 'both');
%     padded_stddev_sig=wavelet_denoise( padded_stddev_sig );
%     filt_stddev_sig = padded_stddev_sig(padding_amt+1:end-padding_amt);
%     
%     figure(200); plot(timeBase,std_dev_sig,timeBase,filt_stddev_sig);
    
    mean_sig = mean_sub(:,j);
%     padding_amt = ceil((2^(nextpow2(length(mean_sig)))-length(mean_sig)) /2);
%     padded_mean_sig = padarray(mean_sig, [0  padding_amt],'symmetric', 'both');
%     padded_mean_sig=wavelet_denoise( padded_mean_sig );
%     filt_mean_sig = padded_mean_sig(padding_amt+1:end-padding_amt);
%             
%     figure(201); plot(timeBase,mean_sig,timeBase,filt_mean_sig);    
%     drawnow;
    
    if ~all( isnan(std_dev_sig) )
%         figure(2);clf; plot(timeBase,std_dev_sig);
        fitData = modelFit_beta(timeBase, std_dev_sig, []);
        fitAmp(j) = fitData.amplitude;
%         pause(1);

%         figure(2);clf; plot(timeBase,mean_sig);
        fitData = modelFit_beta(timeBase, mean_sig, [] );
%         pause(1);
        fitMean(j) = fitData.amplitude;        

    end
end


%%
figure(1); clf; hold on;
for j=1:size(std_dev_sub,2)
    if ~isnan(fitAmp(j))
        plot( fitAmp(j),fitMean(j),'k.');        
    end
end
ylabel('Mean response amplitude');
xlabel('Reflectance response amplitude');
title('Mean vs std reflectance response amplitude')
hold off;