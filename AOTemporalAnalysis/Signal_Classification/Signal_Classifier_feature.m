
clear;
close all force;


controlBaseDir = fullfile(pwd,'control'); %uigetdir(pwd);
stimBaseDir = fullfile(pwd,'stim');

controlDataNames = read_folder_contents(controlBaseDir ,'mat');
stimDataNames = read_folder_contents(stimBaseDir ,'mat');

wavelet = 'gaus3';

prefilt = designfilt('lowpassiir', 'FilterOrder', 5, 'PassbandFrequency', .06, 'PassbandRipple', .2);
% prefilt = designfilt('bandpassiir', 'FilterOrder', 6, 'PassbandFrequency1', .02, 'PassbandFrequency2', .1, 'PassbandRipple', .2);

stimd=0;
stimless=0;
x45=[];
x45ind=[];
x92=[];
x92ind =[];
x04=[];
x04ind =[];




allsigs=[];

clipind = 11; % Need to clip ends due to artifacts after filtering
stimind = 66;
maxind = 241;


affinitycurve = normpdf(0:maxind-stimind-1,16,33);
affinitycurve = affinitycurve./max(affinitycurve);

cutoff=0;
for j=1:length(controlDataNames)

    load(fullfile(controlBaseDir, controlDataNames{j}));
    
    % Remove the empty cells
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
    control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );
    
    controlcoeffs=[];
    controllabels=[];
    
    alllabels{j} = [];
    allcoeffs{j} = [];
    
    for i=1:length(norm_control_cell_reflectance)
        
        times = control_cell_times{i};
        signal = norm_control_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);        

        interptimes = 0:maxind;
        interpsignal = interp1(times,signal,interptimes,'pchip');
        
        interpsignal = interpsignal(clipind:maxind);

        filtinterpsignal = filtfilt(prefilt, interpsignal);
        
        % Remove the remaining signal before the stimulus delivery
        interpsignal = interpsignal((stimind-clipind)+1:end);
        filtinterpsignal = filtinterpsignal((stimind-clipind)+1:end-1);
        
        allsigs=[allsigs; filtinterpsignal];
                
        during = filtinterpsignal( 1:33 );
        after = filtinterpsignal( 34:end-clipind );

        % Determine distance from max response to stim?
        [maxrespval,maxrespind]=max( abs( diff([during after])) );
        
        stim_affinity = affinitycurve(maxrespind)*maxrespval;
        
        derivduring = diff(during);
        derivafter = diff(after);
        
        [peakvals, peaks] = findpeaks([derivduring derivafter]);
        [troughvals, troughs] = findpeaks(-[derivduring derivafter]);
        
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
        
        SWC = swt(interpsignal,4,'db4');
        
        % Features from Subasi et al 2007
        meanabscoeff = mean(abs(SWC'));
        meanpowercoeff = sum(SWC'.^2);
        stddevcoeff = std(SWC');
        coeffratio = [meanabscoeff(1)/meanabscoeff(2) meanabscoeff(2)/meanabscoeff(3) ...
                      meanabscoeff(3)/meanabscoeff(4) meanabscoeff(4)/meanabscoeff(5)];
        
        % Put together the feature lists
        controlcoeffs = [controlcoeffs; stim_affinity std(pk_pk) max(derivduring)-min(derivduring) meanabscoeff(3:5) meanpowercoeff(3:5) stddevcoeff(3:5) coeffratio(2:4)];        
        controllabels = [controllabels; {'control'}];
        
        
%         figure(1); title('Control cones'); hold on;
% %         plot(D5);
% %         hold on; 
%         plot( diff([during after]) );
%         axis([0 249 -1 1]);
% %         axis([0 250 -10 10]);
%         hold off;
        
        controlreconst(i,:) = filtinterpsignal;


    end
    alllabels{j} = [alllabels{j}; controllabels];
    allcoeffs{j} = [allcoeffs{j}; controlcoeffs];
end




stimless=[];
stimd=[];
for j=1:length(controlDataNames)
    
    load(fullfile(stimBaseDir, stimDataNames{j}));
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );
    stimdlabels=[];
    stimdcoeffs=[];
    
    for i=1:length(norm_stim_cell_reflectance)
        times  = stim_cell_times{i};
        signal = norm_stim_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);
        
        interptimes = 0:maxind;
        interpsignal = interp1(times,signal,interptimes,'pchip');

        interpsignal = interpsignal(clipind:maxind);

        filtinterpsignal = filtfilt(prefilt, interpsignal);
        
        % Remove the remaining signal before the stimulus delivery
        interpsignal = interpsignal((stimind-clipind)+1:end);
        filtinterpsignal = filtinterpsignal((stimind-clipind)+1:end-1);
        
        allsigs=[allsigs; filtinterpsignal];
                
        during = filtinterpsignal( 1:33 );
        after  = filtinterpsignal( 34:end-clipind );

        % Determine distance from max response to stim?
        [maxrespval,maxrespind]=max( abs( diff([during after])) );
        
        stim_affinity = affinitycurve(maxrespind)*maxrespval;
        
        derivduring = diff(during);
        derivafter = diff(after);
        
        [peakvals, peaks] = findpeaks([derivduring derivafter]);
        [troughvals, troughs] = findpeaks(-[derivduring derivafter]);
        
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
        
        SWC = swt(interpsignal,4,'db4');
        
        % Features from Subasi et al 2007
        meanabscoeff = mean(abs(SWC'));
        meanpowercoeff = sum(SWC'.^2);
        stddevcoeff = std(SWC');
        coeffratio = [meanabscoeff(1)/meanabscoeff(2) meanabscoeff(2)/meanabscoeff(3) ...
                      meanabscoeff(3)/meanabscoeff(4) meanabscoeff(4)/meanabscoeff(5)];
        
        stimd = [stimd i];

        if( stim_affinity > .05)
            stimdcoeffs = [stimdcoeffs; stim_affinity std(pk_pk) max(derivduring)-min(derivduring) meanabscoeff(3:5) meanpowercoeff(3:5) stddevcoeff(3:5) coeffratio(2:4)];
            stimdlabels = [stimdlabels; {'stimulus'}];
        end

%         figure(3); title(['Stim cones']); hold on; 
% % %         plot(D5); 
% % %         hold on; 
%         plot( diff([during after]) );
% %         axis([0 length([during after]) -1 1]);
% %         %plot(interpsignal); 
% % %         axis([0 250 -10 10]);
%         axis([0 250 -1 1]);
%         hold off;
        
        
    end
    alllabels{j} = [alllabels{j}; stimdlabels];
    allcoeffs{j} = [allcoeffs{j}; stimdcoeffs];
end



% allcoeffs = allcoeffs( [1:3769 3771:end],: );
% alllabels = alllabels( [1:3769 3771:end],: );


% explained


% projected = allcoeffs*pcacoeffs;

% figure;
% gscatter(projected(:,1), projected(:,3), alllabels);


% Pick a random set to fit from
dataSetInds = randperm(length(allcoeffs));
dataSetInds(1)
SVMModel = fitcsvm(allcoeffs{dataSetInds(1)},alllabels{dataSetInds(1)},'KernelFunction','linear',...
                                             'KernelScale','auto',...
                                             'Standardize',true,...
                                             'BoxConstraint',10,'CrossVal','on','KFold',10,'OutlierFraction',0.1);

kfoldPercentModelLoss = 100*kfoldLoss(SVMModel)

% Aggregate all of the other data
validationcoeff=[];
validationlabels=[];
for o=2:length(dataSetInds)
   
    validationcoeff = [validationcoeff; allcoeffs{dataSetInds(o)}];
    validationlabels = [validationlabels; alllabels{dataSetInds(o)}];
    
end

[pcacoeffs, pcascore, latent, ~, explained] = pca( validationcoeff );

% Estimate the classification error.
classificationPercentloss = 100*loss(SVMModel.Trained{1},validationcoeff,validationlabels)

