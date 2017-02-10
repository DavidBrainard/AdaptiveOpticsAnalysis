
clear;
close all force;


controlBaseDir = fullfile(pwd,'control'); %uigetdir(pwd);
stimBaseDir = fullfile(pwd,'stim');

controlDataNames = read_folder_contents(controlBaseDir ,'mat');
stimDataNames = read_folder_contents(stimBaseDir ,'mat');

wavelet = 'gaus3';

stimd=0;
stimless=0;
x45=[];
x45ind=[];
x92=[];
x92ind =[];
x04=[];
x04ind =[];

controlcoeffs=[];
controllabels=[];
stimdcoeffs=[];
stimdlabels=[];

mintime = 10;
maxtime = 201;

cutoff=0;
for j=1:length(controlDataNames)

    load(fullfile(controlBaseDir, controlDataNames{j}));
    
    % Remove the empty cells
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
    control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );
    
    for i=1:length(norm_control_cell_reflectance)
        
        times = control_cell_times{i};
        signal = norm_control_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);        

        interptimes = mintime:maxtime;
        interpsignal = interp1(times,signal,interptimes,'pchip');
        
        D5 = wden(interpsignal,'sqtwolog','s','one',5,'bior3.5');

        before = ( D5( interptimes>10 & interptimes<=66 ) );
        during = ( D5( interptimes>66 & interptimes<=99 ));
        after = ( D5( interptimes>99 ) );
        
%         deriv5 = diff(D5);
        % Determine distance from max response to stim?
        [~,maxind]=max(abs(D5));
        
        controlcoeffs = [controlcoeffs; std(during) maxind-66 ];
        
        controllabels = [controllabels; {'control'}];
        
%         figure(1); title('Control cones'); 
%         plot(D5); hold on; 
%         plot(interpsignal);
%         axis([0 250 -10 10]); hold off;
%         
%         Nx=length(D5);
%         stationarylen=60;
%         overlap=.5;
%         wind = hamming(stationarylen);

%         figure(2);
%         spectrogram(D5,wind,overlap*stationarylen);

        
%         pause(1);
        controlreconst(i,:) = D5;
%         end

    end
    

end

stimless=[];
stimd=[];
for j=1:length(stimDataNames)
    
    load(fullfile(stimBaseDir, stimDataNames{j}));
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );
    
    
    for i=1:length(norm_stim_cell_reflectance)
        times  = stim_cell_times{i};
        signal = norm_stim_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);
        
        interptimes = mintime:maxtime;
        interpsignal = interp1(times,signal,interptimes,'pchip');

        D5 = wden(interpsignal,'sqtwolog','s','one',5,'bior3.5');

        deriv5 = diff(D5);
        
        before = ( D5( interptimes>10 & interptimes<=66 ) );
        during = ( D5( interptimes>66 & interptimes<=99 ));
        after = ( D5( interptimes>99 ) );

        [~,maxind]=max(abs(D5));

        stimd = [stimd i];

        stimdcoeffs = [stimdcoeffs; std(during) maxind-66 ];
        stimdlabels = [stimdlabels; {'stimulus'}];

        figure(3);  title(['Stim cones: ' num2str(max(D5)-min(D5)) ]); %hold on; 
        plot(D5); hold on; plot(interpsignal); 
        axis([0 250 -10 10]); hold off;

        Nx=length(D5);
        stationarylen=40;
        overlap=.5;
        wind = hamming(stationarylen);

        figure(4);
        spectrogram(interpsignal,wind,overlap*stationarylen);
        caxis([-80 30]);
        pause(1);
    end
    

    alllabels = [controllabels; stimdlabels];
    allcoeffs = [controlcoeffs; stimdcoeffs];

end    

SVMModel = fitcsvm(allcoeffs,alllabels,'KernelFunction','linear',...
                                       'KernelScale','auto',...
                                       'Standardize',true,...
                                       'BoxConstraint',10,'OutlierFraction',0.1);
sv = SVMModel.SupportVectors;

d = 0.02; % Step size of the grid
[x1Grid,x2Grid] = meshgrid( min(allcoeffs(:,1)):d:max(allcoeffs(:,1)),...
                            min(allcoeffs(:,2)):d:max(allcoeffs(:,2)) );
xGrid = [x1Grid(:),x2Grid(:)];        % The grid
[~,scores1] = predict(SVMModel,xGrid); % The scores

figure;
h(1:2) = gscatter(allcoeffs(:,1), allcoeffs(:,2), alllabels);
hold on
h(3) = plot(allcoeffs(SVMModel.IsSupportVector,1),...
            allcoeffs(SVMModel.IsSupportVector,2),'ko','MarkerSize',10);
    % Support vectors
contour(x1Grid, x2Grid, reshape(scores1(:,2), size(x1Grid)),[0 0],'k');
    % Decision boundary
title('Scatter Diagram with the Decision Boundary')
legend({SVMModel.ClassNames{1},SVMModel.ClassNames{2},'Support Vectors'},'Location','Best');
hold off


% figure
% gscatter(allcoeffs(:,1),allcoeffs(:,2),alllabels)
% hold on
% plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
% legend('control','stimulus','Support Vector')
% hold off

% CVSVMModel = crossval(SVMModel);
% classLoss = kfoldLoss(CVSVMModel)


%     cind = randperm(length(norm_control_cell_reflectance), 11);
%     aind = randperm(length(stimless), 11);
%     sind = randperm(length(stimd), 11);
%     
%     x45ind = [x45ind cind]; 
%     x92ind = [x92ind aind]; 
%     x04ind = [x04ind sind];
%     x45 = [x45; control_cell_times(cind) norm_control_cell_reflectance(cind)];    
%     x92 = [x92; stim_cell_times(stimless(aind)) norm_stim_cell_reflectance(stimless(aind))];
%     x04 = [x04; stim_cell_times(stimd(sind)) norm_stim_cell_reflectance(stimd(sind))];
%     imagesc(avg_spect);
