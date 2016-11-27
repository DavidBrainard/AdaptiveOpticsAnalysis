
clear;
close all force;


baseDir = pwd; %uigetdir(pwd);
profileDataNames = read_folder_contents(baseDir ,'mat');

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

maxtime = 175;

for j=1:length(profileDataNames)

    load(fullfile(baseDir, profileDataNames{j}));
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
    control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );
    
    all_spect = zeros(length(norm_stim_cell_reflectance), 256);
    cutoff = 0;
    for i=1:length(norm_control_cell_reflectance)
        
        times = control_cell_times{i};
        signal = norm_control_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);        

        interptimes = 1:249;
        interpsignal = interp1(times,signal,interptimes,'pchip');
        
%         figure(10);
%         plot(times, signal,1:249,interpsignal); %hold on;
%         plot(66:98,ones(33,1)*max(signal), 'r*'); hold off;
        
        %Continuous 8,16,32,64
        D5 = cwt(interpsignal,16,wavelet);
        
%         if length(D5) >= 190

        %Discrete
%         [C, L] = wavedec(signal,5,wavelet);
%         [cD1,cD2,cD3,cD4,cD5] = detcoef(C,L,1:5);
% 
%         D4 = wrcoef('d',C,L,wavelet,4);
%         D5 = wrcoef('d',C,L,wavelet,5);

%         subplot(5,1,1); plot(D1);
%         subplot(5,1,2); plot(D2);
%         subplot(5,1,3); plot(D3);
%         subplot(5,1,4); plot(D4);
%         subplot(5,1,5);
%          plot(D5(1,1:maxtime));
        
%          figure(100);
%          plot([signal' D5'])
         
        D5 = D5(1:maxtime);
        interptimes = interptimes(1:maxtime);
        before = ( D5( interptimes<=66 ) );
        during = ( D5( interptimes>66 & interptimes<=99 ));
        after = ( D5( interptimes>99 ) );
         
        controlcoeffs = [controlcoeffs; std(during) std(after) ];
        
        controllabels = [controllabels; {'control'}];
        
        figure(1); title('Control cones'); %hold on; 
        plot(D5); axis([0 250 -10 10])
        

        controlreconst(i,:) = D5;
%         end

    end
    
    simplethreshold = 3*std(controlreconst);
    
    stimless=[];
    stimd=[];
    for i=1:length(norm_stim_cell_reflectance)
        times  = stim_cell_times{i};
        signal = norm_stim_cell_reflectance{i}(times>=cutoff);
        times = times(times>=cutoff);
        
        interptimes = 1:249;
        interpsignal = interp1(times,signal,interptimes,'pchip');
        
%         figure(10);
%         plot(times, signal); hold on;
%         plot(66:98,ones(33,1)*max(signal), 'r*'); hold off;

        %Continuous
        D5 = cwt(signal,16,wavelet);
        
%         if length(D5) >= 190
%             

            %Discrete
%             [C, L] = wavedec(signal,5,wavelet);
%             [cD1,cD2,cD3,cD4,cD5] = detcoef(C,L,1:5);
            
%             D4 = wrcoef('d',C,L,wavelet,4);
%             D5 = wrcoef('d',C,L,wavelet,5);

        D5 = D5(1:maxtime);
        interptimes = interptimes(1:maxtime);
        before = ( D5( interptimes<=66 ) );
        during = ( D5( interptimes>66 & interptimes<=99 ));
        after = ( D5( interptimes>99 ) );

            if all(D5(1:maxtime)<simplethreshold & D5(1:maxtime)>-simplethreshold)
                stimless = [stimless i];

%                 figure(3);  title('StimLESS cones'); hold on;
%                 plot(D5(1:maxtime)); 
%                 axis([0 250 -10 10])
            else            
                stimd = [stimd i];

                stimdcoeffs = [stimdcoeffs; std(during) std(after) ];
                stimdlabels = [stimdlabels; {'stimulus'}];
                
%                 figure(2);  title('Stim cones'); hold on;
%                 plot(D5(1:maxtime)); 
%                 axis([0 250 -10 10])
                
            end
%         end
    end
    
    
    
    alllabels=[controllabels; stimdlabels];
    allcoeffs=[controlcoeffs; stimdcoeffs];

end    
    %%
SVMModel = fitcsvm(allcoeffs,alllabels,'KernelFunction','linear',...
                                       'KernelScale','auto',...
                                       'Standardize',true,...
                                       'BoxConstraint',10,'OutlierFraction',0.1);
sv = SVMModel.SupportVectors;

d = 0.02; % Step size of the grid
[x1Grid,x2Grid] = meshgrid(min(allcoeffs(:,1)):d:max(allcoeffs(:,1)),...
    min(allcoeffs(:,2)):d:max(allcoeffs(:,2)));
xGrid = [x1Grid(:),x2Grid(:)];        % The grid
[~,scores1] = predict(SVMModel,xGrid); % The scores

figure;
h(1:2) = gscatter(allcoeffs(:,1),allcoeffs(:,2),alllabels);
hold on
h(3) = plot(allcoeffs(SVMModel.IsSupportVector,1),...
    allcoeffs(SVMModel.IsSupportVector,2),'ko','MarkerSize',10);
    % Support vectors
contour(x1Grid,x2Grid,reshape(scores1(:,2),size(x1Grid)),[0 0],'k');
    % Decision boundary
title('Scatter Diagram with the Decision Boundary')
legend({'-1','1','Support Vectors'},'Location','Best');
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
