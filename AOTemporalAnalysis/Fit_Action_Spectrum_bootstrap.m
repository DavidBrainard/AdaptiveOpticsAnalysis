clear;
clc
close all;

rootDir = uigetdir(pwd);

fPaths = read_folder_contents(rootDir,'mat');

alldata = {};

% Load and sort all of the mat files we can find in this folder.
for i=1:length(fPaths)

    fname=fPaths{i};
    
    fpieces = textscan(fname,'%s','Delimiter','_','MultipleDelimsAsOne',1);
    fpieces = fpieces{1};        
    
    id = [fpieces{1} '_' fpieces{2}];
    stim_wavelen = [ 'w' fpieces{3}(1:3) ];
    stim_intensity = ['i' strrep(fpieces{4}(1:end-2),'.','p') ];
    stim_time = ['t' fpieces{5}];
    
    load(fPaths{i});
    
    % Build structs as needed to hold the data.
    idInd = -1;
    for j=1:length(alldata)
       if strcmp(alldata{j}.ID, id)
          idInd = j;
          break;
       end
    end
    
    if idInd == -1
        newstruct = struct();
        newstruct.ID = id;
        newstruct.(stim_wavelen).(stim_intensity)=all_amps;

        alldata = [alldata; {newstruct}];
    else        
        alldata{idInd}.(stim_wavelen).(stim_intensity)=all_amps;        
    end
end

fit_type = 'sigmoid';

%% Load CIE data
ciefunc = dlmread( fullfile(pwd,'linCIE2008v2e_5.csv') );
% conefunda = dlmread( fullfile(pwd,'linss2_10e_5.csv') );



%% Create the fits.
legendIDs = cell(length(alldata),1);
% luminosity = cell(1,length(wavelengths));

for i=1:length(alldata)
        
    subdata = alldata{i};
    
    fields = fieldnames(subdata);
    wavelengths = fields(cellfun(@(s) strcmp(s(1),'w'), fields)); % If it contains a 'w' it's a wavelength.
    
    irradiances = cell(1,length(wavelengths));    
    datavals = cell(1,length(wavelengths));
    
    
    
    for w=1:length(wavelengths)

        powers = fieldnames(subdata.(wavelengths{w}));
        
        % Clean up the string, and make the irradiances meaningful.
        irradiances(w) = {cellfun(@(s) str2double(strrep( strrep(s,'p','.') ,'i','')), powers)};
        
        [irradiances{w}, sortind] = sort(irradiances{w});
        
        for p=1:length(powers)            
            datavals{w} = [datavals{w}; { subdata.(wavelengths{w}).(powers{sortind(p)}) } ];
        end        
    end
    % Remove the garbage from the front of the wavelengths
    wavelengths = cellfun(@(s) str2double(strrep(s,'w','')), wavelengths);
    %%        
    parfor f=1:1000
        rng('shuffle');
        
        ampvals = cell(length(datavals),1);
        % Construct a "normal" dataval matrix.
        for d=1:length(datavals)
            inds = randi(1000, length(datavals{d}),1);
            
            somedata = cell2mat(datavals{d});
            somedata = somedata(:,inds); % When we do this, only the diagonal indices are valid.
            somedata = somedata(logical(eye(length(inds)))); % Use an identity matrix to extract the values we want.
            
            ampvals{d} = [somedata];
        end
        
        allamps{f} = ampvals;
        [rel_shifts{i,f}, datafit{i,f} ] = action_spect_fit(wavelengths, irradiances, ampvals', ciefunc, fit_type);

    end
    %%
   
    sumamps = cell(1,5);
    sumfits = cell(1,5);
    labels = cell(1,10);
    irrads = -10:.1:10;

    figure(1);clf; hold on;
    for w=1:length(wavelengths)
        sumamps{w} = zeros(length(allamps{1}{w}), 1000);
        for f=1:1000
            thisamp = allamps{f}{w};
            sumamps{w}(:,f) = thisamp;
            
        end
        meanamps{w} = mean(sumamps{w},2);
%         errorbar(log10(irradiances{w}), meanamps{w},1.96*std(sumamps{w},[],2)*sqrt(1+1/1000) );
        errorbar(log10(irradiances{w}), meanamps{w}, prctile(sumamps{w},5,2)-meanamps{w}, prctile(sumamps{w},95,2)-meanamps{w});
        labels{w} = num2str(wavelengths(w));
        
    end
    
    [meanshifts, datafit] = action_spect_fit(wavelengths, irradiances, meanamps, ciefunc, fit_type);
    figure(1);
    for w=1:length(wavelengths)
        plot(irrads, datafit{w},'-.');
        labels{w+length(wavelengths)} = [num2str(wavelengths(w)) 'Fit'];
    end
    axis([0 5 0 2.5]); legend(labels,'Location','NorthWest');
    title(subdata.ID);
    drawnow;
    saveas(gcf,[subdata.ID '_fit_werr.svg']);
    
    
end


% For getting the irradiance response as a function of wavelength
% for w=1:length(wavelengths)
%     figure(wavelengths(w));
%     legend(legendIDs);
%     saveas(gcf, fullfile(pathname, [num2str(wavelengths(w)) 'irr_response.svg']));
% end

figure; hold on;
for i=1:length(alldata)
    legendIDs{i} = alldata{i}.ID(4:end);
    subshifts = cell2mat(rel_shifts(i,:)');
    meanshifts(i,:) = mean( log10(subshifts) );
%     stdshifts(i,:) = std(subshifts);

    errorbar(wavelengths, mean( log10(subshifts) ), prctile(log10(subshifts),5)-mean( log10(subshifts) ),...
                                                    prctile(log10(subshifts),95)-mean( log10(subshifts) ),...
                                                    'MarkerSize',15);

end
plot(ciefunc(:,1), log10(ciefunc(:,2)),'k')
% plot(conefunda(:,1), log10(conefunda(:,2)), 'r');
% plot(conefunda(:,1), log10(conefunda(:,3)), 'g');
% set(gca,'yscale','log')
legend(legendIDs);
axis([450 700 -3 1]);
hold off;

figure; errorbar(wavelengths, mean( meanshifts ), std(meanshifts)); hold on;
% title('Mean Wavelength vs Normalized Irradiance'); xlabel('Wavelength'); ylabel('Action (relative)');
plot(ciefunc(:,1), log10(ciefunc(:,2)),'k')
% % set(gca,'yscale','log')
axis([450 700 -3 1]);
% hold off;
