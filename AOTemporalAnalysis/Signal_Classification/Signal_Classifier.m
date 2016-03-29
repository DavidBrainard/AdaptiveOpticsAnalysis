
clear;
close all force;


baseDir = uigetdir(pwd);
profileDataNames = read_folder_contents(baseDir,'mat');



for j=1:length(profileDataNames)

    load(fullfile(baseDir, profileDataNames{j}));
    
    % Remove the empty cells
    norm_stim_cell_reflectance = norm_stim_cell_reflectance( ~cellfun(@isempty,norm_stim_cell_reflectance) );
    stim_cell_times            = stim_cell_times(  ~cellfun(@isempty,stim_cell_times) );
    norm_control_cell_reflectance = norm_control_cell_reflectance( ~cellfun(@isempty,norm_control_cell_reflectance)  );
    control_cell_times            = control_cell_times( ~cellfun(@isempty,control_cell_times) );
    
    all_spect = zeros(length(norm_stim_cell_reflectance), 256);
    for i=1:length(norm_stim_cell_reflectance)
        
        times  = control_cell_times{i};
        signal = norm_control_cell_reflectance{i};
        
        
        % The profile itself
        figure(1);
        plot(times, signal); hold on;
        plot(66:98,ones(33,1)*max(signal), 'r*'); hold off;
%         
%         % Calculate the profile's log power spectrum
%         signal_len = length( signal );
%         diff_spec  = 256-signal_len;
%         
%         padded_sig = padarray(signal, [0 ceil(diff_spec/2)], 0, 'pre');
%         padded_sig = padarray(padded_sig, [0 floor(diff_spec/2)], 0, 'post');
         figure(2); spectrogram(signal,16,'yaxis');
        all_spect(i,times) = signal; %fftshift( log( abs(fft( padded_sig )).^2 ) );
    end
    imagesc(all_spect); colormap gray;
end