% []=Random_Trial_Splitter()
% This script takes in a folder and randomly splits the total number of
% files in the folder in to as many splits as are desired.
%
% Created by Robert F Cooper 2019-02-25
%

clear;

rootDir = uigetdir(pwd, 'Select the directory containing the profiles you wish to split.');

profileDataNames = read_folder_contents(rootDir,'mat');

NUM_SPLITS = 2;

numpersplit = floor(size(profileDataNames, 1) / NUM_SPLITS);


randFileSort = randperm(size(profileDataNames,1));
outdirs=cell(NUM_SPLITS,1);

for n=1:NUM_SPLITS
    splitindstart = (n-1)*numpersplit +1;
    splitindend = n*numpersplit;
    outdirs{n} = fullfile(rootDir, ['Split ' num2str(n)]);
    % Need to capture remainder
    if mkdir(outdirs{n})
        for f=splitindstart:splitindend
            copyfile( fullfile(rootDir, profileDataNames{randFileSort(f)}), fullfile(outdirs{n}, profileDataNames{randFileSort(f)}) );
        end
    end
end

if splitindend<size(profileDataNames, 1)
    f=splitindend;
    n=1;
    while f<size(profileDataNames, 1)        
        copyfile( fullfile(rootDir, profileDataNames{randFileSort(f)}), ...
                  fullfile(outdirs{mod(n,NUM_SPLITS)}, profileDataNames{randFileSort(f)}) );
        f=f+1;
        n=n+1;
    end
end
f