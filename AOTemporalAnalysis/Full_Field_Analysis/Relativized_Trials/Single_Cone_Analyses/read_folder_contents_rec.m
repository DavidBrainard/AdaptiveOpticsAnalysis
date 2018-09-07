function [ file_paths ] = read_folder_contents_rec( basedir, extension )

file_paths={};

[file_list, isdir, numfiles] = read_folder_contents(basedir);

for i=1:numfiles   
    if isdir{i}
        file_paths = [file_paths; read_folder_contents_rec(fullfile(basedir,file_list{i}), extension)];
    elseif strcmp(file_list{i}(end-2:end),extension)
        file_paths = [file_paths; {fullfile(basedir,file_list{i})}];
    end    
end

end

