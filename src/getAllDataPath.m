function [fullpath] =  getAllDataPath(dir, allDirNames, windowSize)

% get the directories for all starting time points 
fullpath = cell(length(allDirNames),1);
for i = 1 : length(allDirNames)
    fullpath{i} = sprintf(fullfile(dir, allDirNames(i).name, 'WindowSize', windowSize, '/'));
end

end