function [ filename, subjIDs  ] = getFileNames( datadir, DATA_TYPE )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% dynamically create the file names for each subjects 
listing = dir(strcat(datadir, 's*_', DATA_TYPE, '.mat'));
subjIDs = nan(length(listing),1);
for s = 1 : length(listing)
    % get file name corresponds to the X for each participants
    filename.data{s} = sprintf('%s\n', listing(s).name);
    % ... and their IDs
    subjIDs(s) = sscanf(filename.data{s},'%*c%f%*c',[1 Inf]);
end
% get the corresponding metadata
filename.metadata = strcat('metadata_', DATA_TYPE, '.mat');
end

