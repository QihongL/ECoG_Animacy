function [ nSubjs, allSubjIDs ] = getSubjInfoFromResults( results )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
nSubjs = length(results); % assume all timePts have the same # of subj
allSubjIDs = nan(nSubjs,1);
for i = 1 : nSubjs
    allSubjIDs(i) = results{i}.subjID;
end

end

