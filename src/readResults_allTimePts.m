function [ accuracies, subjIDs ] = readResults_allTimePts( results )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nSubj = length(results);

% run analysis over all subjects
subjIDs = nan(nSubj,1);
accuracies.lasso_min = nan(nSubj,1);
accuracies.ridge_min = nan(nSubj,1);
for i = 1 : nSubj
    subjIDs(i) = results{i}.subjID;
    % mean accuracy over cvb 
    accuracies.lasso_min(i) = mean(results{i}.lasso.accuracy.min);
    accuracies.ridge_min(i) = mean(results{i}.ridge.accuracy.min);
end
end

