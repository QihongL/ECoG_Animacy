function [ nSubjs, allSubjIDs, subjIdx ] = getSubjIdx_ref_raw( nSubjs, allSubjIDs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% common subjects
allSubjIDs.common = intersect(allSubjIDs.raw,allSubjIDs.ref);
nSubjs.common = length(allSubjIDs.common);
% get index for the common subects (across dataType)
[~,subjIdx.common.raw]=ismember(allSubjIDs.common, allSubjIDs.raw);
[~,subjIdx.common.ref]=ismember(allSubjIDs.common, allSubjIDs.ref);

% all subjects, collapse raw and ref, prefer ref
allSubjIDs.all = union(allSubjIDs.raw,allSubjIDs.ref);
nSubjs.all = length(allSubjIDs.all);
% get subject ids for the all subects
[~, subjIdx.all.ref] = ismember(allSubjIDs.ref,allSubjIDs.ref);
[~, subjIdx.all.raw] = ismember(setdiff(allSubjIDs.raw, allSubjIDs.ref), allSubjIDs.raw);

% subjects who have only ref or raw data 
[~, subjIdx.all.refOnly] = ismember(setdiff(allSubjIDs.ref, allSubjIDs.common), allSubjIDs.ref);
[~, subjIdx.all.rawOnly] = ismember(setdiff(allSubjIDs.raw, allSubjIDs.common), allSubjIDs.raw);

end

