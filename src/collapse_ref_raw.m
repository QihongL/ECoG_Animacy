function [ collapsedData ] = collapse_ref_raw(data, allSubjIDs, subjIdx, missingSubjIdx )


% assume only one subject is missing, which is true for our ECOG data set
idx = subjIdx.all.ref;
cols = allSubjIDs.ref(idx) - (allSubjIDs.ref(idx) > missingSubjIdx);
collapsedData(:,cols) = data.ref(:,idx);

idx = subjIdx.all.raw;
cols = allSubjIDs.raw(idx) - (allSubjIDs.raw(idx) > missingSubjIdx);
collapsedData(:,cols) = data.raw(:,idx);


end

