%% electode y dimention movement analysis
function [ weightedCoord ] = iltpHelper_getMeanY(metadata, allTps, curSubjID, runningSubjID, numCVB, yCoords, featureSelected)
    weightedCoord = cell(numCVB,1);
    
    for t = 1 : length(allTps)
        numSensors_bySubj = length(metadata(curSubjID).coords.labels);
        ycoords_centered = metadata(curSubjID).coords.xyz(:,2) - yCoords.median;
        % preallocate: number of times each electrode got selected
        electrodeWeights = zeros(numSensors_bySubj,1);
        % find out the start and end idx for each electrode window
        idx_t0 = metadata(curSubjID).WindowSizeInMilliseconds * (0 : numSensors_bySubj-1);
        idx_t50 = metadata(curSubjID).WindowSizeInMilliseconds * (1 : numSensors_bySubj);
        for e = 1 : numSensors_bySubj
            electrodeWeights(e) = sum(featureSelected{runningSubjID}{t} > idx_t0(e) & featureSelected{runningSubjID}{t}<= idx_t50(e));
        end
        % consider
        for i = 1 : numCVB
            thres_electrodeWeights = electrodeWeights >= i;
            weightedCoord{i}(t) = mean(ycoords_centered .* thres_electrodeWeights);
        end
        
    end


end

