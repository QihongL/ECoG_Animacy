function [featureSelected, numIterations, wtsVal] = ...
    iltpHelper_getWtsVals_t(results, nSubjs, curTimeIdx, numCVB, featureSelected)

numIterations = nan(nSubjs,1);
for s = 1 : nSubjs
    
    % preallocate number of iterations
    numIters = length(results{s}.lasso);
    numIterations(s) = numIters;
    
    wtsMagnitude_overIter = 0;
    wtsMagnitude_nz_overIter = 0;
    for i = 1 : numIters
        allcoef_acrossCV = [];
        for c = 1 : numCVB
            % record feature selected, collapse over CVB
            featureSelected{s}{curTimeIdx} = union(featureSelected{s}{curTimeIdx}, results{s}.lasso{i}.voxSel{c});
            
            % get weights over all cross validation
            wts = results{s}.lasso{i}.coef_min{c};
            wts(1) = [];    % remove the intercept
            % collect weight vector across CVB
            allcoef_acrossCV = vertcat(allcoef_acrossCV, results{s}.lasso{i}.coef_min{c});
            allcoef_acrossCV_nz = allcoef_acrossCV(allcoef_acrossCV ~= 0);
        end
        % accumulate average weight magnitude across iteration
        wtsMagnitude_overIter = wtsMagnitude_overIter + mean(abs(allcoef_acrossCV));
        wtsMagnitude_nz_overIter = wtsMagnitude_nz_overIter + mean(abs(allcoef_acrossCV_nz));
    end
    % put weight magnitude into a matrix
    wtsVal.all(s) = wtsMagnitude_overIter;
    wtsVal.nz(s) = wtsMagnitude_nz_overIter;
end

end