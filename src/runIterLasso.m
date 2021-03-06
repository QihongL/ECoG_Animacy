function [result] = runIterLasso(X, y, cvidx, options, chance)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
printInfoToConsole = true;
[~,N] = size(X);
alpha = .05;
numCVB = max(cvidx); % BAD WAY OF DECIDING numCVB...


%% iterative lasso
allUsedVox = cell(numCVB,1);
unusedVoxIdx = repmat({(1:N)'},[numCVB,1]);
iter = 1;
while true
    if printInfoToConsole
        fprintf('Iter%d: ', iter);
    end
    
    % loop over CV blocks
    for c = 1: numCVB
        
        %% feature processing
        % choose a cv index
        testIdx = cvidx == c;
        % hold out the test set
        X_train = X(~testIdx,:);
        X_test = X(testIdx,:);
        y_train = y(~testIdx);
        y_test = y(testIdx);
        
        % mask out selected voxel!
        if iter > 1
            X_train(:,allUsedVox{c}) = [];
            X_test(:,allUsedVox{c}) = [];
            % mask out the voxel selected last time;
            unusedVoxIdx{c} = setdiff(unusedVoxIdx{c}, result.lasso{iter-1}.voxSel{c});
        end
        
        %% fit logistic lasso
        cvfit = cvglmnet(X_train, y_train, 'binomial', options);
        
        % compute the test set accuracy
        y_hat = myStepFunction(cvglmnetPredict(cvfit, X_test,cvfit.lambda_min));
        accuracy = sum(y_hat == y_test) / length(y_test);
        
        %% recording
        % record model parameters, performance etc.
        coef = cvglmnetCoef(cvfit, 'lambda_min');
        
        % record lambda, beta, cv accuracy
        min_lambdas(c) = cvfit.lambda_min;
        min_lambda_coefs{c} = coef;
        accuracies(c) = accuracy;
        voxSel{c} = find(coef(2:end) ~= 0);
        
        % record voxels ever selected
        if iter == 1
            allUsedVox{c} = unusedVoxIdx{c}(voxSel{c});
            correctedVoxSel{c} = voxSel{c};
        elseif iter > 1
            if ~isempty(intersect(allUsedVox{c}, unusedVoxIdx{c}(voxSel{c})))
                error('ERROR: The same voxel selected twice!');
            end
            correctedVoxSel{c} = unusedVoxIdx{c}(voxSel{c});
            allUsedVox{c} = union(allUsedVox{c}, correctedVoxSel{c});
        else
            error('Iteration number must >= 1')
        end
        
    end
    
    %% record the measures
    result.lasso{iter}.lambda_min = min_lambdas;
    result.lasso{iter}.coef_min = min_lambda_coefs;
    result.lasso{iter}.accuracy = accuracies;
    result.lasso{iter}.voxSel = correctedVoxSel;
    
    for c = 1 : numCVB
        numVox(c) = length(result.lasso{iter}.voxSel{c});
    end
    
    %% stopping criterion
    significant = ttest(result.lasso{iter}.accuracy, chance, 'Tail', 'right', 'Alpha', alpha);
    if mean(result.lasso{iter}.accuracy) > chance * 1.1 && significant
        if printInfoToConsole
            fprintf('Accuracy = %.3f*, featureNum: %d~%d (mean = %.2f), continue...\n',...
                mean(result.lasso{iter}.accuracy), min(numVox),max(numVox), mean(min(numVox)));
        end
        iter = iter + 1;
    else
        if printInfoToConsole
            fprintf('Accuracy = %.3f, featureNum: %d~%d (mean = %.2f), no better than chance, terminate.\n',...
                mean(result.lasso{iter}.accuracy), min(numVox),max(numVox), mean(min(numVox)));
        end
        result.iterNum = iter;
        break;
    end
    
    %     if ~significant
    %         if printInfoToConsole
    %             fprintf('Accuracy = %.3f, featureNum: %d ~ %d (mean = %d), no better than chance, terminate.\n',...
    %                 mean(result.lasso{iter}.accuracy), min(numVox),max(numVox), mean(min(numVox)));
    %         end
    %         result.iterNum = iter;
    %         break;
    %     else
    %         if printInfoToConsole
    %             fprintf('Accuracy = %.3f*, featureNum: %d ~ %d (mean = %.3f), continue...\n',...
    %                 mean(result.lasso{iter}.accuracy), min(numVox),max(numVox), mean(min(numVox)));
    %         end
    %         iter = iter + 1;
    %     end
end % end of iterative lasso

end

