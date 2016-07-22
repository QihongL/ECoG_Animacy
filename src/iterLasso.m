%% Do some quick analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; close all;

% path parameters
BOXCAR = '010';
WIND_START = '0200';
WIND_SIZE = '1000';

% specify path information
% point to the project dir
% DIR.PROJECT = '/Users/Qihong/Dropbox/github/ECOG_Manchester';
% point to the directory for the data
DIR.DATA = strcat('/Users/Qihong/Dropbox/github/ECOG_Manchester/data/ECoG/data/avg/BoxCar/', ...
    BOXCAR, '/WindowStart/', WIND_START, '/WindowSize/', WIND_SIZE, '/');
% point to the output directory
DIR.OUT = '/Users/Qihong/Dropbox/github/ECOG_Manchester/results/allTimePts/';

% specify parameters
DATA_TYPE = 'ref'; % 'ref' OR 'raw'
CVCOL = 1;      % use the 1st column of cv idx for now
numCVB = 10;
options.nlambda = 100;

% constant
chance = .5;
alpha = .05;
% maxIter = 10;

% get filenames and IDs for all subjects
[filename, subjIDs] = getFileNames(DIR.DATA, DATA_TYPE);
numSubjs = length(subjIDs); % number of subjects

for i = 1 : length(filename.data)
    fprintf(filename.data{i})
end

%% start the mvpa analysis
results = cell(numSubjs,1);
% loop over subjects
for i = 1 : numSubjs
% for i = 1
    s = subjIDs(i);
    fprintf('Subj%.2d: \n', s);
    % load the data
    load(strcat(DIR.DATA,filename.metadata))
    load(strcat(DIR.DATA,filename.data{i}))
    y = metadata(s).targets(1).target;  % target 1 = label
    [M,N] = size(X);
    % read data parameters
    cvidx = metadata(s).cvind(:,CVCOL);
    
    
    % preallocate: results{i} for the ith subject
    results{i}.subjID = s;
    results{i}.dataType = DATA_TYPE;
    
    
    %% iterative lasso
    allUsedVox = cell(numCVB,1);
    unusedVoxIdx = repmat({(1:N)'},[numCVB,1]);
    iter = 1;
    while true
        fprintf('Iter%d: ', iter);

        % loop over CV blocks
        for c = 1: numCVB
            fprintf('%d ', c);
            
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
                unusedVoxIdx{c} = setdiff(unusedVoxIdx{c}, results{i}.lasso{iter-1}.voxSel{c});
            end
            
            
            % fit lasso
            options.alpha = 1; % 1 == lasso, 0 == ridge
            cvfit = cvglmnet(X_train, y_train, 'binomial', options);
            
            % compute the test set accuracy
            y_hat = myStepFunction(cvglmnetPredict(cvfit, X_test,cvfit.lambda_min));
            accuracy = sum(y_hat == y_test) / length(y_test);
            
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
                error('?')
            end
            
        end
        
        
        %% record the measures
        results{i}.lasso{iter}.lambda_min = min_lambdas;
        results{i}.lasso{iter}.coef_min = min_lambda_coefs;
        results{i}.lasso{iter}.accuracy = accuracies;
        results{i}.lasso{iter}.voxSel = correctedVoxSel;
        
        %% stopping criterion
        significant = ttest(results{i}.lasso{iter}.accuracy , chance, 'Alpha', alpha);
        if ~significant
            fprintf('|  Accuracy = %.4f, with , no better than chance, terminate.\n',...
                mean(results{i}.lasso{iter}.accuracy));
            results{i}.iterNum = iter;
            break;
        else
            fprintf('|  Accuracy = %.4f*, continue...\n',...
                mean(results{i}.lasso{iter}.accuracy));
            iter = iter + 1;
        end
    end % end of iterative lasso
    
    
    
end
% save the data
saveFileName = sprintf(strcat('results_ilasso_', DATA_TYPE,'_bc',BOXCAR, ...
    '_wStart',WIND_START, 'wSize', WIND_SIZE, '.mat'));

save(strcat(DIR.OUT,saveFileName), 'results')
