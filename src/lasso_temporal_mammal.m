%% Temporal analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; close all;

% specify parameters
DATA_TYPE = 'raw'; % 'ref' OR 'raw'
CVCOL = 1;      % use the 1st column of cv idx for now
numCVB = 10;
options.nlambda = 100;

BOXCAR = '001';
WIND_SIZE = '0050';

% specify path information
DIR.PROJECT = '/Users/Qihong/Dropbox/github/ECOG_Manchester/';
% point to the directory for the data
DIR.WIND_START = strcat(DIR.PROJECT, '/data/ECoG/data/avg/BoxCar/', BOXCAR, '/WindowStart/');
WIND_START = getAllDirNames(DIR.WIND_START);
DIR.DATA = getAllDataPath(DIR.WIND_START, WIND_START, WIND_SIZE);

% point to the output directory
DIR.OUT = strcat(DIR.PROJECT, 'results/temporal/');

% 
TEMPPATH.basicLabels = fullfile(DIR.PROJECT, 'data/labels', 'labels_allObjs.mat');
load(TEMPPATH.basicLabels)


%% loop over all time points
numTimePts = length(WIND_START);
for t = 1 : numTimePts
    fprintf('T = %d:', t)
    % get filenames and IDs for all subjects
    [filename, subjIDs] = getFileNames(DIR.DATA{t}, DATA_TYPE);
    numSubjs = length(subjIDs); % number of subjects
    % load metadata
    load(strcat(DIR.DATA{t},filename.metadata))
    
    %% start the mvpa analysis
    results = cell(numSubjs,1);
    % loop over subjects
    for i = 1 : numSubjs
        s = subjIDs(i);
        fprintf(' %d', s)
        % preallocate: results{i} for the ith subject
        results{i}.title = 'lasso_lambda_min';
        results{i}.subjID = s;
        results{i}.dataType = DATA_TYPE;
        results{i}.boxcar = BOXCAR;
        results{i}.windowStart = WIND_START(t).name;
        results{i}.windowSize = WIND_SIZE;
        
        % preallocate
        results{i}.accuracy = nan(numCVB,1);
        results{i}.lambda_min = nan(numCVB,1);
        results{i}.coef = cell(numCVB,1);
        
        % load the data
        load(strcat(DIR.DATA{t},filename.data{i}))
%         y = metadata(s).targets(1).target;  % target 1 = label
        y = labels_allObjs.isMammal(labels_allObjs.isAnimal);
        if sum(labels_allObjs.isMammal(labels_allObjs.isAnimal))~= 25 
            error('???')
        end
        X = X(labels_allObjs.isAnimal, :);
        [M,N] = size(X);
        
        % read data parameters
        cvidx = metadata(s).cvind(:,CVCOL);
        
        % loop over CV blocks
        for c = 1: numCVB
            % choose a cv index
            testIdx = cvidx == c;
            testIdx = testIdx(labels_allObjs.isAnimal);
            % fit logistic models
            result = runLasso(X, y, testIdx, options);
            % save performance
            results{i}.accuracy(c) = result.lasso_accuracy_lambda_min;
            results{i}.lambda_min(c) = result.lasso_lambda_min;
            results{i}.coef{c} = result.lasso_coef_lambda_min;

        end 
    end
    fprintf('\n');
    %% save the results 
%     get final output directory
    finalOutDir = fullfile(DIR.OUT, 'BoxCar', BOXCAR, ...
        'WindowStart', WIND_START(t).name, 'WindowSize', WIND_SIZE, '/');
%     if exist(finalOutDir,'dir') == 7
% %         warning('Final out dir exists!')
%     else
%         mkdir(finalOutDir)
%     end
    % save the data
    saveFileName = sprintf( strcat('results_basic_', DATA_TYPE, '.mat'));
    save(strcat(finalOutDir,saveFileName), 'results')
    
end