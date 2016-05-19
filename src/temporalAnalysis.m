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


%% loop over all time points
numTimePts = length(WIND_START);
for t = 1 : numTimePts
    fprintf('T = %d:', t)
    % get filenames and IDs for all subjects
    [filename, subjIDs] = getFileNames(DIR.DATA{t}, DATA_TYPE);
    numSubjs = length(subjIDs); % number of subjects
    
    %% start the mvpa analysis
    results = cell(numSubjs,1);
    % loop over subjects
    for i = 1 : numSubjs
        s = subjIDs(i);
        fprintf(' %d', s)
        % preallocate: results{i} for the ith subject
        results{i}.subjID = s;
        results{i}.dataType = DATA_TYPE;
        results{i}.boxcar = BOXCAR;
        results{i}.windowStart = WIND_START(t).name;
        results{i}.windowSize = WIND_SIZE;
        
        % preallocate
        results{i}.lasso.accuracy.onese = nan(numCVB,1);
        results{i}.lasso.accuracy.min = nan(numCVB,1);
%         results{i}.ridge.accuracy.onese = nan(numCVB,1);
%         results{i}.ridge.accuracy.min = nan(numCVB,1);
        
        % load the data
        load(strcat(DIR.DATA{t},filename.metadata))
        load(strcat(DIR.DATA{t},filename.data{i}))
        y = metadata(s).targets(1).target;  % target 1 = label
        [M,N] = size(X);
        % read data parameters
        cvidx = metadata(s).cvind(:,CVCOL);
        
        % loop over CV blocks
        for c = 1: numCVB
            % choose a cv index
            testIdx = cvidx == c;
            
            % fit logistic models
            result = runRegularizedLogisticRegression(X, y, testIdx, options);
            % save performance
            results{i}.lasso.accuracy.onese(c) = result.lasso.accuracy.onese;
            results{i}.lasso.accuracy.min(c) = result.lasso.accuracy.min;
%             results{i}.ridge.accuracy.onese(c) = result.ridge.accuracy.onese;
%             results{i}.ridge.accuracy.min(c) = result.ridge.accuracy.min;
            
        end 
    end
    fprintf('\n');
    % get final output directory
    finalOutDir = fullfile(DIR.OUT, 'BoxCar', BOXCAR, ...
        'WindowStart', WIND_START(t).name, 'WindowSize', WIND_SIZE, '/');
    if exist(finalOutDir,'dir') == 7
%         warning('Final out dir exists!')
    else
        mkdir(finalOutDir)
    end
    % save the data
    saveFileName = sprintf( strcat('results_', DATA_TYPE, '.mat'));
    save(strcat(finalOutDir,saveFileName), 'results')
    
end