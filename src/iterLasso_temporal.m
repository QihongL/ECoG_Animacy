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
CHANCE = .5;

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
for t = 1 : 3 : numTimePts-1
    fprintf('T = %d: \n', t)
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
        fprintf('Subj%d: \n', s)
        % load the data
        load(strcat(DIR.DATA{t},filename.data{i}))
        y = metadata(s).targets(1).target;  % target 1 = label
        % read data parameters
        cvidx = metadata(s).cvind(:,CVCOL);
        
        %% fit iterative lasso 
        [results{i}] = runIterLasso(X, y, cvidx, options, CHANCE);
        results{i}.title = 'iterative_lasso_with_lambda_min';
        results{i}.subjID = s;
        results{i}.dataType = DATA_TYPE;
        results{i}.boxcar = BOXCAR;
        results{i}.windowStart = WIND_START(t).name;
        results{i}.windowSize = WIND_SIZE;
        
    end
    fprintf('\n');
    %% save the results
    % get final output directory
    finalOutDir = fullfile(DIR.OUT, 'BoxCar', BOXCAR, ...
        'WindowStart', WIND_START(t).name, 'WindowSize', WIND_SIZE, '/');
    if exist(finalOutDir,'dir') == 7
        %         warning('Final out dir exists!')
    else
        mkdir(finalOutDir)
    end
    % save the data
    saveFileName = sprintf(strcat('results_ilasso_', DATA_TYPE, '.mat'));
    save(strcat(finalOutDir,saveFileName), 'results')
    
end