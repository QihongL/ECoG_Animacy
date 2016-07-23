%% Iterative lasso for the Manchster ECoG data
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
options.alpha = 1; % 1 == lasso, 0 == ridge

% constant
CHANCE = .5;

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
    
    % read data parameters
    cvidx = metadata(s).cvind(:,CVCOL);
    
    %% run iterative lasso 
    [results{i}] = runIterLasso(X, y, cvidx, options, CHANCE);
    results{i}.subjID = s;
    results{i}.dataType = DATA_TYPE;
    
end
% save the data
saveFileName = sprintf(strcat('results_ilasso_', DATA_TYPE,'_bc',BOXCAR, ...
    '_wStart',WIND_START, 'wSize', WIND_SIZE, '.mat'));

% save(strcat(DIR.OUT,saveFileName), 'results')
