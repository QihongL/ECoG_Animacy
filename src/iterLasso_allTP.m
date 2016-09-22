%% Iterative lasso for the Manchster ECoG data
clear variables; clc; close all;
%% specify parameters
DATA_TYPE = 'ref'; % 'ref' OR 'raw'
CVCOL = 1;      % use the 1st column of cv idx for now
numCVB = 10;
options.nlambda = 100;
options.alpha = 1; % 1 == lasso, 0 == ridge
saveFileName_prefix = 'results_mam_obj_';
saveResultsFile = 1; 

% constant
CHANCE = .5;

%% specify path information

% path parameters
BOXCAR = '010';
WIND_START = '0200';
WIND_SIZE = '1000';

% point to the project dir
% DIR.PROJECT = '/Users/Qihong/Dropbox/github/ECOG_Manchester';
% point to the directory for the data
DIR.DATA = strcat('/Users/Qihong/Dropbox/github/ECOG_Manchester/data/ECoG/data/avg/BoxCar/', ...
    BOXCAR, '/WindowStart/', WIND_START, '/WindowSize/', WIND_SIZE, '/');
% point to the output directory
DIR.OUT = '/Users/Qihong/Dropbox/github/ECOG_Manchester/results/allTimePts/';


% basic level classification labels 
TEMPPATH.basicLabels = fullfile('/Users/Qihong/Dropbox/github/ECOG_Manchester/data/labels/labels_allObjs.mat');
load(TEMPPATH.basicLabels)
labels_allObjs.num.isMammal = find(labels_allObjs.isMammal==1);
labels_allObjs.num.isAnimal = find(labels_allObjs.isAnimal==1);
labels_mam_obj = vertcat(labels_allObjs.num.isMammal, [51:75]');

%% get filenames and IDs for all subjects
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
    
%     % select basic level chuck 
%     X = X(~metadata(s).targets(1).target,:);
%     y = labels_allObjs.isMammal(~metadata(s).targets(1).target);
%     cvidx = cvidx(~metadata(s).targets(1).target);
    
    % compare mammal versus objects 
    y = vertcat(ones(25,1), zeros(25,1));
    X = X(labels_mam_obj, :);
    cvidx = cvidx(labels_mam_obj);
    
    
    %% run iterative lasso 
    [results{i}] = runIterLasso(X, y, cvidx, options, CHANCE);
    results{i}.subjID = s;
    results{i}.dataType = DATA_TYPE;
    
end
% save the data
saveFileName = sprintf(strcat(saveFileName_prefix, DATA_TYPE,'_bc',BOXCAR, ...
    '_wStart',WIND_START, 'wSize', WIND_SIZE, '.mat'));

finalOutDir = '/Users/Qihong/Dropbox/github/ECOG_Manchester/results/allTimePts';
if saveResultsFile 
    save(fullfile(finalOutDir,saveFileName), 'results')
end
