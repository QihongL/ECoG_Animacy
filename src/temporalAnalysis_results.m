%% Temporal analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; close all;

% specify parameters
DATA_TYPE = 'ref'; % 'ref' OR 'raw'
CVCOL = 1;      % use the 1st column of cv idx for now
numCVB = 10;
options.nlambda = 100;

BOXCAR = '001';
WIND_SIZE = '0050';

% specify path information
DIR.PROJECT = '/Users/Qihong/Dropbox/github/ECOG_Manchester';
% point to the directory for the data
DIR.WIND_START = fullfile(DIR.PROJECT, 'results/temporal/BoxCar', BOXCAR, '/WindowStart/');
WIND_START = getAllDirNames(DIR.WIND_START);
DIR.DATA = getAllDataPath(DIR.WIND_START, WIND_START, WIND_SIZE);


nTimePts = length(DIR.DATA);
load(strcat(DIR.DATA{1}, 'results_', DATA_TYPE, '.mat'))
nSubjs = length(results); % assume all timePts have the same # of subj
% preallocate
accuracy.lasso.onese = zeros(nTimePts, nSubjs);
accuracy.lasso.min = zeros(nTimePts, nSubjs);


for t = 1 : nTimePts
    load(strcat(DIR.DATA{t}, 'results_', DATA_TYPE, '.mat'))
    for s = 1 : nSubjs
        accuracy.lasso.onese(t,s) = mean(results{s}.lasso.accuracy.onese);
        accuracy.lasso.min(t,s) = mean(results{s}.lasso.accuracy.min);
    end
end


%% plot the results 
FS = 14;
LW = 2; 

hold on 
plot(mean(accuracy.lasso.onese,2), 'linewidth', LW)
plot(mean(accuracy.lasso.min,2), 'linewidth', LW)
plot([1 nTimePts],[.5 .5], 'k--')
hold off

title_text = sprintf('Accuracy over time, data: %s',DATA_TYPE);
title(title_text, 'fontsize' , FS)
ylabel('Holdout accuracy', 'fontsize' , FS)
xlabel('Time (unit of 10ms)', 'fontsize' , FS)
leg = legend({'lasso(1se)','lasso(min)', 'chance level'}, 'location', 'southeast');
set(leg,'FontSize',FS);
