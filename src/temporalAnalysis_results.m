%% Temporal analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; close all;

% specify parameters
DATA_TYPES = {'raw','ref'};
DATA_TYPE = 'raw'; % 'ref' OR 'raw'
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

% get some common parameters
load(strcat(DIR.DATA{1}, 'results_', DATA_TYPE, '.mat'))

nTimePts = length(DIR.DATA);
nSubjs = length(results); % assume all timePts have the same # of subj




%% collecting accuracy scores over time 
% preallocate
accuracy.lasso.onese = nan(nTimePts, nSubjs);
accuracy.lasso.min = nan(nTimePts, nSubjs);

% loop over time 
for t = 1 : nTimePts-2
    load(strcat(DIR.DATA{t}, 'results_', DATA_TYPE, '.mat'))
    % loop over subjects
    for i = 1 : nSubjs
        accuracy.lasso.onese(t,i) = mean(results{i}.lasso.accuracy.onese);
        accuracy.lasso.min(t,i) = mean(results{i}.lasso.accuracy.min);
    end
end


%% plot the average accuracy
FS = 14;
LW = 2;
figure(1)
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

%% individual subject
figure(2)

for i = 1 : nSubjs
    subplot(2, nSubjs/2, i);
    hold on
    plot(accuracy.lasso.onese(:,i))
    plot(accuracy.lasso.min(:,i))
    plot([1 nTimePts],[.5 .5], 'k--')
    hold off
    title_text = sprintf('Accuracy, data: %s, subjID: %d', DATA_TYPE, results{i}.subjID);
    title(title_text, 'fontsize' , FS)
    ylabel('Holdout accuracy', 'fontsize' , FS)
    xlabel('Time (unit of 10ms)', 'fontsize' , FS)
    leg = legend({'lasso(1se)','lasso(min)', 'chance level'}, 'location', 'southeast');
    set(leg,'FontSize',FS);
    
end



