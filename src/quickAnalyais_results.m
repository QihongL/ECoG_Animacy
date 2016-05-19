%% Do some quick analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; 
% specify path information
% point to the output directory
DIR.OUT = '/Users/Qihong/Dropbox/github/ECOG_Manchester/results/allTimePts/';
% mvpa parameters
BOXCAR = '010';
WIND_START = '0000';
WIND_SIZE = '1000';

% read results, raw data
DATA_TYPE = 'raw';
resultName = strcat('results_', DATA_TYPE,'_bc',BOXCAR, ...
    '_wStart',WIND_START, 'wSize', WIND_SIZE, '.mat');
load(strcat(DIR.OUT, resultName))
[rawAcc, subjIDs.raw] = readResults_allTimePts(results);

% read results, raw data
DATA_TYPE = 'ref';
resultName = strcat('results_', DATA_TYPE,'_bc',BOXCAR, ...
    '_wStart',WIND_START, 'wSize', WIND_SIZE, '.mat');
load(strcat(DIR.OUT, resultName))
[refAcc, subjIDs.ref] = readResults_allTimePts(results);

subjIDs.common = intersect(subjIDs.raw,subjIDs.ref);
nSubjs.common = length(subjIDs.common);

% get index for the common subects (across dataType)
[~,subjIdx.raw]=ismember(subjIDs.common, subjIDs.raw);
[~,subjIdx.ref]=ismember(subjIDs.common, subjIDs.ref);


%% plot, all subjects 
figure(1);
subplot(1,2,1)
batMat = horzcat(rawAcc.lasso_min, rawAcc.ridge_min);

bar(batMat)
ylim([.5 1])
ax = gca;
ax.XTickLabel = subjIDs.raw;
FS = 14;

title_text = sprintf('Accuracy (chance = .5), data: %s \n [boxcar: %s, windowStart: %s, windowSize: %s]', ...
    'raw', BOXCAR, WIND_START, WIND_SIZE);
title(title_text, 'fontsize' , FS)
ylabel('Holdout classifcation accuracy', 'fontsize' , FS)
xlabel('Subject ID', 'fontsize' , FS)
leg = legend({'lasso(min)','ridge(min)'}, 'location', 'northwest');
set(leg,'FontSize',FS);

% 
subplot(1,2,2)
batMat = horzcat(refAcc.lasso_min, refAcc.ridge_min);
bar(batMat)
ylim([.5 1])
ax = gca;
ax.XTickLabel = subjIDs.ref;

FS = 14;
title_text = sprintf('Accuracy (chance = .5), data: %s \n [boxcar: %s, windowStart: %s, windowSize: %s]', ...
    'ref', BOXCAR, WIND_START, WIND_SIZE);
title(title_text, 'fontsize' , FS)
ylabel('Holdout classifcation accuracy', 'fontsize' , FS)
xlabel('Subject ID', 'fontsize' , FS)
leg = legend({'lasso(min)','ridge(min)'}, 'location', 'northwest');
set(leg,'FontSize',FS);

%% 
figure(2)
batMat = horzcat(rawAcc.lasso_min(subjIdx.raw), refAcc.lasso_min(subjIdx.ref),...
                rawAcc.ridge_min(subjIdx.raw), refAcc.ridge_min(subjIdx.ref));
bar(batMat)
ylim([.5 1])
ax = gca;
ax.XTickLabel = subjIDs.raw;
FS = 14;

 title_text = sprintf('Accuracy (chance = .5), \n [boxcar: %s, windowStart: %s, windowSize: %s]', ...
    BOXCAR, WIND_START, WIND_SIZE);
title(title_text, 'fontsize' , FS)
ylabel('Holdout classifcation accuracy', 'fontsize' , FS)
xlabel('Subject ID', 'fontsize' , FS)
leg = legend({'lasso(raw)','lasso(ref)', 'ridge(raw)','ridge(ref)'}, 'location', 'northwest');
set(leg,'FontSize',FS);


