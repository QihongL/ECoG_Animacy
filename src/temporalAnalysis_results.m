%% Temporal analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; close all;

% specify mvpa parameters
DATA_TYPES = {'raw','ref'};
DATA_TYPE = 'raw'; % 'ref' OR 'raw'
CVCOL = 1;      % use the 1st column of cv idx for now
numCVB = 10;
options.nlambda = 100;
BOXCAR = '001';
WIND_SIZE = '0050';
TRIM = 2; % TODO fix this! 

% specify path information
DIR.PROJECT = '/Users/Qihong/Dropbox/github/ECOG_Manchester';
% point to the directory for the data
DIR.WIND_START = fullfile(DIR.PROJECT, 'results/temporal/BoxCar', BOXCAR, '/WindowStart/');
WIND_START = getAllDirNames(DIR.WIND_START);
DIR.DATA = getAllDataPath(DIR.WIND_START, WIND_START, WIND_SIZE);

% get some common parameters, for data processing purpose 
nTimePts = length(DIR.DATA);
% read subject info for the raw data
load(strcat(DIR.DATA{1}, 'results_', DATA_TYPES{1}, '.mat'))
[nSubjs.raw, allSubjIDs.raw] = getSubjInfoFromResults(results);
% read subject info for the ref data
load(strcat(DIR.DATA{1}, 'results_', DATA_TYPES{2}, '.mat'))
[nSubjs.ref, allSubjIDs.ref] = getSubjInfoFromResults(results);
allSubjIDs.common = intersect(allSubjIDs.raw,allSubjIDs.ref);
nSubjs.common = length(allSubjIDs.common);

% get index for the common subects (across dataType)
[~,subjIdx.raw]=ismember(allSubjIDs.common, allSubjIDs.raw);
[~,subjIdx.ref]=ismember(allSubjIDs.common, allSubjIDs.ref);


%% collecting accuracy scores over time
% preallocate
rawAccuracy.lasso_min = nan(nTimePts-TRIM, nSubjs.raw);
refAccuracy.lasso_min = nan(nTimePts-TRIM, nSubjs.ref);

% loop over time
for t = 1 : nTimePts-TRIM
    % loop over subjects
    load(strcat(DIR.DATA{t}, 'results_', 'raw', '.mat'))
    for s = 1 : nSubjs.raw
        rawAccuracy.lasso_min(t,s) = mean(results{s}.lasso.accuracy.min);
    end
    load(strcat(DIR.DATA{t}, 'results_', 'ref', '.mat'))
    for s = 1 : nSubjs.ref
        refAccuracy.lasso_min(t,s) = mean(results{s}.lasso.accuracy.min);
    end
end



%% plot 
FS = 14;
LW = 2;

%% accuracy for individual subject
figure(1)
for i = 1 : nSubjs.common
    subplot(2, ceil(nSubjs.common/2), i);
    hold on
    plot(rawAccuracy.lasso_min(:,subjIdx.raw(i)))
    plot(refAccuracy.lasso_min(:,subjIdx.ref(i)))
    plot([1 nTimePts],[.5 .5], 'k--')
    hold off
    
    title_text = sprintf('Accuracy, data: %s, subjID: %d', DATA_TYPE, allSubjIDs.common(i));
    title(title_text, 'fontsize' , FS)
    ylabel('Holdout accuracy', 'fontsize' , FS)
    xlabel('Time (unit of 10ms)', 'fontsize' , FS)
    leg = legend({'Raw','Ref', 'chance(.5)'}, 'location', 'southeast');
    set(leg,'FontSize',FS); 
end


%% plot the average accuracy
figure(2)
xrange = [0,170];
yrange = [.4,.9];

% plot with error bars 
subplot(2,1,1)
alpha = .975; 
tval = tinv(.975,nSubjs.common-1);
average.raw = mean(rawAccuracy.lasso_min(:,subjIdx.raw),2);
average.ref = mean(refAccuracy.lasso_min(:,subjIdx.ref),2);
stdev.raw = std(rawAccuracy.lasso_min(:,subjIdx.raw)');
stdev.ref = std(refAccuracy.lasso_min(:,subjIdx.ref)');
xlim(xrange)
ylim(yrange)
hold on 
errorbar(1:nTimePts-TRIM,average.raw,stdev.raw)
errorbar(1:nTimePts-TRIM,average.ref,stdev.ref)
plot([1 nTimePts],[.5 .5], 'k--')
hold off 

% text 
title_text = sprintf('Lasso accuracy over time, averaged across %d common subjects', nSubjs.common);
title(title_text, 'fontsize' , FS)
ylabel('Holdout accuracy', 'fontsize' , FS)
xlabel('Time (unit of 10ms)', 'fontsize' , FS)
leg = legend({'Raw','Ref', 'chance(.5)'}, 'location', 'southeast');
set(leg,'FontSize',FS);


% plot withour error bars 
subplot(2,1,2)
hold on
plot(mean(rawAccuracy.lasso_min(:,subjIdx.raw),2), 'linewidth', LW)
plot(mean(refAccuracy.lasso_min(:,subjIdx.ref),2), 'linewidth', LW)
plot([1 nTimePts],[.5 .5], 'k--')
hold off
xlim(xrange)
ylim(yrange)
title_text = sprintf('Lasso accuracy over time, averaged across %d common subjects', nSubjs.common);
% title(title_text, 'fontsize' , FS)
ylabel('Holdout accuracy', 'fontsize' , FS)
xlabel('Time (unit of 10ms)', 'fontsize' , FS)
leg = legend({'Raw','Ref', 'chance(.5)'}, 'location', 'southeast');
set(leg,'FontSize',FS);




