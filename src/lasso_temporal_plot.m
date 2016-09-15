%% Temporal analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; clf;

% specify mvpa parameters
DATA_TYPES = {'raw','ref'};

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

% point to the directory for the metadata
DIR.META_WIND_START = strcat(DIR.PROJECT, '/data/ECoG/data/avg/BoxCar/', BOXCAR, '/WindowStart/');
DIR.METADATA = getAllDataPath(DIR.META_WIND_START, WIND_START, WIND_SIZE);


%% get some common parameters, for data processing purpose 
nTimePts = length(DIR.DATA);
% read subject info for the raw data
load(strcat(DIR.DATA{1}, 'results_basic_', DATA_TYPES{1}, '.mat'))
[nSubjs.raw, allSubjIDs.raw] = getSubjInfoFromResults(results);
% read subject info for the ref data
load(strcat(DIR.DATA{1}, 'results_basic_', DATA_TYPES{2}, '.mat'))
[nSubjs.ref, allSubjIDs.ref] = getSubjInfoFromResults(results);

% get subject idx, number with difference reference frame
[nSubjs, allSubjIDs, subjIdx] = getSubjIdx_ref_raw( nSubjs, allSubjIDs);


%% collecting accuracy scores over time
% preallocate
rawAccuracy.lasso_min = nan(nTimePts-TRIM, nSubjs.raw);
refAccuracy.lasso_min = nan(nTimePts-TRIM, nSubjs.ref);

% loop over time
for t = 1 : nTimePts-TRIM
    % process RAW data 
    
    DATA_TYPE = DATA_TYPES{1}; 
    % load metadata
    load(strcat(DIR.METADATA{t}, 'metadata_', DATA_TYPE, '.mat'))
    % loop over all subjects
    load(strcat(DIR.DATA{t}, 'results_basic_', DATA_TYPE, '.mat'))
    for s = 1 : nSubjs.raw
        rawAccuracy.lasso_min(t,s) = mean(results{s}.accuracy);
%         fprintf('RAW: beta dim: %d, num electrodes: %d \n ', ...
%             (round(length(results{s}.coef{1}) - 1)/ 100), ...
%             size(metadata(s).coords.xyz,1))
    end

    
    % process REF data 
    DATA_TYPE = DATA_TYPES{2};
    % load metadata
    load(strcat(DIR.METADATA{t}, 'metadata_', DATA_TYPE, '.mat'))
    % loop over subjects 
    load(strcat(DIR.DATA{t}, 'results_basic_', DATA_TYPE, '.mat'))
    for s = 1 : nSubjs.ref
        refAccuracy.lasso_min(t,s) = mean(results{s}.accuracy);
%         fprintf('REF: beta dim: %d, num electrodes: %d \n', ...
%             (length(results{s}.coef{1}) - 1)/ 100,size(metadata(s).coords.xyz,1))
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FS = 20;
LW = 2;
xrange = [0,170];

%% accuracy for individual subject
figure(1)
for i = 1 : nSubjs.common
    subplot(2, ceil(nSubjs.common/2), i);
    hold on
    plot(rawAccuracy.lasso_min(:,subjIdx.common.raw(i)), 'linewidth', LW)
    plot(refAccuracy.lasso_min(:,subjIdx.common.ref(i)), 'linewidth', LW)
    cor_raw_ref = corr(refAccuracy.lasso_min(:,subjIdx.common.ref(i)),rawAccuracy.lasso_min(:,subjIdx.common.raw(i)));
    plot([1 nTimePts],[.5 .5], 'k--')
    hold off
    ylim([0 1])
    xlim(xrange)
    title_text = sprintf('subjID: %d \n Correlation(raw,ref): %.2f', ...
        allSubjIDs.common(i), cor_raw_ref);
    title(title_text, 'fontsize' , FS)
    ylabel('Holdout accuracy', 'fontsize' , FS)
    xlabel('Time (unit of 10ms)', 'fontsize' , FS)
    leg = legend({'Raw','Ref', 'chance(.5)'}, 'location', 'southeast');
    set(leg,'FontSize',FS); 
end



%% plot the average accuracy
figure(3)
yrange = [0,1];

% compute the mean and the standard error 
alpha = .05; 
tval = tinv(1 - alpha/2, nSubjs.common-1);
average.raw = mean(rawAccuracy.lasso_min(:,subjIdx.common.raw),2);
average.ref = mean(refAccuracy.lasso_min(:,subjIdx.common.ref),2);
se.raw = std(rawAccuracy.lasso_min(:,subjIdx.common.raw)') * tval / sqrt(nSubjs.common);
se.ref = std(refAccuracy.lasso_min(:,subjIdx.common.ref)') * tval / sqrt(nSubjs.common);

% plot with error bars 
subplot(2,1,1)
hold on 
errorbar(1:nTimePts-TRIM,average.raw,se.raw)
errorbar(1:nTimePts-TRIM,average.ref,se.ref)
plot([1 nTimePts],[.5 .5], 'k--')
hold off 
xlim(xrange)
ylim(yrange)

% text 
title_text = sprintf('Lasso accuracy over time, averaged across %d common subjects', nSubjs.common);
title(title_text, 'fontsize' , FS)
ylabel('Holdout accuracy', 'fontsize' , FS)
xlabel('Time (unit of 10ms)', 'fontsize' , FS)
leg = legend({'Raw','Ref', 'chance(.5)'}, 'location', 'southeast');
set(leg,'FontSize',FS);

% 
% % plot withour error bars 
% subplot(3,1,2)
% hold on
% plot(mean(rawAccuracy.lasso_min(:,subjIdx.common.raw),2), 'linewidth', LW)
% plot(mean(refAccuracy.lasso_min(:,subjIdx.common.ref),2), 'linewidth', LW)
% plot([1 nTimePts],[.5 .5], 'k--')
% hold off
% xlim(xrange)
% ylim(yrange)
% title_text = sprintf('Logistic Lasso accuracy, averaged across %d common subjects', nSubjs.common);
% % title(title_text, 'fontsize' , FS)
% ylabel('Holdout accuracy', 'fontsize' , FS)
% xlabel('Time (unit of 10ms)', 'fontsize' , FS)
% leg = legend({'Raw','Ref', 'chance(.5)'}, 'location', 'southeast');
% set(leg,'FontSize',FS);

%% mean accuracy over time, collapse raw and ref, prefer ref
% figure(4)
% collapse raw and ref 
collapAccData = horzcat(refAccuracy.lasso_min(:,subjIdx.all.ref),...
    rawAccuracy.lasso_min(:,subjIdx.all.raw));
% compute mean and std 
average.collapAcc = mean(collapAccData,2);
se.collapAcc = std(collapAccData,0,2) * tval / sqrt(nSubjs.all);

tval = tinv(1 - alpha/2, nSubjs.all-1);

% plot with error bar 
subplot(2,1,2)
hold on 
errorbar(1:nTimePts-TRIM, average.collapAcc,se.collapAcc,'LineWidth',1.5)
plot([1 nTimePts],[.5 .5], 'k--')
hold off 
xlim(xrange)
ylim(yrange)

title_text = sprintf('Logistic Lasso accuracy, averaged across all %d subjects \n collapsed across ref & raw (prefer ref)', nSubjs.all);
title(title_text, 'fontsize' , FS)
ylabel('Holdout accuracy', 'fontsize' , FS)
xlabel('Time (unit of 10ms)', 'fontsize' , FS)
leg = legend({'collapsed accuracy ', 'chance(.5)'}, 'location', 'northeast');
set(leg,'FontSize',FS);


% %% plot on the poster (big font, bold, etc...)
% figure(5)
% hold on 
% errorbar(1:nTimePts-TRIM,average.collapAcc,se.collapAcc,'LineWidth',1.5, 'color', 'b')
% plot([1 nTimePts],[.5 .5], 'k--','LineWidth',3)
% plot(1:nTimePts-TRIM,average.collapAcc,'LineWidth',3, 'color', 'b')
% 
% hold off 
% xlim([0, 163])
% ylim(yrange)
% 
% title_text = sprintf('Logistic Lasso accuracy, averaged across all %d subjects \n collapsed across ref & raw (prefer ref)', nSubjs.all);
% % title(title_text, 'fontsize' , FS)
% ylabel('Holdout accuracy', 'fontsize' , FS, 'FontWeight','bold')
% xlabel('Time (unit: 10ms)', 'fontsize' , FS, 'FontWeight','bold')
% leg = legend({'Accuracy ', 'Chance(.5)'}, 'location', 'northeast');
% set(leg,'FontSize',FS);
% set(gca,'FontSize',FS)
% 
% % % plot withour error bar 
% % subplot(2,1,2)
% % hold on 
% % plot(average.collapAcc, 'linewidth', LW)
% % plot([1 nTimePts],[.5 .5], 'k--')
% % hold off 
% % 
% % xlim(xrange)
% % ylim(yrange)



