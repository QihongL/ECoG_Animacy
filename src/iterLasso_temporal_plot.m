%% Temporal analysis for the Manchster ECoG data using iterative lasso
% plot the data from iterative lasso
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
load(strcat(DIR.DATA{1}, 'results_', DATA_TYPES{1}, '.mat'))
[nSubjs.raw, allSubjIDs.raw] = getSubjInfoFromResults(results);
% read subject info for the ref data
load(strcat(DIR.DATA{1}, 'results_', DATA_TYPES{2}, '.mat'))
[nSubjs.ref, allSubjIDs.ref] = getSubjInfoFromResults(results);

% get subject idx, number with difference reference frame
[nSubjs, allSubjIDs, subjIdx] = getSubjIdx_ref_raw( nSubjs, allSubjIDs);


%% specifiy time points with data
disp('All time points computed:')
allTps = 1 : 3:  163

%% resource preallocation
rawAccuracy = nan(nTimePts-TRIM, nSubjs.raw);
refAccuracy = nan(nTimePts-TRIM, nSubjs.ref);

% preallocate: number of lasso iterations (time X subject)
numIterations.raw = nan(length(allTps), nSubjs.raw);
numIterations.ref = nan(length(allTps), nSubjs.ref);
% preallocate: feature selected (time X 1)
featureSelected.raw = cell(nSubjs.raw,1);
featureSelected.ref = cell(nSubjs.ref,1);
for s = 1 : nSubjs.ref
    featureSelected.ref{s} = cell(length(allTps),1);
end
for s = 1 : nSubjs.raw
    featureSelected.raw{s} = cell(length(allTps),1);
end
% preallocate: weight magnitude (time X subject)
wtsVal.ref.all = nan(length(allTps), nSubjs.ref);
wtsVal.ref.nz = nan(length(allTps), nSubjs.ref);
wtsVal.raw.all = nan(length(allTps), nSubjs.raw);
wtsVal.raw.nz = nan(length(allTps), nSubjs.raw);

% 
nFeatureSelected.ref = nan(length(allTps), nSubjs.ref);
propFeatureSelected.ref = nan(length(allTps), nSubjs.ref);

%% gather data over time
for t = allTps
    curTimeIdx = find(allTps==t);
    
    %
    %     %% process RAW data
    %     DATA_TYPE = DATA_TYPES{1};
    %     % load metadata
    %     load(strcat(DIR.METADATA{t}, 'metadata_', DATA_TYPE, '.mat'))
    %     % loop over subjects
    %     load(strcat(DIR.DATA{t}, 'results_ilasso_', DATA_TYPE, '.mat'))
    
    
    
    
    %% process REF data
    DATA_TYPE = DATA_TYPES{2};
    % load metadata
    load(strcat(DIR.METADATA{t}, 'metadata_', DATA_TYPE, '.mat'))
    % loop over subjects
    load(strcat(DIR.DATA{t}, 'results_ilasso_', DATA_TYPE, '.mat'))
    
    
    %% gather weights values and number of iterations
    [featureSelected.ref, numIterations_t.ref, wtsVal_t.ref] = ...
        iltpHelper_getWtsVals_t(results, nSubjs.ref, curTimeIdx, numCVB, featureSelected.ref);
    % unpack data
    numIterations.ref(curTimeIdx,:) = numIterations_t.ref;
    wtsVal.ref.all(curTimeIdx,:) = wtsVal_t.ref.all;
    wtsVal.ref.nz(curTimeIdx,:) = wtsVal_t.ref.nz;
    
    
    %% get number/proportion of features & selected features
    for s = 1 : nSubjs.ref
        % total number of columns (features)
        numFeatures.ref(s) = metadata(allSubjIDs.ref(s)).ncol;
        % feature selection measures
        nFeatureSelected.ref(curTimeIdx,s) = length(featureSelected.ref{s}{curTimeIdx});
        propFeatureSelected.ref(curTimeIdx,s) = nFeatureSelected.ref(curTimeIdx,s)/numFeatures.ref(s);
    end
end


%% compute median y coordinate over all subjects (antierior-postierior)
yCoords.all = [];
for s = 1 : nSubjs.ref
    yCoords.all = vertcat(yCoords.all, metadata(allSubjIDs.ref(s)).coords.xyz(:,2));
end
yCoords.median = median(yCoords.all);

% % show the distribution of y coordinate
% hist(yCoords.all)
% xlabel('Postierior <---> Antierior', 'fontsize', 14)
% ylabel('Number of electrodes', 'fontsize', 14)
% title('Distribution of electrodes across subjects', 'fontsize', 14)

%% gather data for plotting
legend_subjIDs = cell(nSubjs.ref,1);

weightedCoord = cell(numCVB,1);
% weightedCoord = nan(length(allTps),nSubjs.ref);
for s = 1 : nSubjs.ref
    curSubjID = allSubjIDs.ref(s);
    % get the legend for subject id
    legend_subjIDs{s} = num2str(curSubjID);
    
    %% gather the Y coordinates for the selected features, thresholded by nCVB
    weightedCoord_singleSub = iltpHelper_getMeanY(metadata, allTps, curSubjID,...
        s, numCVB, yCoords, featureSelected.ref);
    % unpack the data
    for c = 1 : numCVB
        weightedCoord{c}(:,s) = weightedCoord_singleSub{c};
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.FS = 14;
p.LW = 2;


%% proportion of features selected over time
figure(1)

subplot(3,1,1)
plot(allTps, propFeatureSelected.ref, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Proportion features selected over time, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('% features selected', 'fontsize', p.FS)


subplot(3,1,2)
plot(allTps, movingmean(propFeatureSelected.ref, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('% features selected', 'fontsize', p.FS)


subplot(3,1,3)
se = tinv(.975, nSubjs.ref-1) * std(propFeatureSelected.ref, 0, 2) / sqrt(nSubjs.ref);
hold on
errorbar(allTps, mean(propFeatureSelected.ref,2), se, 'linewidth', p.LW)
plot([0 allTps(end)],[0 0], '--k')
hold off
xlim([0 allTps(end)])
title('Averaged across subjects', 'fontsize', p.FS);
ylabel('% features selected', 'fontsize', p.FS)
xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)


%% number of lasso iterations over time (indicate multi-collinearity)
figure(2)


subplot(3,1,1)
plot(allTps, numIterations.ref, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Feature multicollinearity over time , dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Number of lasso iterations', 'fontsize', p.FS)


subplot(3,1,2)
plot(allTps, movingmean(numIterations.ref, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Number of lasso iterations', 'fontsize', p.FS)

subplot(3,1,3)
se = tinv(.975, nSubjs.ref-1) * std(numIterations.ref, 0, 2) / sqrt(nSubjs.ref);
hold on
errorbar(allTps, mean(numIterations.ref,2), se, 'linewidth', p.LW)
plot([0 allTps(end)],[0 0], '--k')
hold off
xlim([0 allTps(end)])
title('Averaged across subjects', 'fontsize', p.FS);
ylabel('Number of lasso iterations', 'fontsize', p.FS)
xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)



%% average weight magnitude
figure(3)

subplot(3,2,1)
plot(allTps, wtsVal.ref.all, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Weights magnitude over time, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of all features', 'fontsize', p.FS)


subplot(3,2,3)
plot(allTps, movingmean(wtsVal.ref.all, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of all features', 'fontsize', p.FS)

subplot(3,2,5)
se = tinv(.975, nSubjs.ref-1) * std(wtsVal.ref.all, 0, 2) / sqrt(nSubjs.ref);
hold on
errorbar(allTps, mean(wtsVal.ref.all,2), se, 'linewidth', p.LW)
plot([0 allTps(end)],[0 0], '--k')
hold off
xlim([0 allTps(end)])
title('Averaged across subjects', 'fontsize', p.FS);
ylabel('Average weights of all features', 'fontsize', p.FS)
xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)

subplot(3,2,2)
plot(allTps, wtsVal.ref.nz, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Weights magnitude over time, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of non-zero features', 'fontsize', p.FS)


subplot(3,2,4)
plot(allTps, movingmean(wtsVal.ref.nz, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of non-zero features', 'fontsize', p.FS)

subplot(3,2,6)
se = tinv(.975, nSubjs.ref-1) * std(wtsVal.ref.nz, 0, 2) / sqrt(nSubjs.ref);
hold on
errorbar(allTps, mean(wtsVal.ref.nz,2), se, 'linewidth', p.LW)
plot([0 allTps(end)],[0 0], '--k')
hold off
xlim([0 allTps(end)])
title('Averaged across subjects', 'fontsize', p.FS);
ylabel('Average weights of non-zero features', 'fontsize', p.FS)
xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)

%% senor distribution along the y axis
tempSubPlotIdx = 1;
figure(4)
increment = 3;
for i = 2 : 2: numCVB
    subplot(5,increment,tempSubPlotIdx)
    
    plot(allTps, weightedCoord{i}, 'linewidth', p.LW)
    title_text = sprintf('CVB threshold = %d', i);
    title(title_text, 'fontsize', p.FS);
    ylabel('Pos <---> Ant', 'fontsize', p.FS)
    
    xlim([0 allTps(end)])
    tempSubPlotIdx = tempSubPlotIdx + increment;
end
xlabel('Time (unit of 10 ms)', 'fontsize', p.FS)

movingWindowSize = 5;
tempSubPlotIdx = 2;
for i = 2 : 2: numCVB
    subplot(5,increment,tempSubPlotIdx)
    
    plot(allTps, movingmean(weightedCoord{i}, 5), 'linewidth', p.LW)
    title_text = sprintf('Moving window size = %d', movingWindowSize);
    if i == 2
        title(title_text, 'fontsize', p.FS);
    end
    %     ylabel('Pos <---> Ant', 'fontsize', p.FS)
    
    xlim([0 allTps(end)])
    tempSubPlotIdx = tempSubPlotIdx + increment;
end
xlabel('Time (unit of 10 ms)', 'fontsize', p.FS)

tempSubPlotIdx = 3;
for i = 2 : 2: numCVB
    subplot(5,increment,tempSubPlotIdx)
    
    se = std(weightedCoord{i},0,2) * tinv(.975, nSubjs.ref - 1);
    plot(allTps, mean(weightedCoord{i},2), 'linewidth', p.LW)
    %     errorbar(allTps, mean(weightedCoord{i},2), se, 'linewidth', p.LW)
    title_text = sprintf('Average');
    if i == 2
        title(title_text, 'fontsize', p.FS);
    end
    %     ylabel('Pos <---> Ant', 'fontsize', p.FS)
    
    xlim([0 allTps(end)])
    tempSubPlotIdx = tempSubPlotIdx + increment;
end
xlabel('Time (unit of 10 ms)', 'fontsize', p.FS)
% legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
% suptitle('Feature locations')





