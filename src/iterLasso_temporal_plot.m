%% Temporal analysis for the Manchster ECoG data using iterative lasso
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
% read subject info for the ref data
load(strcat(DIR.DATA{1}, 'results_ilasso_', DATA_TYPES{2}, '.mat'))
[nSubjs.ref, allSubjIDs.ref] = getSubjInfoFromResults(results);



%% collecting accuracy scores over time
% preallocate
refAccuracy = nan(nTimePts-TRIM, nSubjs.ref);

% loop over time
disp('All time points computed:')
allTps = 1 : 3:  163

featureSelected = cell(nSubjs.ref,1);
numIterations = nan(length(allTps), nSubjs.ref);
for s = 1 : nSubjs.ref
    featureSelected{s} = cell(length(allTps),1);
end

wtsMagnitude_t_subj = nan(length(allTps), nSubjs.ref);
wtsMagnitude_nz_t_subj = nan(length(allTps), nSubjs.ref);
for t = allTps
    curTimeIdx = find(allTps==t);
    
    % process REF data
    DATA_TYPE = DATA_TYPES{2};
    % load metadata
    load(strcat(DIR.METADATA{t}, 'metadata_', DATA_TYPE, '.mat'))
    % loop over subjects
    load(strcat(DIR.DATA{t}, 'results_ilasso_', DATA_TYPE, '.mat'))
    
    for s = 1 : nSubjs.ref
        % total number of columns (features)
        numFeatures(s) = metadata(allSubjIDs.ref(s)).ncol;
        % number of iterations
        numIters = length(results{s}.lasso);
        numIterations(curTimeIdx,s) = numIters;
        
        wtsMagnitude_overIter = 0;
        wtsMagnitude_nz_overIter = 0;
        for i = 1 : numIters
            allcoef_acrossCV = [];
            for c = 1 : numCVB
                % record feature selected, collapse over CVB
                featureSelected{s}{curTimeIdx} = union(featureSelected{s}{curTimeIdx}, results{s}.lasso{i}.voxSel{c});
                
                % get weights over all cross validation
                wts = results{s}.lasso{i}.coef_min{c};
                wts(1) = [];    % remove the intercept
                % collect weight vector across CVB
                allcoef_acrossCV = vertcat(allcoef_acrossCV, results{s}.lasso{i}.coef_min{c});
                allcoef_acrossCV_nz = allcoef_acrossCV(allcoef_acrossCV ~= 0);
            end
            % accumulate average weight magnitude across iteration
            wtsMagnitude_overIter = wtsMagnitude_overIter + mean(abs(allcoef_acrossCV));
            wtsMagnitude_nz_overIter = wtsMagnitude_nz_overIter + mean(abs(allcoef_acrossCV_nz));
        end
        % put weight magnitude into a matrix
        wtsMagnitude_t_subj(curTimeIdx, s) = wtsMagnitude_overIter;
        wtsMagnitude_nz_t_subj(curTimeIdx, s) = wtsMagnitude_nz_overIter;
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
numFeatureSelected = nan(length(allTps), nSubjs.ref);
percentFeatureSelected = nan(length(allTps), nSubjs.ref);
legend_subjIDs = cell(nSubjs.ref,1);
numSensors_bySubj = zeros(nSubjs.ref,1);

weightedCoord = cell(numCVB,1);
% weightedCoord = nan(length(allTps),nSubjs.ref);
for s = 1 : nSubjs.ref
    curSubjID = allSubjIDs.ref(s);
    % get the legend for subject id
    legend_subjIDs{s} = num2str(curSubjID);
    for t = 1 : length(allTps)
        % feature selection measures
        numFeatureSelected(t,s) = length(featureSelected{s}{t});
        percentFeatureSelected(t,s) = numFeatureSelected(t,s)/numFeatures(s);
        
        %% electode y dimention movement analysis
        numSensors_bySubj = length(metadata(curSubjID).coords.labels);
        ycoords = metadata(curSubjID).coords.xyz(:,2) - yCoords.median;
        % preallocate: number of times each electrode got selected
        electrodeWeights = zeros(numSensors_bySubj,1);
        % find out the start and end idx for each electrode window
        idx_t0 = metadata(curSubjID).WindowSizeInMilliseconds * (0 : numSensors_bySubj-1);
        idx_t50 = metadata(curSubjID).WindowSizeInMilliseconds * (1 : numSensors_bySubj);
        for e = 1 : numSensors_bySubj
            electrodeWeights(e) = sum(featureSelected{s}{t} > idx_t0(e) & featureSelected{s}{t}<= idx_t50(e));
        end
        % consider
        for i = 1 : numCVB
            thres_electrodeWeights = electrodeWeights >= i;
            weightedCoord{i}(t,s) = mean(ycoords .* thres_electrodeWeights);
        end
    end
end

%% plot

p.FS = 14;
p.LW = 2;


%% proportion of features selected over time
figure(1)

subplot(3,1,1)
plot(allTps, percentFeatureSelected, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Proportion features selected over time, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('% features selected', 'fontsize', p.FS)


subplot(3,1,2)
plot(allTps, movingmean(percentFeatureSelected, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('% features selected', 'fontsize', p.FS)


subplot(3,1,3)
se = tinv(.975, nSubjs.ref-1) * std(percentFeatureSelected, 0, 2) / sqrt(nSubjs.ref);
hold on
errorbar(allTps, mean(percentFeatureSelected,2), se, 'linewidth', p.LW)
plot([0 allTps(end)],[0 0], '--k')
hold off
xlim([0 allTps(end)])
title('Averaged across subjects', 'fontsize', p.FS);
ylabel('% features selected', 'fontsize', p.FS)
xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)


% %% number of lasso iterations over time (indicate multi-collinearity)
% figure(2)
%
%
% subplot(3,1,1)
% plot(allTps, numIterations, 'linewidth', p.LW)
% xlim([0 allTps(end)])
% legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
% title_text = sprintf('Feature multicollinearity over time , dataType = %s', DATA_TYPE);
% title(title_text, 'fontsize', p.FS);
% ylabel('Number of lasso iterations', 'fontsize', p.FS)
%
%
% subplot(3,1,2)
% plot(allTps, movingmean(numIterations, 3), 'linewidth', p.LW)
%
% legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
% xlim([0 allTps(end)])
% title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
% title(title_text, 'fontsize', p.FS);
% ylabel('Number of lasso iterations', 'fontsize', p.FS)
%
% subplot(3,1,3)
% se = tinv(.975, nSubjs.ref-1) * std(numIterations, 0, 2) / sqrt(nSubjs.ref);
% hold on
% errorbar(allTps, mean(numIterations,2), se, 'linewidth', p.LW)
% plot([0 allTps(end)],[0 0], '--k')
% hold off
% xlim([0 allTps(end)])
% title('Averaged across subjects', 'fontsize', p.FS);
% ylabel('Number of lasso iterations', 'fontsize', p.FS)
% xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)



%% average weight magnitude
figure(3)

subplot(3,2,1)
plot(allTps, wtsMagnitude_t_subj, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Weights magnitude over time, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of all features', 'fontsize', p.FS)


subplot(3,2,3)
plot(allTps, movingmean(wtsMagnitude_t_subj, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of all features', 'fontsize', p.FS)

subplot(3,2,5)
se = tinv(.975, nSubjs.ref-1) * std(wtsMagnitude_t_subj, 0, 2) / sqrt(nSubjs.ref);
hold on
errorbar(allTps, mean(wtsMagnitude_t_subj,2), se, 'linewidth', p.LW)
plot([0 allTps(end)],[0 0], '--k')
hold off
xlim([0 allTps(end)])
title('Averaged across subjects', 'fontsize', p.FS);
ylabel('Average weights of all features', 'fontsize', p.FS)
xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)

subplot(3,2,2)
plot(allTps, wtsMagnitude_nz_t_subj, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Weights magnitude over time, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of non-zero features', 'fontsize', p.FS)


subplot(3,2,4)
plot(allTps, movingmean(wtsMagnitude_nz_t_subj, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Average weights of non-zero features', 'fontsize', p.FS)

subplot(3,2,6)
se = tinv(.975, nSubjs.ref-1) * std(wtsMagnitude_nz_t_subj, 0, 2) / sqrt(nSubjs.ref);
hold on
errorbar(allTps, mean(wtsMagnitude_nz_t_subj,2), se, 'linewidth', p.LW)
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





