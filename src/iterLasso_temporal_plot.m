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
allTps = 1 : 3:  163

featureSelected = cell(nSubjs.ref,1);
numIterations = nan(length(allTps), nSubjs.ref);
for s = 1 : nSubjs.ref
    featureSelected{s} = cell(length(allTps),1);
end
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
        for i = 1 : numIters
            for c = 1 : numCVB
                featureSelected{s}{curTimeIdx} = union(featureSelected{s}{curTimeIdx}, results{s}.lasso{i}.voxSel{c});         
            end
        end
    end
end


%% gather data for plotting 
numFeatureSelected = nan(length(allTps), nSubjs.ref);
percentFeatureSelected = nan(length(allTps), nSubjs.ref);
legend_subjIDs = cell(nSubjs.ref,1);

for s = 1 : nSubjs.ref
    % get the legend for subject id 
    legend_subjIDs{s} = num2str(allSubjIDs.ref(s));
    for t = 1 : length(allTps)
        % feature selection measures
        numFeatureSelected(t,s) = length(featureSelected{s}{t});
        percentFeatureSelected(t,s) = numFeatureSelected(t,s)/numFeatures(s);
    end
end

%% plot 

p.FS = 14;
p.LW = 2;

% proportion of features selected over time 
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


%% number of lasso iterations over time (indicate multi-collinearity)
figure(2)


subplot(3,1,1)
plot(allTps, numIterations, 'linewidth', p.LW)
xlim([0 allTps(end)])
legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
title_text = sprintf('Feature multicollinearity over time , dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Number of lasso iterations', 'fontsize', p.FS)


subplot(3,1,2)
plot(allTps, movingmean(numIterations, 3), 'linewidth', p.LW)

legend(legend_subjIDs, 'fontsize', p.FS, 'location', 'NW')
xlim([0 allTps(end)])
title_text = sprintf('Moving average window size = 3, dataType = %s', DATA_TYPE);
title(title_text, 'fontsize', p.FS);
ylabel('Number of lasso iterations', 'fontsize', p.FS)

subplot(3,1,3)
se = tinv(.975, nSubjs.ref-1) * std(numIterations, 0, 2) / sqrt(nSubjs.ref);
hold on 
errorbar(allTps, mean(numIterations,2), se, 'linewidth', p.LW)
plot([0 allTps(end)],[0 0], '--k')
hold off
xlim([0 allTps(end)])
title('Averaged across subjects', 'fontsize', p.FS);
ylabel('Number of lasso iterations', 'fontsize', p.FS)
xlabel('Time (Unit of 10ms)', 'fontsize', p.FS)
