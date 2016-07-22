%% Do some quick analysis for the Manchster ECoG data
% a single iteration of logistic regression with lasso or ridge penalty
clear variables; clc; clf;

% specify path information
% point to the output directory
DIR.OUT = '/Users/Qihong/Dropbox/github/ECOG_Manchester/results/allTimePts/';


% specify parameters
DATA_TYPE = 'ref'; % 'ref' OR 'raw'
resultsFileName = strcat('results_ilasso_', DATA_TYPE,'_bc010_wStart0200wSize1000');

% constant
chance = .5;
alpha = .05;
numCVB = 10;
p.FS = 14;
p.LW = 2;

load(fullfile(DIR.OUT, resultsFileName))

% loop over subjects
for s = 1 : length(results)
    
    accuracies = nan(results{s}.iterNum,numCVB);
    numVoxSelected = nan(results{s}.iterNum,numCVB);
    % loop over lasso iterations
    for i = 1 : results{s}.iterNum
        accuracies(i,:) = results{s}.lasso{i}.accuracy;
        for c = 1 : numCVB
            numVoxSelected(i,c) = length(results{s}.lasso{i}.voxSel{c});
        end
    end
    
    figure(1)
    subplot(2, ceil(length(results)/2), s)
    se = tinv(1 - alpha/2, numCVB - 1) * std(accuracies,0,2) / sqrt(numCVB);
    errorbar(1:results{s}.iterNum, mean(accuracies,2),se, 'linewidth', p.LW);
    title_text = sprintf('Subject ID = %d, dataType = %s', results{s}.subjID, DATA_TYPE);
    title(title_text, 'fontsize', p.FS)
    ylabel('Mean CV accuracy', 'fontsize', p.FS)
    xlabel('Lasso iterations', 'fontsize', p.FS)
    xlim([1, results{s}.iterNum])
    ylim([.5 1])
    
    figure(2)
    subplot(2, ceil(length(results)/2), s)
    se = tinv(1 - alpha/2, numCVB - 1) * std(numVoxSelected,0,2) / sqrt(numCVB);
    errorbar(1:results{s}.iterNum, mean(numVoxSelected,2),se, 'linewidth', p.LW);
    title_text = sprintf('Subject ID = %d, dataType = %s', results{s}.subjID, DATA_TYPE);
    title(title_text, 'fontsize', p.FS)
    ylabel('Mean numVoxel selected', 'fontsize', p.FS)
    xlabel('Lasso iterations', 'fontsize', p.FS)
    xlim([1, results{s}.iterNum])
    
end

