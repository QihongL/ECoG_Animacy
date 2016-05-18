%% Do some quick analysis for the Manchster ECoG data
clear variables; clc; 
% specify path information

% point to the output directory
DIR.OUT = '/Users/Qihong/Dropbox/github/ECOG_Manchester/results/';
DATA_TYPE = 'raw';

% load data 
load(strcat(DIR.OUT, 'results_', DATA_TYPE, '.mat'))
nSubj = length(results);

subjIDs = nan(nSubj,1);
accuracies.lasso.onese = nan(nSubj,1);
accuracies.lasso.min = nan(nSubj,1);
accuracies.ridge.onese = nan(nSubj,1);
accuracies.ridge.min = nan(nSubj,1);
for i = 1 : nSubj
    subjIDs(i) = results{i}.subjID;
    % mean accuracy over cvb 
    accuracies.lasso.onese(i) = mean(results{i}.lasso.accuracy.onese);
    accuracies.lasso.min(i) = mean(results{i}.lasso.accuracy.min);
    accuracies.ridge.onese(i) = mean(results{i}.ridge.accuracy.onese);
    accuracies.ridge.min(i) = mean(results{i}.ridge.accuracy.min);
end

%%
batMat = horzcat(accuracies.lasso.onese, accuracies.lasso.min, ...
    accuracies.ridge.onese, accuracies.ridge.min);

bar(batMat)

ax = gca;
ax.XTickLabel = subjIDs;

FS = 14;
title_text = sprintf('Accuracy, all time pts, data: %s',DATA_TYPE);
title(title_text, 'fontsize' , FS)
ylabel('Classifcation accuracy', 'fontsize' , FS)
xlabel('Subject ID', 'fontsize' , FS)
leg = legend({'lasso(1se)','lasso(min)','ridge(1se)','ridge(min)'}, 'location', 'northwest');
set(leg,'FontSize',FS);