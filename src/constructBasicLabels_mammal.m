clear all; clc;
% specify path information

% specify parameters
DATA_TYPE = 'ref'; % 'ref' OR 'raw'
basicLabel_filename = 'Leuven_and_ECog_Items.xlsx';

DIR.PROJ = '/Users/Qihong/Dropbox/github/ECOG_Manchester/data';

DIR.LABEL = fullfile(DIR.PROJ, 'labels');
PATH.LABEL = fullfile(DIR.LABEL, basicLabel_filename); 
matadata_filename = strcat('metadata_', DATA_TYPE,'.mat');
PATH.METADATA = fullfile(DIR.PROJ, '/ECoG/data', matadata_filename);
PATH.OUT = fullfile(DIR.LABEL, 'basicLabels.mat');

%% load labels 
% sheet = 4;
% xlRange = 'C2:C51';
% [~,lab.bas.txt, ~] = xlsread(PATH.LABEL, sheet, xlRange);
% xlRange = 'D2:D51';
% [lab.bas.num,~, ~] = xlsread(PATH.LABEL, sheet, xlRange);
% xlRange = 'E2:E51';
% [~,lab.sub.txt, ~] = xlsread(PATH.LABEL, sheet, xlRange);
% 
% % save the basic level labels 
% save(PATH.OUT, 'lab')
% clear lab

%% testing 
load(PATH.OUT)
lab.bas



%% load metadata 
load(PATH.METADATA)
allObjs = metadata(1).stimuli;

% construct the label vector 
isMammal = false(length(allObjs),1);
for i = 1 : length(allObjs)
    for j = 1 : length(lab.sub.txt)
        if (strcmp(allObjs{i}, lab.sub.txt{j}) || strcmp(allObjs{i}, 'gout')) && lab.bas.num(j) == 1 
            isMammal(i) = true; 
            break; % assume non-redundency 
        end
        
    end
end

table(allObjs, isMammal)
allObjs(isMammal)
%% 
labels_allObjs.isAnimal = ~logical(metadata(1).targets(1).target);
labels_allObjs.isMammal = isMammal;
save(fullfile(DIR.LABEL, 'labels_allObjs.mat'), 'labels_allObjs')

%% 
% for a = 1 : 10
%     for b = 1 : 10 
%     if sum(metadata(a).targets(1).target == metadata(b).targets(1).target) ~= 100 
%         disp('???')
%     end
%     end
% end