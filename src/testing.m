clear all; clc;
% specify path information
datadir = '/Users/Qihong/Dropbox/github/ECOG_Manchester/data/labels';
basicLabelFileName = 'Leuven_and_ECog_Items.xlsx';

x = xlsread(fullfile(datadir, basicLabelFileName));


animals = metadata(1).stimuli(~logical(metadata(1).targets(1).target));

