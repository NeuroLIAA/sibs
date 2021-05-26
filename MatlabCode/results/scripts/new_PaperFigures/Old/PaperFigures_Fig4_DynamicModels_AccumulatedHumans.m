clear all
close all
clc

% addpath('~/Dropbox/my_functions/')
addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
% load('../matrix/info_per_subj_additional_data.mat'); % loads images, img_order, templates, template_order
% addpath(genpath('../static_models/utils/'));
% 
% load('../matrix/initial_fixations.mat');

load('../matrix/info_per_subj_final.mat');
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);

a=dir('../matrix/images/*.mat');    filenames_img = {a.name}'; Nimg = length(filenames_img);
a=dir('../matrix/subjects/*.mat');  filenames_subj = {a.name};

% models = fun_define_models;
models = fun_define_models_tmp(3);

%% Human data
nntrthr     = 20;

% Average number of fixations needed by the observers to find the target
Nfix_img_mean   = nan(Nimg,1);
Nfix_img_std    = nan(Nimg,1);
Nfix_img_nsuj   = nan(Nimg,1);
for ind_img=1:Nimg % images
    path = char(strcat('../matrix/images/info_per_subj_img_', num2str(ind_img), '.mat')); fprintf('%s\n',path);
    load(path);
    if (ind_img ~= 132)
        founds = [info_per_subj.target_found];
        fixs = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_subj);
        Nfix_img_mean(ind_img)  = mean(fixs(founds));
        Nfix_img_std(ind_img)   = std(fixs(founds));
        Nfix_img_nsuj(ind_img)  = sum(founds);
    end
end
nntrfilt    = (Nfix_img_nsuj>nntrthr);

% Proportion of targets found as function of the number of saccades allowed for the humans
target_found_img        = nan(Nsubj,4);
target_found_img_nsuj   = nan(Nsubj,4);
for ind_subj=1:Nsubj
    path = char(strcat('../matrix/subjects/info_per_subj_subj_', num2str(ind_subj), '.mat')); fprintf('%s\n',path);
    load(path);
    
    founds              = [info_per_subj.target_found];
    nsaccades_allowed   = arrayfun(@(x) x.exp_data.nsaccades_thr, info_per_subj);
    NSACCADES_allowed   = [2 4 8 12]; %unique(nsaccades_allowed);
    for isacc=1:length(NSACCADES_allowed)
        target_found_img(ind_subj,isacc)        = sum(founds(nsaccades_allowed==NSACCADES_allowed(isacc)));
        target_found_img_nsuj(ind_subj,isacc)   = sum(nsaccades_allowed==NSACCADES_allowed(isacc));
    end
    if size(target_found_img,2)>4; break; end
end
P_target_found = target_found_img./target_found_img_nsuj;
