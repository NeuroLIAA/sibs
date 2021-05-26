% if needed
clear all
close all

% tres structs que me modulen todo:
%   - exp que es un struc de parametros generales para no meter la pata
%   - human_data para guardar todos los calculos 


%% Load data and configurations

clc
% load paths
addpath('./utils/')
addpath('../data_analysis/utils/')
addpath('../dynamic_models/utils/')
addpath('../data_analysis/utils/heatmap_code/')
addpath('../compare_models')

% load data --> info_per_subj_final and exp struct
exp             = [];
humans          = [];
exp.src_path = '../new_data/new_matrix/';
load(strcat(exp.src_path,'info_all_subj.mat')); %varname: info_per_subj_final
[exp.subjects, ~, subj_order] = unique({info_per_subj_final(:).subj});
[exp.images_files, ~, images_order] = unique({info_per_subj_final(:).image_name});
a=dir('../matrix/images/*.mat');    filenames_img = {a.name}';
a=dir('../matrix/subjects/*.mat');  filenames_subj = {a.name};

% exp conf
exp.delta        = 32;
exp.min_fix      = 2;
exp.max_fix      = 13;
exp.image_size   = info_per_subj_final(1).image_size; % [768 1024];
exp.grid_size    = exp.image_size/exp.delta;
exp.nsacc        = [2,4,8,12];
exp.images_order = images_order';
exp.subj_order   = subj_order'; % {info_per_subj_final.subj}
exp.Nsubj        = length(unique(subj_order));
exp.Ntr          = length(info_per_subj_final);
exp.Nimg         = length(filenames_img);

% load models
models = fun_define_models_tmp(5);

%% Apply filter minimum saccades per image


%% Proportion of targets found as function of the number of saccades allowed for the humans

[humans.target_found_mean, humans.target_found, humans.target_found_ntrials, ...
    humans.target_proportion] = get_human_mean_fix_target_found(exp, info_per_subj_final);

%% Model data number of fixations predicted by the models

models = get_models_mean_fix_target_found(exp, models);


