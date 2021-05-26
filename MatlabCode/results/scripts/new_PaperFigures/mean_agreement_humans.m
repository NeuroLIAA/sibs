
clear all
close all
clc
%%
addpath('./utils/')
addpath('../data_analysis/utils/')
addpath('../dynamic_models/utils/')
addpath('../data_analysis/utils/heatmap_code/')
addpath('../compare_models')

src_path = '../new_data/new_matrix/';

load(strcat(src_path,'info_all_subj.mat'));
[subjects, ~, subj_order] = unique({info_per_subj_final(:).subj});
[images_files, ~, images_order] = unique({info_per_subj_final(:).image_name});
subj_order  = subj_order'; % {info_per_subj_final.subj}
images_order = images_order';
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj_final);

a=dir('../matrix/images/*.mat');    filenames_img = {a.name}'; Nimg = length(filenames_img);
a=dir('../matrix/subjects/*.mat');  filenames_subj = {a.name};

% Para variar priors case = 6
% Para variar serchers case = 5
models = fun_define_models_tmp(5);
delta = 32;
image_size = info_per_subj_final(1).image_size;

%% creacion de un indice de imagenes
% asumimos que las imagenes son las que vio el sujeto 1
% aux_path = char(strcat(src_path,'sinfo_subj/info_per_subj_1.mat'));
% load(aux_path)
% 
% 
% %%
% del info_per_subj aux_path
%% armado del tensor booleano
nsacc = [2,4,8,12];
bool_subj_trials = nan(Nsubj, Nimg, 4);
for ind_trial=1:length(info_per_subj_final)
    ind_subj = subj_order(ind_trial);
    ind_img = images_order(ind_trial);
    ind_nsacc = find(nsacc==info_per_subj_final(ind_trial).nsaccades_allowed);
    if info_per_subj_final(ind_trial).nsaccades_allowed==16 || ...
            info_per_subj_final(ind_trial).nsaccades_allowed==64
        ind_nsacc = 4;
    end
    if info_per_subj_final(ind_trial).nsaccades_allowed==3
        ind_nsacc = 2;
    end
    if length(ind_nsacc) ~= 1
        keyboard
    end
    if info_per_subj_final(ind_trial).target_found
        bool_subj_trials(ind_subj, ind_img, ind_nsacc) = 1;
    else
        bool_subj_trials(ind_subj, ind_img, ind_nsacc) = 0;
    end
end

%% calculo de distancias
subj_agreement = nan(Nsubj, Nsubj);
for subj_idx = 1:Nsubj
    for subj_idy = subj_idx+1:Nsubj
        subj1 = squeeze(bool_subj_trials(subj_idx,:,:));
        subj2 = squeeze(bool_subj_trials(subj_idy,:,:));
        diff  = abs(subj1-subj2);
        subj_agreement(subj_idx,subj_idy) = 1-nanmean(diff,'all');
    end
end
notnan_idx = isfinite(subj_agreement);
mean_agreement = mean(subj_agreement(notnan_idx));
std_agreement = std(subj_agreement(notnan_idx));

%% visualización

hist(subj_agreement(notnan_idx))
hold on
xline(mean_agreement,'-r')
%% chequeo

for i=1:Nsubj
   path = char(strcat(src_path,'sinfo_subj/info_per_subj_', num2str(i), '.mat'));
   load(path);
   if length(info_per_subj) ~=134
       keyboard
   end
end