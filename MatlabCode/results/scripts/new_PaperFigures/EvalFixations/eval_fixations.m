clc
clear

%%
addpath('../utils/')
addpath('../utils/saliency_mit/code_forMetrics/')
addpath('../../data_analysis/utils/')
addpath('../../dynamic_models/utils/')
addpath('../../data_analysis/utils/heatmap_code/')
image_ind  = 2;
src_path   = '../../new_data/new_matrix/';
model      = 'correlation';
prior      = 'deepgaze';
delta      = 32;
model_path = char(strcat('../../out_models/', prior, '/', model, '/a_3_b_4_tam_celda_32/probability_map/'));

trials_tmp                      = load(strcat(src_path,'info_all_subj.mat'));
[images_files, ~, images_order] = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                    = images_order';
image_size = trials_tmp.info_per_subj_final(1).image_size;
trials     = new_subjMapFixationToMatrix(trials_tmp.info_per_subj_final, '', delta, image_size);
clear trials_tmp

%%
nfix = [2,3,4,5,6,7,8,9];
metrics = nan(3,length(nfix));
tot_subj = nan(1,length(nfix));
for ind_fix = 1:length(nfix) %length fixations
    % para esa imagen y esa nfix, cargo el mapa de prioridad del modelo
    [fixations_map, fixations, tot_subj(ind_fix)] = loadProbabilityMap(images_files(image_ind), image_ind, nfix(ind_fix), trials);
    % para esa imagen y esa nfix, cargo todas esas fijaciones
    [probability_map, ~] =  loadProbabilityMap(images_files(image_ind), image_ind, nfix(ind_fix), trials, model_path);
    heatmap(probability_map)
    %keyboard
    % calculate AUC Judd, NSS and KLdiv
    metrics (1,ind_fix) = AUC_Judd(probability_map, fixations_map); % second coord of metrics should be nfix
    metrics (2,ind_fix) = NSS(probability_map, fixations_map);
    metrics (3,ind_fix) = KLdiv(probability_map, fixations_map);
end

%%
%resultados = array2table((nfix-1)',metrics', tot_subj');