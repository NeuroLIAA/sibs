% Figure 4 - MultiMatch comparison
clc
clear all

%%
addpath('./utils/')
addpath('./MultiMatchToolbox/')

src_path = '../../data_subjects/data_final/';

aux=dir('../../matrix/images/*.mat');    filenames_img = {aux.name}'; Nimg = length(filenames_img);
aux=dir('../../matrix/subjects/*.mat');  filenames_subj = {aux.name};
clear aux

load(strcat(src_path, 'info_all_subj.mat'));
[ids_subjs, ~, subj_order]  = unique({info_per_subj_final(:).subj});
subj_order                  = subj_order'; % {info_per_subj_final.subj}
Nsubj                       = length(unique(subj_order));
Ntr                         = length(info_per_subj_final);

% save metrics
guardar = 1;

trials_tmp                      = load(strcat(src_path, 'info_all_subj.mat'));
[ids_images, ~, images_order]   = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                    = images_order';
image_size                      = trials_tmp.info_per_subj_final(1).image_size;
delta                           = 32;
grid_size                       = image_size/delta;
trials                          = reduce_scanpaths(trials_tmp.info_per_subj_final, delta, image_size);
clear trials_tmp

%% 1 - Calculate MM between subjects for each image (BH)

min_fix = 2;
max_fix = 13;
imagenes_revisar = [];
casos_revisar = {};
mean_shape_sim = nan(134,57); % para comparar con nuestro scanpath distance
index = [];
distance = [];

for ind_img=1:Nimg
    img_id = ids_images{ind_img};
    info_per_img = trials(strcmp({trials.image_name}, img_id));

    % filter 
    ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                        (length(x.fixations_matrix_reduced(:,1)) >= min_fix && ...
                                        length(x.fixations_matrix_reduced(:,1)) <= max_fix) ), info_per_img);
    % debug
    size(ind_subj_target_found)

    % calculate the distance between subjects for each image
    distance_img = [];
    index_img    = [];
    % faltan algunos sujetos saque el Nsubj
    for subj_i=1:length(info_per_img)
        %subj_i_id = ids_subjs{subj_i};
        tmp_shape = [];
        if ind_subj_target_found(subj_i)
            % if target was found and fix is in range for subj_i
            for subj_j=subj_i:length(info_per_img) % CHECK antes estaba desde subj_i en adelante, se lo saque
                if subj_j ~= subj_i
                    if ind_subj_target_found(subj_j)
                        % if target was found and fix is in range for subj_j
                        fprintf('\n img_id: %d, subj1_id: %d, subj2_id: %d\n', ind_img, subj_i, subj_j)

                        % reduced grid
                        [data1, data2] = create_trial_array_mm(subj_i, subj_j, info_per_img);
                        tmp_distance = doComparison(data1, data2, grid_size+2, 0);

                        % original size
%                         [data1, data2] = create_trial_array_mm(subj_i, subj_j, info_per_img, false);
%                         tmp_distance = doComparison(data1,data2, image_size, 0);

                        if ~isnan(tmp_distance)
                            %distance(subj_i, subj_j) = tmp_distance;
                            index = [index; subj_i subj_j ind_img];
                            distance = [distance; tmp_distance'];
                            
                            index_img = [index_img;subj_i subj_j ind_img];
                            distance_img = [distance_img; tmp_distance'];
                            
                            tmp_shape = [tmp_shape;1-tmp_distance(1)];
                        else
                            %debug
                            %keyboard
                            nan_entry = [subj_i subj_j img_id];
                            casos_revisar = [casos_revisar; nan_entry];
                        end
                    end
                end
            end
        end
        mean_shape_sim(ind_img,subj_i) = mean(tmp_shape); 
    end

    if length(distance) <= 1
        imagenes_revisar = [imagenes_revisar, ind_img];
    else
        mean_dist_img(ind_img,:)    = nanmean(distance,1)';
        std_dist_img(ind_img,:)     = nanstd(distance,0,1)';
    end    
end

%% save 

multimatch_bh = [index, distance];

if guardar
    save('results_metrics/mm_bh_reduced.mat', 'multimatch_bh', 'mean_dist_img', 'std_dist_img')
end
    
%% 2 - Calculate MM between subjects and model for each image (HM)
% Verificar siempre que concidan las variables en ambos casos si se carga
% el archivo de BH

delta       = 32;
min_fix     = 2;
max_fix     = 13;
image_size  = [768 1024];
models = fun_define_models('searchers-deepgaze');
    
% for ind_model=1:length(models)
%     models(ind_model).mean_dist_img   = nan(Nimg,Nsubj);
%     models(ind_model).std_dist_img    = nan(Nimg,Nsubj);
% 
%     models(ind_model).eliminados_subj = nan(Nimg,1);
%     models(ind_model).adentro_subj    = nan(Nimg,1);
% end

%% 2 - b

ind_model = 1;
for ind_img=1:Nimg
    img_id = ids_images{ind_img};
    info_per_img = trials(strcmp({trials.image_name}, img_id));

    % filter 
    ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                        (length(x.fixations_matrix_reduced(:,1)) >= min_fix && ...
                                        length(x.fixations_matrix_reduced(:,1)) <= max_fix) ), info_per_img);
    % load variable "cfg"
    path_cfg = sprintf('../results_models/%s/%s/%s/cfg/cfg_%d.mat',...
                    models(ind_model).prior,...
                    models(ind_model).searcher,...
                    models(ind_model).params,...
                    ind_img);
    load(path_cfg)

    % load variable "scanpath"
    path_scanpath = sprintf('../results_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                    models(ind_model).prior,...
                    models(ind_model).searcher,...
                    models(ind_model).params,...
                    ind_img);
    load(path_scanpath)
    
    % check they are the same video
    s

    if (size(scanpath,1) <= max_fix)
        distance = nan(Nsubj,5);
        for subj_i=1:length(info_per_img)
            if ind_subj_target_found(subj_i)

                fprintf('\n img_id: %d, subj1_id: %d, model: %s\n', ind_img, subj_i,  models(ind_model).name)

                sdata_model = [scanpath, 100*ones(length(scanpath),1)]; %con los filtros length deberia ser igual o mayor a 2

                sdata_subj = [info_per_img(subj_i).x_grid,...
                                info_per_img(subj_i).y_grid,...
                                100*ones(length(info_per_img(subj_i).x_grid),1)];
                
                tmp_distance = doComparison(sdata_model, sdata_subj, [768 1024], 0);

                if tmp_distance(1,1) ~= 0
                    distance(subj_i,:) = tmp_distance;
                end
            end
        end

        models(ind_model).adentro_subj(ind_img)   = sum(ind_subj_target_found);
        % los eliminados que lo encontraron pero tienen mas de
        % max_fix y menos de min_fix
        models(ind_model).eliminados_subj(ind_img)= sum(~ind_subj_target_found & [info_per_img.target_found]);

        %guardo la distancia por imagen 
        models(ind_model).mean_dist_img(ind_img,:)    = nanmean(distance)';
        models(ind_model).std_dist_img(ind_img,:)     = nanstd(distance)';
    else
        models(ind_model).mean_dist_img(ind_img,:)    = nan(1,5);
        models(ind_model).std_dist_img(ind_img,:)     = nan(1,5)';
    end
end

%% Summary

mm_mean = nanmean(models(ind_model).mean_dist_img);
mm_std  = nanstd(models(ind_model).mean_dist_img);
gold_standard = mean(mean_dist_img);