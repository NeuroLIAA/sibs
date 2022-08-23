% Auxiliar - MultiMatch calculations
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
guardar = 0;

trials_tmp                    = load(strcat(src_path, 'info_all_subj.mat'));
[ids_images, ~, images_order] = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                  = images_order';
image_size                    = trials_tmp.info_per_subj_final(1).image_size;
delta                         = 32;
grid_size                     = image_size/delta;
trials                        = reduce_scanpaths(trials_tmp.info_per_subj_final, delta, image_size);
multimatch_names              = {'vectorSim','directionSim','lengthSim','positionSim','durationSim'};
clear trials_tmp

%% 1 - Calculate MM between subjects for each image (BH)

min_fix = 2;
max_fix = 13;
imagenes_revisar = [];
casos_revisar = {};
mean_shape_sim = nan(134,57); % para comparar con nuestro scanpath distance
index = [];
distance = [];
adentro_subj = nan(1,Nimg);
eliminados_subj = nan(1,Nimg);

SubjectsData = [];
row = 1;
for ind_img=1:Nimg
    img_id = ids_images{ind_img};
    info_per_img = trials(strcmp({trials.image_name}, img_id));
    %data = [info_per_img, data];
    % filter 
    ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                        (length(x.fixations_matrix_reduced(:,1)) >= min_fix && ...
                                        length(x.fixations_matrix_reduced(:,1)) <= max_fix) ), info_per_img);
    % debug
    size(ind_subj_target_found);    

    % calculate the distance between subjects for each image
    distance_img = [];
    index_img    = [];
    % faltan algunos sujetos saque el Nsubj
    for subj_i=1:length(info_per_img)
        if ind_subj_target_found(subj_i)
            % if target was found and fix is in range for subj_i
            [data1, data2] = create_trial_array_mm(subj_i, subj_i, info_per_img);
            SubjectsData(row).img_id = info_per_img(subj_i).image_name;
            SubjectsData(row).subj_id = info_per_img(subj_i).subj;
            SubjectsData(row).scanpath_x = data1(:,2);
            SubjectsData(row).scanpath_y = data1(:,1);
            SubjectsData(row).target_found = true;
            SubjectsData(row).grid_pad_y = 26;
            SubjectsData(row).grid_pad_x = 34;            
            row = row + 1;
        end
    end    
end

%%

fieldnames(SubjectsData)
%%
%saveJSONfile(SubjectsData,'SubjectsData.json')
saveAsJSON(SubjectsData,'SubjectsData.json')
%%
d = jsonencode(SubjectsData);
fid = fopen('ss.json','wt');
fprintf(fid, d);
fclose(fid);