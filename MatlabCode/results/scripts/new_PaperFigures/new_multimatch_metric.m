clc
clear all

% Script para probar métrica MultiMatch
% 21/1/21: Utilizar el doComparison para calcular todas las métricas. Para
% eso hay que adaptar al formato que necesita multimatch
% 3/3/21: Se agregaron variables para guardar y graficar y se chequeo que
% este funcionando correctamente hasta 2:
%   TODO: Verificar  
%%
addpath('./utils/')
addpath('./MultiMatchToolbox/')
addpath('../data_analysis/utils/')
addpath('../dynamic_models/utils/')
addpath('../data_analysis/utils/heatmap_code/')
addpath('../compare_models')

src_path = '../new_data/new_matrix/';

aux=dir('../matrix/images/*.mat');    filenames_img = {aux.name}'; Nimg = length(filenames_img);
aux=dir('../matrix/subjects/*.mat');  filenames_subj = {aux.name};
clear aux

load(strcat(src_path,'info_all_subj.mat'));
[ids_subjs, ~, subj_order]   = unique({info_per_subj_final(:).subj});
subj_order                  = subj_order'; % {info_per_subj_final.subj}
Nsubj                       = length(unique(subj_order));
Ntr                         = length(info_per_subj_final);

%debug
guardar = 1;
graficar = 1;

trials_tmp                      = load(strcat(src_path,'info_all_subj.mat'));
[ids_images, ~, images_order]   = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                    = images_order';
image_size                      = trials_tmp.info_per_subj_final(1).image_size;
delta                           = 32;
grid_size                       = image_size/delta;
trials                          = new_subjMapFixationToMatrix(trials_tmp.info_per_subj_final, '', delta, image_size);
clear trials_tmp

%% 1 - Calculate MM between subjects for each image (BH)

min_fix = 2;
max_fix = 13;
imagenes_revisar = [];
casos_revisar = {};
mean_shape_sim = nan(134,57); % para comparar con nuestro scanpath distance
for ind_img=1:Nimg
    if ind_img ~= 132 % lo haría para todas las imagenes
        %path = char(strcat(src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat')); 
        %fprintf('%s\n',path);
        %load(path);
        
        img_id = ids_images{ind_img};
        info_per_img = trials(strcmp({trials.image_name}, img_id));
        
        % filter 
        ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                            (length(x.fixations_matrix_reduced(:,1)) > min_fix && ...
                                            length(x.fixations_matrix_reduced(:,1)) <= max_fix) ), info_per_img);
        % debug
        size(ind_subj_target_found)
        
        % calculate the distance between subjects for each image
        %dist_table  = cell2table(cell(0,4), 'VariableNames', {'Subj1', 'Subj2', 'ImageInd', 'mm'});
        distance = [];
        index = [];
        % faltan algunos sujetos saque el Nsubj
        for subj_i=1:length(info_per_img)
            %subj_i_id = ids_subjs{subj_i};
            tmp_shape = [];
            if ind_subj_target_found(subj_i)
                % if target was found and fix is in range for subj_i
                for subj_j=1:length(info_per_img) % CHECK antes estaba desde subj_i en adelante, se lo saque
                    %subj_j_id = ids_subjs{subj_j};
                    if subj_j ~= subj_i
                        if ind_subj_target_found(subj_j)
                            % if target was found and fix is in range for subj_j
                            fprintf('\n img_id: %d, subj1_id: %d, subj2_id: %d\n', ind_img, subj_i, subj_j)
                            
                            %tmp_distance = scanpathDistance(info_per_img(subj_i).fixations_matrix_reduced,...
                            %                                            info_per_img(subj_j).fixations_matrix_reduced,...
                            %                                            grid_size);
                            
                            
                            %[data1, data2] = createTrialArrayMM(subj_i, subj_j, img_id, info_per_img);
                            
                            [data1, data2] = createTrialArrayMM(subj_i, subj_j, img_id, info_per_img, false);
                            
                            % FIX poner un try y agarrar ambos problemas
                            % para ya solucionar, 1 - menos de una fix, 2 -
                            % no coincidencia de los tamaños entre las
                            % grillas reducidas y el scanpath original
                            
                            %tmp_distance = doComparison(data1,data2, [25 33], 0);
                            tmp_distance = doComparison(data1,data2, [1024,768], 0);
                            
                            
                            if ~isnan(tmp_distance)
                                %distance(subj_i, subj_j) = tmp_distance;
                                index = [index; subj_i subj_j ind_img];
                                distance = [distance; tmp_distance'];
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
end

%% save 
if guardar
    save('log/subsj_distance_reduced.mat', 'mean_dist_img', 'std_dist_img', 'mean_shape_sim')
end
    
%% 2 - Calculate MM between subjects and model for each image (HM)
% Verificar siempre que concidan las variables en ambos casos si se carga
% el archivo de BH

delta       = 32;
min_fix     = 2;
max_fix     = 13;
image_size  = [768 1024];
% priors case = 6, searchers case = 5
models = fun_define_models_tmp(6);
    
% for ind_model=1:length(models)
%     models(ind_model).mean_dist_img   = nan(Nimg,Nsubj);
%     models(ind_model).std_dist_img    = nan(Nimg,Nsubj);
% 
%     models(ind_model).eliminados_subj = nan(Nimg,1);
%     models(ind_model).adentro_subj    = nan(Nimg,1);
% end

%% 2 bis

ind_model = 1;
for ind_img=1:Nimg
    if ind_img ~= 132
        path = char(strcat(src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat'));
        fprintf('%s\n',path);
        load(path)

        % (2019-09-18) Map fixations into the Grid
        % Update fixations_matrix_reduced, timefix_matrix_reduced, mascara_matrix_reduced and delta
        % fields in the info_per_subj structure.
        info_per_img = new_subjMapFixationToMatrix( info_per_img, path, delta, image_size );
        % Calculo la distancia de todos contra todos para los sujetos que
        % tengan min_fix < #fix < max_fix y que ambos hayan encontrado el
        % target 

        % Calculate the distance between subjects for each image
        ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                            (length(x.fixations_matrix_reduced(:,1)) > min_fix && ...
                                            length(x.fixations_matrix_reduced(:,1)) <= max_fix) ), info_per_img);
    %indice aclarado mas arriba
    %for ind_model=1:length(models) 
        
        path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                        models(ind_model).prior,...
                        models(ind_model).searcher,...
                        models(ind_model).params,...
                        ind_img);
        load(path)

        if (size(scanpath,1) < max_fix)
            distance = nan(Nsubj,5);
            for subj_i=1:length(info_per_img)
                if ind_subj_target_found(subj_i)
                    %try
                        fprintf('\n img_id: %d, subj1_id: %d, model: %s\n', ind_img, subj_i,  models(ind_model).name)
%                         sdata_model = [scanpath, 100*ones(length(scanpath),1)]; %con los filtros length deberia ser igual o mayor a 2
%                         sdata_subj = [info_per_img(subj_i).fixations_matrix_reduced(:,2),...
%                                         info_per_img(subj_i).fixations_matrix_reduced(:,1),...
%                                         100*ones(length(info_per_img(subj_i).fixations_matrix_reduced(:,2)),1)];
%                         tmp_distance = doComparison(sdata_model, sdata_subj, [25 33], 0);

                        sdata_model = [scanpath, 100*ones(length(scanpath),1)]; %con los filtros length deberia ser igual o mayor a 2
                        
                        sdata_subj = [info_per_img(subj_i).x',...
                                        info_per_img(subj_i).y',...
                                        100*ones(length(info_per_img(subj_i).fixations_matrix_reduced(:,2)),1)];
                        
                        tmp_distance = doComparison(sdata_model, sdata_subj, [768 1024], 0);
                        
                        if tmp_distance(1,1) ~= 0
                            distance(subj_i,:) = tmp_distance;
                        end
                            
%                         distance(subj_i) = scanpathDistance(info_per_img(subj_i).fixations_matrix_reduced,...
%                                                                 [scanpath(:,2) scanpath(:,1)],...
%                                                                 grid_size);
                    %catch e
                    %    fprintf(2,"Identificador de error:\n%s \n",e.identifier);
                    %    fprintf(1,"Mensaje del error:\n%s \n", e.message);
                    %    keyboard
                        %break
                    %end
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
end

%% Summary

mm_mean = nanmean(models(ind_model).mean_dist_img);
mm_std  = nanstd(models(ind_model).mean_dist_img);
gold_standard = mean(mean_dist_img);

%% Plot
metric_ind = 1;
x = 1-models(ind_model).mean_dist_img(:,metric_ind);
y = 1-mean_dist_img(:,metric_ind);
not_nan = isfinite(x) & isfinite(y); 
figure(1)
    scatter(x,y,[], models(ind_model).cols, 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0 0.25])
    ylim([0 0.25])
    refline(1,0)

%% Comparacion MM vs SD

load('humans_data.mat') %var humans
x_MM = nanmean(mean_shape_sim,2);
y_SD = nanmean(humans.mean_dist_img,2);
not_nan = isfinite(x_MM) & isfinite(y_SD); 
figure(2)
    scatter(x_MM,y_SD,[], 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0 0.25])
    ylim([0 0.25])
    refline(1,0)
    xlabel('MultiMatch Shape')
    ylabel('Scanpath Dissimilarity')
    title('Between Humans Metrics')
    
%%
figure(3)
    boxplot([x_MM;y_SD],[ones(size(x_MM));2*ones(size(y_SD))]);
    set(gca,'XTickLabel', {'MM Shape','SD'})

%%
x_MM_sorted = sort(x_MM);
y_SD_sorted = sort(y_SD);
not_nan = isfinite(x_MM) & isfinite(y_SD); 
figure(4)
    scatter(x_MM_sorted, y_SD_sorted,[], 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0 0.25])
    ylim([0 0.25])
    refline(1,0)
    
%%
figure(5)
    boxplot([1-mean_dist_img(:,1)'; 1-mean_dist_img_singrilla(:,1)'],[ones(size(mean_dist_img(:,1)')),2*ones(size(mean_dist_img(:,1)'))])
    ylim([0 0.25])
    set(gca,'XTickLabel', {'MMs Sin Grilla','MMs Con Grilla'})
