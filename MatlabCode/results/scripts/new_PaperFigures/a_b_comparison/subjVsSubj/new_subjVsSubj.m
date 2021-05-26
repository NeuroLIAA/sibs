clear all
close all
clc
% para actualizar con la misma version de los datos, por las dudas
% 23/2/21: actualizar para que use los datos curados nuevos
% 26/2/21: actualizado, genera un archivo en la carpeta con el max_fix

%% En cada imagen, por cada persona calcula la distancia promedio de su scanpath al resto de las personas. 

delta = 32;
min_fix = 3;
max_fix = 9;
image_size = [768 1024];
addpath('../../../compare_models')
addpath('../../../dynamic_models/utils')
addpath('../../utils')
    
mean_dist_img = [];
std_dist_img = [];

eliminados_subj = [];
adentro_subj = [];
    
for img=1:134
    if img ~= 132
        path = char(strcat('../../../new_data/new_matrix/sinfo_img/info_per_img_', num2str(img), '.mat')); fprintf('%s\n',path);
        load(path)
        
        elim = 0;
        aden = 0;
        
        info_per_subj = new_subjMapFixationToMatrix( info_per_img, path, delta, image_size );
        % Calculo la distancia de todos contra todos para los sujetos que
        % tengan min_fix < #fix < max_fix y que ambos hayan encontrado el
        % target 
        
        distance = nan(57,57);
        for subj_i=1:length(info_per_subj)
            
            for k=1:length(info_per_subj(subj_i).x)
                
                if info_per_subj(subj_i).target_found && ...
                        length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) >= min_fix && ...
                        length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) <= max_fix

                    for subj_j=(subj_i+1):length(info_per_subj) 
                        if info_per_subj(subj_j).target_found && ...
                                length(info_per_subj(subj_j).fixations_matrix_reduced(:,1)) >= min_fix && ...
                                length(info_per_subj(subj_j).fixations_matrix_reduced(:,1)) <= max_fix

                            distance(subj_i, subj_j) = scanpathDistance(info_per_subj(subj_i).fixations_matrix_reduced,info_per_subj(subj_j).fixations_matrix_reduced);
                        end
                    end
                end
            end
        end
        
        for subj_i=1:length(info_per_subj)
            if info_per_subj(subj_i).target_found && ...
                    length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) >= min_fix && ...
                    length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) <= max_fix
                aden = aden + 1;
            else
                if info_per_subj(subj_i).target_found
                    elim = elim + 1;
                end
            end
        end
        adentro_subj = [adentro_subj aden];
        eliminados_subj = [eliminados_subj elim];
        
        mean_dist = [];
        std_dist = [];
        
        % Calculo la distancia promedio de todos contra todos
        full_distance = [distance distance'];
        for i=1:57
            mean_dist = [mean_dist nanmean(full_distance(i,:))];
            std_dist = [std_dist nanstd(full_distance(i,:))];
        end
        
        %guardo la distancia por imagen 
        mean_dist_img = [mean_dist_img; mean_dist];
        std_dist_img = [std_dist_img; std_dist];
    end
end

save_name = char(strcat('subjVsSubj_delta_',num2str(delta),'_max_fix_',num2str(max_fix)));
save(save_name, 'adentro_subj', 'eliminados_subj', 'delta', 'max_fix', 'mean_dist_img', 'min_fix', 'std_dist_img')

%% Cantidad de sujetos que encuentran el objeto con X fijaciones
delta = 32;
min_fix = 0;
max_fix = 18;

count_subj = zeros(1,18);

for img=1:134
    if img ~= 132
        path = char(strcat('../../../new_data/new_matrix/sinfo_img/info_per_img_', num2str(img), '.mat')); fprintf('%s\n',path);
        load(path)
        
        for subj_i=1:length(info_per_subj)
            if info_per_subj(subj_i).target_found && ...
                    length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) > min_fix && ...
                    length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) <= max_fix
                
                count_fix = length(info_per_subj(subj_i).fixations_matrix_reduced(:,1));
                
                count_subj(1,count_fix) = count_subj(1,count_fix) + 1;
                
            end
        end
    end
end

%% Cantidad de imágenes encontradas para cada fijación 
delta = 32;
min_fix = 0;
max_fix = 18;

count_img = zeros(134,18);

for img=1:134
    if img ~= 132
        path = char(strcat('../../../new_data/new_matrix/sinfo_img/info_per_img_', num2str(img), '.mat')); fprintf('%s\n',path);
        load(path)
        
        for subj_i=1:length(info_per_subj)
            if info_per_subj(subj_i).target_found && ...
                    length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) > min_fix && ...
                    length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) <= max_fix
                
                fix = length(info_per_subj(subj_i).fixations_matrix_reduced(:,1));
                count_img(img,fix) = 1;
                
            end
        end
    end
end

count_img = sum(count_img);

%% Plot Cantidad de sujetos que encuentran el objeto con X fijaciones

figure;
subplot(2,1,1)
    bar(count_subj)
    set(gca,'XLim',[1 18])
    title('Cantidad de sujetos que encontraron el target por fijación')
    xlabel('Cantidad de fijaciones')
    grid on
subplot(2,1,2)
    siz = (count_img .* 100) ./ 133;
    bar(siz)
    set(gca,'XLim',[1 18])
    title('Porcentaje de imágenes en las que se encontró el target por fijación')
    xlabel('Cantidad de fijaciones')
    grid on