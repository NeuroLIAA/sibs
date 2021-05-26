%% Chequeo de orden de fijaciones - 2/4 - Gaston

clear all
close all
clc

%% Chequeo nuevo


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

% Para variar priors case = 6
% Para variar serchers case = 5
% models = fun_define_models_tmp(5);
trials_tmp                      = load(strcat(src_path,'info_all_subj.mat'));
[ids_images, ~, images_order] = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                    = images_order';
image_size                      = trials_tmp.info_per_subj_final(1).image_size;
delta                           = 32;
grid_size                       = image_size/delta;
trials                          = new_subjMapFixationToMatrix(trials_tmp.info_per_subj_final, '', delta, image_size);
clear trials_tmp

%% Chequeo de funcionamiento de new_subjMapFixationToMatrix
% - coincidan la cantidad de fijaciones reducidas vs fix
% - identificar los casos con una sola fijacion y porque
% - en base a lo anterior, estamos teniendo en cuenta las primera fijacion
% forzada?
check = [];
n_
for i=1:length(trials)
    % fix_red vs fix
    if length(trials(i).x)~=length(trials(i).fixatio) 


%% VIEJOS CHEQUEOS
addpath('../data_analysis/utils/')
addpath('../dynamic_models/utils/')
addpath('../data_analysis/utils/heatmap_code/')
addpath('./utils/')
src_path =  '../new_data/new_matrix/';

info_aux = load(strcat(src_path,'info_all_subj.mat'));
[subjects, ~, subj_order] = unique({info_aux.info_per_subj_final(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_aux.info_per_subj_final);

img_dir = dir(strcat(src_path,'sinfo_img/*.mat'));    
filenames_img = {img_dir.name}'; 
Nimg = length(filenames_img);
subj_dir = dir(strcat(src_path,'sinfo_subj/*.mat'));
filenames_subj = {subj_dir.name};
%Nsubj = length(filenames_subj);
% definimos los modelos a comparar - en fig 4 case 3

delta = 32;
image_size = size(imread('../images/grayscale_100_oliva.jpg'));
% Cargamos la reducciones de las imagenes
img_ind = 1;
info_img = load(strcat(src_path,'sinfo_img/info_per_img_',num2str(img_ind),'.mat'));
info_per_subj = new_subjMapFixationToMatrix(info_img.info_per_img, '', delta, image_size);

%% CHEQUEO 1 - Orden de las fijaciones en FIXATIONS
clc
errores_x = [];
errores_y = [];
errores_x_ind = [];
errores_y_ind = [];

% cuales son los indicies de las fijaciones - a testear
ind_x = 1;
ind_y = 2;

% image size de esa imagen es [768 1024] es decir [ylim xlim] 
for ind_img=1:length(info_per_subj)
    
    x = max(info_per_subj(ind_img).fixations(:,ind_x));
    if x > image_size(2)
        %disp(i);
        %disp(x);
        errores_x =[errores_x x];
        errores_x_ind =[errores_x_ind ind_img];
    end
    
    y = max(info_per_subj(ind_img).fixations(:,ind_y));
    if y > image_size(1)
        %disp(i);
        %disp(y);
        errores_y =[errores_y y];
        errores_y_ind =[errores_y_ind ind_img];
    end
    
    % ploteo
    set(gca, 'ydir', 'reverse')
    scatter(x,y)
    hold on
    if ind_img > 1000
        xline(image_size(2),'-.r');
        yline(image_size(1),'-.b');
        break
    end
end

% 2/4 - RESULTADO en info_per_subj las FIXATIONS tienen el orden [X, Y]
%                   x_ind = 1;   y_ind = 2;
% Ojo que los modelos no, por eso estaba dando tan mal - HAY QUE AJUSTARLO

%% CHEQUEO 2 - Orden las fijaciones en matrix_reduced
clc
grid_size   = image_size/delta;
errores_x = [];
errores_y = [];
errores_x_ind = [];
errores_y_ind = [];

% cuales son los indicies de las fijaciones - a testear
ind_x = 1;
ind_y = 2;

% image size de esa imagen es [768 1024] es decir [ylim xlim] 
for ind_img=1:length(info_per_subj)
    %try
        x = max(info_per_subj(ind_img).fixations_matrix_reduced(:,ind_x));
    %catch e
    %    keyboard
    %end
    if x > grid_size(2)
        %disp(i);
        %disp(x);
        errores_x =[errores_x x];
        errores_x_ind =[errores_x_ind ind_img];
    end
    
    y = max(info_per_subj(ind_img).fixations_matrix_reduced(:,ind_y));
    if y > grid_size(1)
        %disp(i);
        %disp(y);
        errores_y =[errores_y y];
        errores_y_ind =[errores_y_ind ind_img];
    end
    
    % ploteo
    set(gca, 'ydir', 'reverse')
    set(gca,'XLim',[-1 33], 'YLim', [-1, 25])
    scatter(x,y)
    hold on
    if ind_img > 1000
        xline(image_size(2),'-.r');
        yline(image_size(1),'-.b');
        break
    end
end

% 2/4 - RESULTADO: Si se mantiene el orden pero tuve un error
% ind = 617 info_per_subj no tiene fijaciones

%% CHEQUEO 3 MODELOS

clc
models = fun_define_models_tmp(5);

%% Chequeo como estan guardados los scanpaths para los modelos

Nimg = 134;
grid_size   = image_size/delta;
errores_x = [];
errores_y = [];
errores_x_ind = [];
errores_y_ind = [];
% cuales son los indicies de las fijaciones - a testear
ind_x = 2;
ind_y = 1;

for ind_model=1:length(models)
    x_all = [];
    y_all = [];
    for ind_img=1:Nimg
        path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                                    models(ind_model).prior,...
                                    models(ind_model).searcher,...
                                    models(ind_model).params,...
                                    ind_img);
        fprintf('%s\n',path);
        load(path)
        x = max(scanpath(:,ind_x));
        y = max(scanpath(:,ind_y));
        if x > grid_size(2)
            %disp(i);
            %disp(x);
            errores_x =[errores_x x];
            errores_x_ind =[errores_x_ind ind_img];
        end
        if y > grid_size(1)
            %disp(i);
            %disp(y);
            errores_y =[errores_y y];
            errores_y_ind =[errores_y_ind ind_img];
        end
        x_all = [x_all x];
        y_all = [y_all y];
    end
    %if i > 2000
    %    xline(image_size(2),'-.r');
    %    yline(image_size(1),'-.b');
    %break
    scatter(x_all,y_all)
    xline(grid_size(2),'-.r');
    yline(grid_size(1),'-.b');
    keyboard
    hold off
end

% 6/4 - Entonces, el problema es que cuando cargamos los modelos las
% variables fixations estan al reves:
% humanos fixations = [X, Y]
% modelos fixations = [Y, X]
% Voy a ver donde lo invierto