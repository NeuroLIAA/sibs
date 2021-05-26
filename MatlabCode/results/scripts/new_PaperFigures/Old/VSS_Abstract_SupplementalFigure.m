%% Figura compuesta de la 1,3 y 4 para la informacion suplementaria de VSS 2020
% sacado de fig3_saliencyMaps_aux

clear all
close all
clc

addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
addpath('../PaperFigures/utils/')
addpath(genpath('../static_models/utils/'));
load('../matrix/info_per_subj_additional_data.mat'); % loads images, img_order, templates, template_order
addpath(genpath('../static_models/utils/'));

load('../matrix/initial_fixations.mat');

load('../matrix/info_per_subj_final.mat');
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);

cfg=[];
cfg.modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center', 'noise'};
cfg.modelos_name    = {'MLNet', 'SAM-ResNet', 'SAM-VGG', 'DeepGaze', 'icf', 'Humans', 'Center', 'Noise'};
cfg.modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 

cfg.salientPercentThr = 1:7:100;
cfg.minFix  = 3;
cfg.maxFix  = 3;
cfg.radii   = 3;

cfg.NN      = 12;

%% Figura 0 - Mapeo a la grilla (de Figura 1 Paper)

delta       = 32;
image_size  = [768 1024];
grid_size   = image_size/delta;

a = dir('../matrix/subjects/*.mat'); filenames = {a.name};
info_per_subj = [];
for i=1:length(filenames)
%     tmp = load(['../matrix/subjects/' filenames{i}]);
%     if ~strcmp(tmp.info_per_subj(1).subj,'MVA')
%         info_per_subj = [info_per_subj tmp.info_per_subj];
%     end
    tmp = load(['../matrix/subjects/' filenames{i}]);
    tmp.info_per_subj = subjMapFixationToMatrix( tmp.info_per_subj, path, delta, image_size );
    info_per_subj = [info_per_subj tmp.info_per_subj];
end
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);

%% COMPLETAR

%% Figura A - AUC vs Fixation Rank - Calculo datos generales (de Figura 3C Paper)

% FixAll = load(sprintf('../matrix/fixated_area_fix_%d_to_%d.mat', 3, 12));
% % [FixAll.mlnet.meanSalientPercent, FixAll.mlnet.meanTPR, FixAll.mlnet.errTPR, ~, ~] = ...
% %         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'mlnet'); 
% % [FixAll.sam_resnet.meanSalientPercent, FixAll.sam_resnet.meanTPR, FixAll.sam_resnet.errTPR, ~, ~] = ...
% %         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'sam_resnet'); 
% [FixAll.sam_vgg.meanSalientPercent, FixAll.sam_vgg.meanTPR, FixAll.sam_vgg.errTPR, ~, ~] = ...
%         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'sam_vgg'); 
% [FixAll.deepgaze.meanSalientPercent, FixAll.deepgaze.meanTPR, FixAll.deepgaze.errTPR, ~, ~] = ...
%         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'deepgaze'); 
% [FixAll.icf.meanSalientPercent, FixAll.icf.meanTPR, FixAll.icf.errTPR, ~, ~] = ...
%         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'icf'); 
% [FixAll.humans.meanSalientPercent, FixAll.humans.meanTPR, FixAll.humans.errTPR, ~, ~] = ...
%         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'Humans3'); 
% [FixAll.center.meanSalientPercent, FixAll.center.meanTPR, FixAll.center.errTPR, ~, ~] = ...
%         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'Center'); 
% [FixAll.noise.meanSalientPercent, FixAll.humans.meanTPR, FixAll.humans.errTPR, ~, ~] = ...
%         tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'noise'); 
%     
% modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center', 'noise'};
% for i=1:length(modelos_list)
%     FixAll.(modelos_list{i}).auc = fun_auc_with_errors(FixAll.(modelos_list{i}).meanSalientPercent/100, FixAll.(modelos_list{i}).meanTPR, FixAll.(modelos_list{i}).errTPR);
% end

%% Fig A. (Calculations) AUC vs Fixation Rank Data (sacado de Fig3_SaliencyMaps_aux.m)

Fix = []; for j = 1:cfg.NN; Fix(j).fixrank = j; end
for j = 1:cfg.NN
    fprintf('*********************************\n',j)
    filename = sprintf('../matrix/fixated_area_fix_%d_to_%d.mat', j, j);
    if ~exist(filename,'file')
        fprintf('Creo archivo de datos de la fijacion %d, \n',j)
        fixatedPoints(cfg.radii, i, i) % RUN THIS CODE TO CREATE FILE BELOW
    end
%     filename = sprintf('../matrix/fixated_area_fix_1_to_%d.mat', j);
%     if ~exist(filename,'file')
%         fixatedPoints(radii, 1, j) % RUN THIS CODE TO CREATE FILE BELOW
%     end
    tmp = load(filename);
    
    fprintf('Cargo datos de la fijacion %d, \n',j)
    fprintf('\ty empiezo a calcular la ROC para todos los modelos.\n')
    
    [Fix(j).mlnet.meanSalientPercent,       Fix(j).mlnet.meanTPR,       Fix(j).mlnet.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'mlnet'); 
    [Fix(j).sam_resnet.meanSalientPercent,  Fix(j).sam_resnet.meanTPR,  Fix(j).sam_resnet.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'sam_resnet'); 
    [Fix(j).sam_vgg.meanSalientPercent,     Fix(j).sam_vgg.meanTPR,     Fix(j).sam_vgg.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'sam_vgg'); 
    [Fix(j).deepgaze.meanSalientPercent,    Fix(j).deepgaze.meanTPR,    Fix(j).deepgaze.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'deepgaze'); 
    [Fix(j).icf.meanSalientPercent,    Fix(j).icf.meanTPR,    Fix(j).icf.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'icf'); 
    [Fix(j).humans.meanSalientPercent,      Fix(j).humans.meanTPR,      Fix(j).humans.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'Humans3'); 
    [Fix(j).center.meanSalientPercent,      Fix(j).center.meanTPR,      Fix(j).center.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'Center'); 
    [Fix(j).noise.meanSalientPercent,      Fix(j).noise.meanTPR,      Fix(j).noise.errTPR, ~, ~] = ...
            tpr_vs_area([], images, [], cfg.salientPercentThr, tmp.fixatedArea, 'Noise'); 

    fprintf('\tahora calculo la AUC para todos los modelos.\n')
    for i=1:length(cfg.modelos_list)
        Fix(j).(cfg.modelos_list{i}).auc = fun_auc_with_errors(Fix(j).(cfg.modelos_list{i}).meanSalientPercent/100, Fix(j).(cfg.modelos_list{i}).meanTPR, Fix(j).(cfg.modelos_list{i}).errTPR);
    end
end


auccurve = [];
for j = 1:cfg.NN
    for i=1:length(cfg.modelos_list)
        auccurve.(cfg.modelos_list{i}).mean(j)  = Fix(j).(cfg.modelos_list{i}).auc.mean;
        auccurve.(cfg.modelos_list{i}).lower(j) = Fix(j).(cfg.modelos_list{i}).auc.err(1);
        auccurve.(cfg.modelos_list{i}).upper(j) = Fix(j).(cfg.modelos_list{i}).auc.err(2);
    end
end

% save('roccurves_i_to_i.mat','Fix','auccurve')
% save('roccurves_i_to_10_20191002.mat','Fix','auccurve','cfg')
save('roccurves_i_to_i_20191205.mat','Fix','auccurve','cfg')

%% Fig A. AUC vs Fixation Rank - (sacado de Fig3_SaliencyMaps.m)
%load('roccurves_i_to_i_20191205.mat')
addpath('../PaperFigures/utils/')
modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
modelos_name    = {'MLNet (Cornia et al., 2016)', ...
    'SAM-ResNet (Cornia et al., 2018)', 'SAM-VGG (Cornia et al., 2018)',...
    'DeepGaze2 (Kümmerer et al., 2017)', 'ICF (Kümmerer et al., 2017)', ...
    'Humans (3rd Fixation)', 'Center'};
modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 
NN = cfg.NN;
figure(4); clf
    set(gcf,'Color','w')
        hold on
            H = [];
            for i=1:length(modelos_list)
                [errorPatch,H(i)] = niceBars(1:NN, auccurve.(modelos_list{i}).mean, ...
                                        auccurve.(modelos_list{i}).lower,...
                                        auccurve.(modelos_list{i}).upper,...
                                        modelos_col{i},0.7);
                
%                 H(i) = errorbar([1:NN] + 0.01*i, auccurve.(modelos_list{i}).mean, ...
%                                         auccurve.(modelos_list{i}).lower-auccurve.(modelos_list{i}).mean,...
%                                         auccurve.(modelos_list{i}).mean-auccurve.(modelos_list{i}).upper,...
%                                         '-', 'Color', modelos_col{i});
            end
%             plot([0 NN+1], [0.5 0.5], 'k-')
        hold off
        set(gca,'XLim',[0.5 NN+0.5])
        set(gca,'YLim',[0.5 1])
        legend(H,modelos_name,'Location','NorthEast')
        legend boxoff
        xlabel('Fixation Rank')
        ylabel('AUC')
%         box on

%% Fig B. Proportion of targets found according max saccades (Sacado de Fi4_DynamicModels_GeneralComparison_v2)
load('../matrix/info_per_subj_final.mat');
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);

a=dir('../matrix/images/*.mat');    filenames_img = {a.name}'; Nimg = length(filenames_img);
a=dir('../matrix/subjects/*.mat');  filenames_subj = {a.name};

% models = fun_define_models;
models = fun_define_models_tmp(3);
nntrthr     = 20;
%% Fig B - Model Data - Number of fixations predicted by the models

for ind_model = 1:length(models)
    Nfix_img_model = nan(Nimg,1);
    for ind_img=1:Nimg
        path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                                    models(ind_model).prior,...
                                    models(ind_model).searcher,...
                                    models(ind_model).params,...
                                    ind_img);
        fprintf('%s\n',path);
        if ind_img ~= 132
            load(path); % scanpath = ( (x,y) , fixation number )
            Nfix_img_model(ind_img)  = size(scanpath,1) - 1;
        end
    end
    models(ind_model).Nfix_img_model = Nfix_img_model;
end

% Proportion of targets found as function of the number of saccades allowed for the models
exps_thr = nan(Nimg,Nsubj);
for ind_img=1:Nimg % images
    path = char(strcat('../matrix/images/info_per_subj_img_', num2str(ind_img), '.mat')); fprintf('%s\n',path);
    load(path);
    for  ind_subj = 1:length(info_per_subj)
        exps_thr(ind_img,ind_subj)   = info_per_subj(ind_subj).exp_data.nsaccades_thr;
    end
end

for ind_model=1:length(models)
    target_found_img_model     = zeros(Nimg,4);
%     target_notfound_img_model  = zeros(Nimg,4);
    for ind_img=1:Nimg
        path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                                    models(ind_model).prior,...
                                    models(ind_model).searcher,...
                                    models(ind_model).params,...
                                    ind_img);
        fprintf('%s\n',path);    
        if ind_img ~= 132
            load(path); % scanpath = ( (x,y) , fixation number )
            
            % No me convence esto... hacerlo bien. Tengo que empezar con las fijaciones de cero!!!            
            count       = length(scanpath);% - 1;     % Numero de sacadas
            tmpcountthr = exps_thr(ind_img,:);      % Numero de sacadas permitidas (igual que el experimento).
            tmpcountthr = tmpcountthr(~isnan(tmpcountthr));% Numero de sacadas permitidas (igual que el experimento).
            countthr    = [2 4 8 12];
            for j=1:length(countthr)
                tmp = (count <= tmpcountthr+ 1);
                target_found_img_model(ind_img,j) = mean( tmp( tmpcountthr==countthr(j) ) );
            end            
        end        
    end
    models(ind_model).target_found_img_model    = target_found_img_model;
    models(ind_model).P_target_found_model      = mean(target_found_img_model([1:131 133:134],:));
    models(ind_model).P_target_found_model_std  = std(target_found_img_model([1:131 133:134],:));
end

%%






























