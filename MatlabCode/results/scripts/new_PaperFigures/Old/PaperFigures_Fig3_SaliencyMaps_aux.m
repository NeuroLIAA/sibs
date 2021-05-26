clear all
close all
clc

addpath('~/Dropbox/my_functions/')
addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
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

%% Calc
Fix = []; for j = 1:cfg.NN; Fix(j).fixrank = j; end
for j = 1:cfg.NN
    fprintf('*********************************\n',j)
    filename = sprintf('../matrix/fixated_area_fix_%d_to_%d.mat', j, j);
    if ~exist(filename,'file')
        fixatedPoints(cfg.radii, i, i) % RUN THIS CODE TO CREATE FILE BELOW
    end
%     filename = sprintf('../matrix/fixated_area_fix_1_to_%d.mat', j);
%     if ~exist(filename,'file')
%         fixatedPoints(radii, 1, j) % RUN THIS CODE TO CREATE FILE BELOW
%     end
    tmp = load(filename);
    
    fprintf('Cargue datos de la fijacion %d, \n',j)
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
save('roccurves_i_to_i_20191017.mat','Fix','auccurve','cfg')

