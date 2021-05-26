clear all
close all
clc

%addpath('~/Dropbox/my_functions/')
addpath(genpath('utils'))
addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
load('../matrix/info_per_subj_additional_data.mat'); % loads images, img_order, templates, template_order
addpath(genpath('../static_models/utils/'));

% load('../matrix/initial_fixations.mat');

load('../new_data/new_matrix/info_all_subj.mat');

% data subj
[subjects, ~, subj_order] = unique({info_per_subj_final(:).subj});
subj_order  = subj_order'; 
Nsubj       = length(unique(subj_order));
% data img
[images, ~, img_order] = unique({info_per_subj_final(:).image_name});
img_order  = img_order'; 
Nimg       = length(unique(img_order));
Ntr         = length(info_per_subj_final);
%%
clear info_per_subj templates template_order img_order

%% Creacion de fixatedArea - probamos con area 1
% Si quiero que el radio sea menor (o distinto) a 3
minFix = 2;
maxFix = 3;
radii = 0;
%new_fixatedPoints(radii, minFix, maxFix) % RUN THIS CODE TO CREATE FILE BELOW

%%
load(sprintf('../new_data/new_matrix/fixated_area_fix_%d_to_%d_radii_%d.mat', 2, 3, 0))
clear fixatedAreaTrain fixatedAreaTest

%%
% el estandar... con radios menores que parece ser lo que sugiere el MIT no
% parece funcionar bien
%load(sprintf('../new_data/new_matrix/fixated_area_fix_%d_to_%d.mat', 3, 3))
%clear fixatedAreaTrain fixatedAreaTest

%% load saliency maps
clc
addpath('utils/saliency_mit/code_forMetrics')
modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
modelos_name    = {'MLNet', 'SAM-ResNet', 'SAM-VGG', 'DeepGaze', 'icf', 'Humans', 'Center'};
modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 
modelo_ind = 4;
features = nan(768,1024,length(images));
for im = 1:length(images)
    if strcmp(modelos_list{modelo_ind}, 'center')
        features(:,:,im) = imread('../saliency/center.jpg');
    elseif strcmp(modelos_list{modelo_ind}, 'humans')
        features(:,:,im) = imread(['../saliency/humans_fix_3_to_3/' images{im}]);
    else
        features(:,:,im) = imread(['../saliency/' modelos_list{modelo_ind} '/' images{im}]);
    end
end

%% Pipeline para ver la data de AUC
rng(1);
Nsamp = 5;
sample = unidrnd(Nimg,[1,Nsamp]);
sampleNeg = unidrnd(Nimg,[1,20]);
%%
%visualize_AUC(features(:,:,1),fixatedArea(:,:,1))
%% shuffled
img = 11;
otherMap = sum(fixatedArea(:,:,sampleNeg),3);
[score,tp,fp] = AUC_shuffled(features(:,:,img), fixatedArea(:,:,img), otherMap, 100, 0.05, 1);
visualize_NSS(features(:,:,img),fixatedArea(:,:,img))
%%

for i = 1:Nsamp
    
end



%% No correr
cfg_metrics = [];
cfg_metrics.borji.Nsplits = 500;
cfg_metrics.borji.stepSize = 0.05;
cfg_metrics.shuffled.M = 10;
cfg_metrics.shuffled.Nsplits = 300;
cfg_metrics.shuffled.stepSize = 0.05;

score_AUC_borji     = nan(1,length(images));
score_AUC_judd      = nan(1,length(images));
score_AUC_shuffled  = nan(1,length(images)); 
% tp_borji     = nan(2 + round(1/cfg_metrics.borji.stepSize),length(images));
% fp_borji     = nan(2 + round(1/cfg_metrics.borji.stepSize),length(images));
% tp_judd      = [];
% fp_judd      = [];
% tp_shuffled  = nan(2 + round(1/cfg_metrics.borji.stepSize),length(images));
% fp_shuffled  = nan(2 + round(1/cfg_metrics.borji.stepSize),length(images));

borji = [];
judd = [];
shuff = [];
fprintf('Running saliency.mit.edu measures:\n')
for im = 1:length(images)
    % AUC Borji
    [score_AUC_borji(im),borji(im).tp,borji(im).fp] = ...
        AUC_Borji(double(squeeze(features(:,:,im))/255), squeeze(fixatedArea(:,:,im)), ...
        cfg_metrics.borji.Nsplits, cfg_metrics.borji.stepSize, 0);

    % AUC Judd
    [score_AUC_judd(im),judd(im).tp,judd(im).fp] = AUC_Judd(double(squeeze(features(:,:,im))/255), squeeze(fixatedArea(:,:,im)), 1);
   
    % AUC Shuffled
    % M random integers (without reposition) without the analyzed image
    tic
    indrand = 1:length(images); indrand(im)=[]; indrand = indrand(randperm(length(images)-1, cfg_metrics.shuffled.M));
    otherMap = double(sum(features(:,:,indrand),3)/255);
    [score_AUC_shuffled(im),shuff(im).tp,shuff(im).tp] = ...
        AUC_shuffled(double(squeeze(features(:,:,im))/255), squeeze(fixatedArea(:,:,im)),...
        otherMap, cfg_metrics.shuffled.Nsplits, cfg_metrics.shuffled.stepSize, 0);
    toc 

    if mod(im,5); fprintf('%d, '); end
end
fprintf('Done!\n')
save('saliency_mit_AUC_center.mat')

%% AUC para el modelo de DeepGaze - no correr
load('saliency_mit_AUC.mat')
figure(100); clf
    subplot(1,3,1);
        hold on
            title(sprintf('score AUC borji = %0.2f',nanmean(score_AUC_borji)))
            hist(score_AUC_borji,0:0.05:1)
        hold off
    subplot(1,3,2);
        hold on
            title(sprintf('score AUC judd = %0.2f',nanmean(score_AUC_judd)))
            hist(score_AUC_judd,0:0.05:1)
        hold off        
    subplot(1,3,3);
        hold on
            title(sprintf('score AUC shuffled = %0.2f',nanmean(score_AUC_shuffled)))
            hist(score_AUC_shuffled,0:0.05:1)
        hold off
        
%% Correr desde aca - Gaston - INICIALIZACION DE VARIABLES

fprintf('Loading variables...\n')
todo = load(sprintf('../matrix/fixated_area_fix_%d_to_%d.mat', 2, 12)); %fixatedArea
figure; imagesc(todo.fixatedArea(:,:,1))
cfg_metrics = [];
cfg_metrics.borji.Nsplits = 500;
cfg_metrics.borji.stepSize = 0.05;
cfg_metrics.shuffled.M = 10;
cfg_metrics.shuffled.Nsplits = 500;
cfg_metrics.shuffled.stepSize = 0.05;
log_fix_img = [];
sAUC_allfix = [];

%% Cálculo de sAUC
fprintf('Running saliency.mit.edu measures:\n')
for fix = 2:12 
    load(sprintf('../matrix/fixated_area_fix_%d_to_%d.mat', fix, fix)); %fixatedArea
    
    shuff     = [];
    shuffTodo = [];
    fix
    tic
    for im = 1:length(images)
        % AUC Shuffled
        %tic
        % M random integers (without reposition) without the analyzed image
        indrand = 1:length(images); indrand(im)=[]; indrand = indrand(randperm(length(images)-1, cfg_metrics.shuffled.M));
        otherMap    = sum(fixatedArea(:,:,indrand),3) > 0;
        otherMapTodo= sum(todo.fixatedArea(:,:,indrand),3) > 0;
        
        % tiro error en fijacion 6 e imagen 5
        % no fixationMap
        % Output argument "tp" (and maybe others) not assigned during call to "AUC_shuffled".
        % Error in PaperFigures_Fig3_SaliencyMaps_MITSaliencyMeasures (line 132)
        %    [shuff(im).score, shuff(im).tp, shuff(im).fp] = ...
        
        %[shuff(im).score, shuff(im).tp, shuff(im).fp] = ...
        %    AUC_shuffled(   double(squeeze(features(:,:,im))/255), ...
        %                    squeeze(fixatedArea(:,:,im)),...
        %                    otherMap, cfg_metrics.shuffled.Nsplits, cfg_metrics.shuffled.stepSize, 0);
        
        % chequeamos que la imagen tenga fijaciones
        if ~any(fixatedArea(:,:,im))
            im
        else 
            [shuffTodo(im).score, shuffTodo(im).tp, shuffTodo(im).fp] = ...
                    AUC_shuffled(   double(squeeze(features(:,:,im))/255), ...
                            squeeze(fixatedArea(:,:,im)),...
                            otherMapTodo, cfg_metrics.shuffled.Nsplits, cfg_metrics.shuffled.stepSize, 0);
            %toc
        end
        %else
        %    log_fix_img = [log_fix_img [fix;im]];
        %end

        if mod(im,5); fprintf('%d, '); end
    end
    toc
    
    %a.shuff = shuff;
    a.shuffTodo = shuffTodo;    
    
    sAUC_allfix = [sAUC_allfix a];
end
%save('saliency_mit_AUC_shuffled_deepgaze_todo.mat')
%save('saliency_sauc_center.mat')
% para ver cuanto dio arrayfun(@(x) mean([x.shuffTodo.score]), sAUC_allfix)

%% All fixations compared with sAUC (Gaston) - copiado del posta
 
FixAll = load(sprintf('../matrix/fixated_area_fix_%d_to_%d.mat', 3, 12));
[FixAll.mlnet.meanSalientPercent, FixAll.mlnet.meanTPR, FixAll.mlnet.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'mlnet'); 
[FixAll.sam_resnet.meanSalientPercent, FixAll.sam_resnet.meanTPR, FixAll.sam_resnet.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'sam_resnet'); 
[FixAll.sam_vgg.meanSalientPercent, FixAll.sam_vgg.meanTPR, FixAll.sam_vgg.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'sam_vgg'); 
[FixAll.deepgaze.meanSalientPercent, FixAll.deepgaze.meanTPR, FixAll.deepgaze.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'deepgaze'); 
[FixAll.icf.meanSalientPercent, FixAll.icf.meanTPR, FixAll.icf.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'icf'); 
[FixAll.humans.meanSalientPercent, FixAll.humans.meanTPR, FixAll.humans.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'Humans3'); 
[FixAll.center.meanSalientPercent, FixAll.center.meanTPR, FixAll.center.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, FixAll.fixatedArea, 'Center'); 
    
modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
for i=1:length(modelos_list)
    FixAll.(modelos_list{i}).auc = fun_auc_with_errors(FixAll.(modelos_list{i}).meanSalientPercent/100, FixAll.(modelos_list{i}).meanTPR, FixAll.(modelos_list{i}).errTPR);
end
%% Plot
hold on
H = [];
for i=1:length(modelos_list)
    H(i) = errorbar(FixAll.(modelos_list{i}).meanSalientPercent, FixAll.(modelos_list{i}).meanTPR, FixAll.(modelos_list{i}).errTPR, ...
                        '-', 'Color', modelos_col{i});
end
plot([0 100], [0 1], 'k-')
legend(H,modelos_name,'Location','SouthEast')
xlabel('Percent salient')
ylabel('True positive ratio (TPR)')
set(gca,'XTickLabel',[])
axis([0 100 0 1])
