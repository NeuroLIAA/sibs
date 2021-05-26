clear all
close all
clc

addpath('~/Dropbox/my_functions/')
addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
load('../matrix/info_per_subj_additional_data.mat'); % loads images, img_order, templates, template_order
addpath(genpath('../static_models/utils/'));

% load('../matrix/initial_fixations.mat');

load('../matrix/info_per_subj_final.mat');
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);
clear info_per_subj templates template_order img_order

%%
% Si quiero que el radio sea menor (o distinto) a 3
% minFix = 3;
% maxFix = 3;
% radii = 0;
% fixatedPoints(radii, minFix, maxFix) % RUN THIS CODE TO CREATE FILE BELOW
% load(sprintf('../matrix/fixated_area_fix_%d_to_%d_radii_%d.mat', 3, 3, 0))
% figure; imagesc(fixatedArea(:,:,1))

% el estandar... con radios menores que parece ser lo que sugiere el MIT no
% parece funcionar bien
load(sprintf('../matrix/fixated_area_fix_%d_to_%d.mat', 3, 3))

%% load saliency maps
clc
addpath('utils/saliency_mit/code_forMetrics')
modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
modelos_name    = {'MLNet', 'SAM-ResNet', 'SAM-VGG', 'DeepGaze', 'icf', 'Humans', 'Center'};
modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 

features = nan(768,1024,length(images));
for im = 1:length(images)
    features(:,:,im) = imread(['../saliency/' modelos_list{4} '/' images{im}]);
end

%%
cfg_metrics = [];
cfg_metrics.borji.Nsplits = 500;
cfg_metrics.borji.stepSize = 0.05;
cfg_metrics.shuffled.M = 10;
cfg_metrics.shuffled.Nsplits = 500;
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
save('saliency_mit_AUC.mat')

figure(100); clf
    subplot(3,3,1);
        hold on
            title(sprintf('score_AUC_borji = %0.2f',nanmean(score_AUC_borji)))
            hist(score_AUC_borji,0:0.05:1)
        hold off
    subplot(3,3,2);
        hold on
            title(sprintf('score_AUC_judd = %0.2f',nanmean(score_AUC_judd)))
            hist(score_AUC_judd,0:0.05:1)
        hold off        
    subplot(3,3,3);
        hold on
            title(sprintf('score_AUC_shuffled = %0.2f',nanmean(score_AUC_shuffled)))
            hist(score_AUC_shuffled,0:0.05:1)
        hold off        