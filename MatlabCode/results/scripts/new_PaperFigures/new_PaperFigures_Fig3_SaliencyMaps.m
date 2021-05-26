clear all
close all
clc
%%
%addpath('../data_analysis/utils/')
%addpath('../data_analysis/utils/heatmap_code/')

load('../matrix/info_per_subj_additional_data.mat'); % loads images, img_order, templates, template_order
%addpath(genpath('../static_models/utils/'));
addpath('./utils/')

%load('../matrix/initial_fixations.mat');

%load('../matrix/info_per_subj_final.mat');
load('../new_data/new_matrix/info_all_subj.mat') % loads all subjects into info_per_subj_final
[subjects, ~, subj_order] = unique({info_per_subj_final(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj_final);
clear info_per_subj_final
% path donde estan 
src_path = '../new_data/new_matrix/';
out_path = '../new_data/new_matrix/';

salientPercentThr = 1:10:100;
minFix = 3;
maxFix = 3;
radii = 3;

modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
modelos_name    = {'MLNet', 'SAM-ResNet', 'SAM-VGG', 'DeepGaze2', 'ICF', 'Humans', 'Center'};
modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 

%% Fig 3A: Corre bien con el MLNet. Revisarlo para el resto.
if 1
    new_fixatedPoints(radii, minFix, maxFix, src_path, out_path) % RUN THIS CODE TO CREATE FILE BELOW
    load(sprintf('../new_data/new_matrix/fixated_area_fix_%d_to_%d.mat',...
        minFix, maxFix));

    modelo = 'deepgaze'; % 'sam_resnet'; 'sam_vgg'; 
    %trainAndTestModel('all_features', 'allFixationsMap/', 'models/model_all_subjects_fix_3.mat');
    load('../static_models/models/model_all_subjects_fix_3.mat')
    %3
    image_to_print = [1];
    percent_salient_to_print = [2, 4];

    figure(1); clf; set(gcf,'Color','w'); set(gcf,'Position', [60 720 1190 260])
        for im = image_to_print
            subplot(1,3,1);
                imshow(imread(['../images/' images{im}]));
            if strcmp(modelo,'center')
                features = imgaussfilt(imread(['../saliency/center.jpg']), 20);
            else
                features = imgaussfilt(imread(['../saliency/' modelo '/' images{im}]), 20);
            end
            tRange = 0:1:255;

            salientPercentAux = [];        
            for t = tRange
                salientPercentAux(end+1) = sum(sum(features >= t)) * 100 ./ (768*1024);
            end

            idx = 1;
            for p = percent_salient_to_print
                [c, t] = min(abs(salientPercentThr(p)-(salientPercentAux)));
                thr = features >= tRange(t);
                true_positives  = sum(sum(fixatedArea(:,:,im) == 1 & thr == 1));
                false_negatives = sum(sum(fixatedArea(:,:,im) == 1 & thr == 0));

                idx = idx + 1;
                subplot(1, 3, idx);
                    imshow(thr)
                    hold on
                    [xaux, yaux] = find(fixatedArea(:,:,im) == 1 & thr == 1);
                    scatter(yaux, xaux, 'go', 'filled')
                    [xaux, yaux] = find(fixatedArea(:,:,im) == 1 & thr == 0);
                    scatter(yaux, xaux, 'ro', 'filled')
                    title([sprintf('%.2f %% salient, TPR = %.2f', salientPercentAux(t), ...
                        true_positives * 1.0 / (true_positives + false_negatives))])
            end
        end
end

%% Fig 3B. ROC
%% Viene de staticModels.m
% Uso tpr_vs_area.m para calcular lo que necesito para la ROC
% [meanSalientPercent, meanTPR, errTPR, meanFPR, errFPR] = tpr_vs_area(maxFix, images, ...
%                 templates, salientPercentThr, fixatedArea, mode, horizonMode, modelIdx, ...
%                 ModelsJudd, ROCPerformancesJudd, dimsJudd)

% fixatedArea(:,:,test_set): Antes le pasaba un test_set porque entrenaba Judd, 
% pero ahora son todos preentrenados, asi que puedo usar todas las imagenes 

modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};

%% Third fixation (load models)
Fix3 = load(strcat(src_path, sprintf('fixated_area_fix_%d_to_%d.mat', 3, 3))); % Levanta datos
Fix3 = rmfield(Fix3,{'fixatedAreaTrain','fixatedAreaTest','testSubjects'});
[Fix3.mlnet.meanSalientPercent, Fix3.mlnet.meanTPR, Fix3.mlnet.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, Fix3.fixatedArea, 'mlnet'); 
[Fix3.sam_resnet.meanSalientPercent, Fix3.sam_resnet.meanTPR, Fix3.sam_resnet.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, Fix3.fixatedArea, 'sam_resnet'); 
[Fix3.sam_vgg.meanSalientPercent, Fix3.sam_vgg.meanTPR, Fix3.sam_vgg.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, Fix3.fixatedArea, 'sam_vgg'); 
[Fix3.deepgaze.meanSalientPercent, Fix3.deepgaze.meanTPR, Fix3.deepgaze.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, Fix3.fixatedArea, 'deepgaze'); 
[Fix3.icf.meanSalientPercent, Fix3.icf.meanTPR, Fix3.icf.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, Fix3.fixatedArea, 'icf'); 
[Fix3.humans.meanSalientPercent, Fix3.humans.meanTPR, Fix3.humans.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, Fix3.fixatedArea, 'Humans3'); 
[Fix3.center.meanSalientPercent, Fix3.center.meanTPR, Fix3.center.errTPR, ~, ~] = ...
        tpr_vs_area([], images, templates, salientPercentThr, Fix3.fixatedArea, 'Center'); 
    

%% All fixations (load models)
% fixatedPoints(radii, 1, 10) % RUN THIS CODE TO CREATE FILE BELOW
FixAll = load(strcat(src_path,sprintf('fixated_area_fix_%d_to_%d.mat', 3, 12)));
FixAll = rmfield(FixAll,{'fixatedAreaTrain','fixatedAreaTest','testSubjects'});
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
    
%%    
for i=1:length(modelos_list)
    Fix3.(modelos_list{i}).auc = fun_auc_with_errors(Fix3.(modelos_list{i}).meanSalientPercent/100, Fix3.(modelos_list{i}).meanTPR, Fix3.(modelos_list{i}).errTPR);
    FixAll.(modelos_list{i}).auc = fun_auc_with_errors(FixAll.(modelos_list{i}).meanSalientPercent/100, FixAll.(modelos_list{i}).meanTPR, FixAll.(modelos_list{i}).errTPR);
end

%% Third fixation (plot results)
% modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
% modelos_name    = {'MLNet', 'SAM-ResNet', 'SAM-VGG', 'DeepGaze', 'ICF', 'Humans', 'Center'};
% modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 
axpos           = [ 0.10 0.55 0.4 0.4; ...
                    0.55 0.55 0.4 0.4; ...
                    0.10 0.10 0.4 0.4; ...
                    0.55 0.10 0.4 0.4];

figure(3);clf
    set(gcf,'Color','w')
    set(gcf,'Position',[680 360 520 600])
    axes('Position',axpos(1,:))
        hold on
            plot([0 100], [0 1], 'k-')
            H = [];
            for i=1:length(modelos_list)
                [errorPatch, H(i)] = niceBars(Fix3.(modelos_list{i}).meanSalientPercent',...
                                        Fix3.(modelos_list{i}).meanTPR', ...
                                        Fix3.(modelos_list{i}).meanTPR' - Fix3.(modelos_list{i}).errTPR',...
                                        Fix3.(modelos_list{i}).meanTPR' + Fix3.(modelos_list{i}).errTPR',...
                                        modelos_col{i},0.4);
            end
        hold off
        axis square
        box on
        xlabel('Percent salient')
        ylabel('True positive ratio (TPR)')
        set(gca,'YLim',[0 1])
        set(gca,'XLim',[0 100])
        set(gca,'FontName','Arial','FontSize',10)
        
    axes('Position',axpos(2,:))
        hold on
            plot([0 100], [0 1], 'k-')
            H = [];
            for i=1:length(modelos_list)
                [errorPatch, H(i)] = niceBars(FixAll.(modelos_list{i}).meanSalientPercent',...
                                        FixAll.(modelos_list{i}).meanTPR', ...
                                        FixAll.(modelos_list{i}).meanTPR' - FixAll.(modelos_list{i}).errTPR',...
                                        FixAll.(modelos_list{i}).meanTPR' + FixAll.(modelos_list{i}).errTPR',...
                                        modelos_col{i},0.4);
            end
        hold off
        axis square
        box on
        xlabel('Percent salient')
        set(gca,'YLim',[0 1],'YTickLabel',[])
        set(gca,'XLim',[0 100])
        set(gca,'FontName','Arial','FontSize',10)
        
    axes('Position',axpos(3,:))
        hold on
            for i=1:length(modelos_list)
                H=bar(i,Fix3.(modelos_list{i}).auc.mean); set(H,'FaceColor',modelos_col{i})
                errorbar(i,Fix3.(modelos_list{i}).auc.mean, ...
                    Fix3.(modelos_list{i}).auc.mean-Fix3.(modelos_list{i}).auc.err(1),...
                    Fix3.(modelos_list{i}).auc.err(2)-Fix3.(modelos_list{i}).auc.mean,...
                    'Color', modelos_col{i})
            end
        hold off
        ylabel('AUC')
        axis square
        box on
        set(gca,'YLim',[0.5 1.0])
        set(gca,'XLim',[0.5 7.5],'XTickLabel',[])
        set(gca,'XTick',1:7,'XTickLabel',modelos_name,'XTickLabelRotation',45)
        set(gca,'FontName','Arial','FontSize',10)
        
    axes('Position',axpos(4,:))
        hold on
            for i=1:length(modelos_list)
                H=bar(i,FixAll.(modelos_list{i}).auc.mean); set(H,'FaceColor',modelos_col{i})
                errorbar(i,FixAll.(modelos_list{i}).auc.mean, ...
                    FixAll.(modelos_list{i}).auc.mean-FixAll.(modelos_list{i}).auc.err(1),...
                    FixAll.(modelos_list{i}).auc.err(2)-FixAll.(modelos_list{i}).auc.mean,...
                    'Color', modelos_col{i})
            end
        hold off
        axis square
        box on
        set(gca,'YLim',[0.5 1.0],'YTickLabel',[])
        set(gca,'XLim',[0.5 7.5],'XTickLabel',[])
        set(gca,'XTick',1:7,'XTickLabel',modelos_name,'XTickLabelRotation',45)
        set(gca,'FontName','Arial','FontSize',10)
        
%% Fig 3C. (Calculations) AUC vs Fixation Rank
if 0
    % run PaperFigures_Fig3_SaliencyMaps_aux
    modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'humans', 'center'};
%     modelos_list    = {'deepgaze', 'icf', 'humans', 'center', 'noise'};
    NN = 12;
    Fix = []; for j = 1:NN; Fix(j).fixrank = j; end
    parfor j = 1:NN
        fprintf('*********************************\n',j)
        filename = sprintf('../new_data/new_matrix/fixated_area_fix_%d_to_%d.mat', j, j);
        %keyboard
        if ~exist(filename,'file')
            new_fixatedPoints(radii, j, j,src_path,out_path); % RUN THIS CODE TO CREATE FILE BELOW
        end
    %     filename = sprintf('../matrix/fixated_area_fix_1_to_%d.mat', j);
    %     if ~exist(filename,'file')
    %         fixatedPoints(radii, 1, j) % RUN THIS CODE TO CREATE FILE BELOW
    %     end
    %    tmp = load(strcat(src_path,filename));
        tmp = load(filename);
        fprintf('Cargue datos de la fijacion %d, \n',j)
        fprintf('\ty empiezo a calcular la ROC para todos los modelos.\n')

        [Fix(j).mlnet.meanSalientPercent,       Fix(j).mlnet.meanTPR,       Fix(j).mlnet.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'mlnet'); 
        [Fix(j).sam_resnet.meanSalientPercent,  Fix(j).sam_resnet.meanTPR,  Fix(j).sam_resnet.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'sam_resnet'); 
        [Fix(j).sam_vgg.meanSalientPercent,     Fix(j).sam_vgg.meanTPR,     Fix(j).sam_vgg.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'sam_vgg'); 
        [Fix(j).deepgaze.meanSalientPercent,    Fix(j).deepgaze.meanTPR,    Fix(j).deepgaze.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'deepgaze'); 
        [Fix(j).icf.meanSalientPercent,    Fix(j).icf.meanTPR,    Fix(j).icf.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'icf'); 
        [Fix(j).humans.meanSalientPercent,      Fix(j).humans.meanTPR,      Fix(j).humans.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'Humans3'); 
        [Fix(j).center.meanSalientPercent,      Fix(j).center.meanTPR,      Fix(j).center.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'Center'); 
        [Fix(j).noise.meanSalientPercent,      Fix(j).noise.meanTPR,      Fix(j).noise.errTPR, ~, ~] = ...
                tpr_vs_area([], images, [], salientPercentThr, tmp.fixatedArea, 'Noise'); 

        fprintf('\tahora calculo la AUC para todos los modelos.\n')
        for i=1:length(modelos_list)
            Fix(j).(modelos_list{i}).auc = fun_auc_with_errors(Fix(j).(modelos_list{i}).meanSalientPercent/100, Fix(j).(modelos_list{i}).meanTPR, Fix(j).(modelos_list{i}).errTPR);
        end
    end

    auccurve = [];
    for j = 1:NN
        for i=1:length(modelos_list)
            auccurve.(modelos_list{i}).mean(j)  = Fix(j).(modelos_list{i}).auc.mean;
            auccurve.(modelos_list{i}).lower(j) = Fix(j).(modelos_list{i}).auc.err(1);
            auccurve.(modelos_list{i}).upper(j) = Fix(j).(modelos_list{i}).auc.err(2);
        end
    end

    % save('roccurves/roccurves_i_to_i.mat','Fix','auccurve')
    save('roccurves/roccurves_2_to_12_20200512.mat','Fix','auccurve')
else
    % load('roccurves/roccurves_i_to_i.mat'); NN = 10;
    load('roccurves/roccurves_2_to_12_20200512.mat'); NN = 12;
end

%% Fig 3C. AUC vs Fixation Rank
load('roccurves_2_to_12_20200512.mat')
% modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
% modelos_name    = {'MLNet (Cornia et al., 2016)', ...
%     'SAM-ResNet (Cornia et al., 2018)', 'SAM-VGG (Cornia et al., 2018)',...
%     'DeepGaze2 (K�mmerer et al., 2017)', 'ICF (K�mmerer et al., 2017)', ...
%     'Humans (3rd Fixation)', 'Center'};
modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 
% NN = cfg.NN;
figure(4); clf
    set(gcf,'Color','w')
    set(gcf,'Position',[550 550 690 420])
        hold on
            H = [];
            plot([0 NN+1], [0.5 0.5], 'k--')
            plot([3 3], [0.4 0.9], 'k--')
            for i=1:length(modelos_list)
                [errorPatch, H(i)] = niceBars(2:NN, auccurve.(modelos_list{i}).mean(2:end), ...
                                        auccurve.(modelos_list{i}).lower(2:end),...
                                        auccurve.(modelos_list{i}).upper(2:end),...
                                            modelos_col{i},0.6);
            end
        hold off
        set(gca,'XLim',[2 NN])
        set(gca,'YLim',[0.4 1])
        xline(3,'--')
        legend(H,modelos_name,'Location','NorthEastOutside')
        legend boxoff
        xlabel('Fixation Rank')
        ylabel('AUC')
        box on

%% Fig 3B Third fixation (plot results)
% figure(3);clf
%     set(gcf,'Color','w')
% %     subplot(2,2,1)
%         hold on
%             plot([0 100], [0 1], 'k-')
%             H = [];
%             for i=1:length(modelos_list)
%                 [errorPatch,H(i)] = niceBars(Fix(3).(modelos_list{i}).meanSalientPercent',...
%                                         Fix(3).(modelos_list{i}).meanTPR', ...
%                                         Fix(3).(modelos_list{i}).errTPR', ...
%                                         modelos_col{i},0.7);
% %                 H(i) = errorbar(Fix(3).(modelos_list{i}).meanSalientPercent, Fix(3).(modelos_list{i}).meanTPR, Fix(3).(modelos_list{i}).errTPR, ...
% %                                     '-', 'Color', modelos_col{i});
%             end
%         hold off
% %         legend(H,modelos)
%         xlabel('Percent salient')
%         ylabel('True positive ratio (TPR)')
%         set(gca,'XTickLabel',[])
%         axis([0 100 0 1])
% %         box on
        
%     subplot(2,2,2)
%         hold on
%             for i=1:length(modelos_list)
%                 H=bar(i,Fix(3).(modelos_list{i}).auc.mean); set(H,'FaceColor',modelos_col{i})
%                 errorbar(i,Fix(3).(modelos_list{i}).auc.mean, ...
%                     Fix(3).(modelos_list{i}).auc.mean-Fix(3).(modelos_list{i}).auc.err(1),...
%                     Fix(3).(modelos_list{i}).auc.err(2)-Fix(3).(modelos_list{i}).auc.mean,...
%                     'Color', modelos_col{i})
%             end
%         hold off
%         ylabel('AUC')
%         set(gca,'YLim',[0.5 1.0],'XTickLabel',[])
%         set(gca,'XTick',1:6,'XLim',[0.5 7.5],'XTickLabel',[])
%         box on
% %         set(gca,'XTick',1:5,'XTickLabel',modelos,'XTickLabelRotation',45)


%% Fig 3D. AUC vs Fixation Rank (accum)
% % load('roccurves_i_to_10.mat')
% % modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'humans', 'center'};
% % modelos_name    = {'MLNet', 'SAM-ResNet', 'SAM-VGG', 'DeepGaze', 'Humans', 'Center'};
% % modelos_col     = {'r', 'b', 'c', 'g', 'k', [.7 .7 .7]}; 
% 
% load('roccurves_i_to_10_20191002.mat')
% modelos_list    = {'mlnet', 'humans', 'center', 'noise'};
% modelos_name    = {'MLNet', 'Humans', 'Center', 'Noise'};
% modelos_col     = {'r', 'k', [.7 .7 .7], 'm'}; 
% 
% figure(4); clf
%     set(gcf,'Color','w')
%         hold on
%             H = [];
%             for i=1:length(modelos_list)
%                 H(i) = errorbar(1:NN, auccurve.(modelos_list{i}).mean, ...
%                                         auccurve.(modelos_list{i}).lower-auccurve.(modelos_list{i}).mean,...
%                                         auccurve.(modelos_list{i}).mean-auccurve.(modelos_list{i}).upper,...
%                                         '-', 'Color', modelos_col{i});
%             end
% %             plot([0 NN+1], [0.5 0.5], 'k-')
%             set(gca,'XLim',[0.5 NN+0.5])
%             set(gca,'YLim',[0.4 1])
% %             legend(H,modelos_name,'Location','UpperEast')
%             xlabel('Fixation Rank')
%             ylabel('AUC')
%             box on
%         hold off
