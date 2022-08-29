% Figure 4 - MultiMatch graph
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

trials_tmp                    = load(strcat(src_path, 'info_all_subj.mat'));
[ids_images, ~, images_order] = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                  = images_order';
image_size                    = trials_tmp.info_per_subj_final(1).image_size;
delta                         = 32;
grid_size                     = image_size/delta;
trials                        = reduce_scanpaths(trials_tmp.info_per_subj_final, delta, image_size);
clear trials_tmp

mm_names = {'vectorSim','directionSim','lengthSim','positionSim'};
mm_titles = {'MultiMatch Vector', 'MultiMatch Direction', 'MultiMatch Length', 'MultiMatch Position'};

%% Load subjs results and model metrics

addpath('./results_metrics/')
% var_name -> multimatch_bh, mean_dist_img, std_dist_img
load('mm_bh_reduced.mat')

%% Load models

graph_type = 'priors-ibs';
save_plots = 1; % if saving plots will show them and then close them

switch graph_type
    case '5searchers'
        load('mm_hm_reduced_5searchers.mat')
        short_names    = {'cIBS+DGII','sIBS+DGII','IBS+DGII','Greedy+DGII','SaliencyB.'};
        fig_preffix    = 'Fig4_5searcher_boxplot_';
        fig_preffix_1  = 'Fig4_5searcher_mm_all_';
        fig_preffix_2  = 'Fig4_5searcher_scatter_';
        path_to_save   = 'searchers/';
    case 'priors-correlation'
        load('mm_hm_reduced_priors-correlation.mat')
        short_names    = {'cIBS+DGII','cIBS+Center','cIBS+Flat','cIBS+Noisy'};
        fig_preffix    = 'Fig4_priors_boxplot_';
        fig_preffix_1  = 'Fig4_prior_mm_all_';
        fig_preffix_2  = 'Fig4_5searcher_scatter_';
        path_to_save   = 'priors-cIBS/';
    case 'priors-ssim'
        load('mm_hm_reduced_priors-ssim.mat')
        short_names    = {'sIBS+DGII','sIBS+Center','sIBS+Noisy','sIBS+Flat'};
        fig_preffix    = 'Fig4_prior_sibs_boxplot_';
        fig_preffix_1  = 'Fig4_prior_sibs_mm_all_';
        fig_preffix_2  = 'Fig4_prior_sibs_scatter_';
        path_to_save   = 'priors-sIBS/';
    case 'priors-ibs'
        load('mm_hm_reduced_priors-ibs.mat')
        short_names    = {'IBS+DGII','IBS+Center','IBS+Noisy','IBS+Flat'};
        fig_preffix    = 'Fig4_prior_ibs_boxplot_';
        fig_preffix_1  = 'Fig4_prior_ibs_mm_all_';
        fig_preffix_2  = 'Fig4_prior_ibs_scatter_';
        path_to_save   = 'priors-IBS/';
end

%% Calculate subjects statistics

ind_model = 1;
mm_mean   = nanmean(models(ind_model).mean_dist_img);
mm_std    = nanstd(models(ind_model).mean_dist_img);
gold_standard = mean(mean_dist_img);

metric_ind = 1;
x = models(ind_model).mean_dist_img(:,metric_ind);
y = mean_dist_img(:,metric_ind);
not_nan = isfinite(x) & isfinite(y); 

%% BOXPLOT Plot MultiMatch 4 dimensions models comparison

%close all
for mm = 1:length(mm_names)
    figure(mm+400); 
        group = [];
        for mod = 1:length(models)
            group = [group models(mod).(mm_names{mm})];
        end
        boxplot(group, 'Color','k','Notch','on','Labels', short_names) 
        title(mm_titles{mm})
        ylabel(mm_names{mm})
        
        h = findobj(gca,'Tag','Box');
        colors = {models.cols};
        % findobj gets boxplot reversed so I reverse colors with 5-j
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors{length(models)+1-j},'FaceAlpha',.5)
        end
        
        % add humans data
        hold on
        betweenh_mean = nanmean(multimatch_bh.(mm_names{mm}));
        yline(betweenh_mean,'b--')
        hold off
        
        set(gcf,'name', strcat(fig_preffix, mm_names{mm}))
        
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        
        pause(0.2)
end

%% Save figs

if save_plots
    disp('guardando figuras...')
    FolderFigsName = './../figs/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    SAVEPATH = strcat(FolderFigsName, path_to_save, 'supp/');
    if ~exist(SAVEPATH,'dir') 
        mkdir(SAVEPATH);
    end
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, strcat(SAVEPATH, FigName, '.fig'));
      %saveas(FigHandle, strcat(SAVEPATH, FigName, '.svg'));
      print(FigHandle,strcat(SAVEPATH, FigName,'.png'),'-dpng','-r400')
    end
    close all
end

%% SCATTERPLOT MultiMatch 4 dimensions models comparison

for ind_metric=1:4
    figure(410+ind_metric)
    set(gcf,'Renderer', 'painters', 'Position', [100 100 800 450])
    set(gcf,'Color','w')
        Nimg_thr = 15;
            hold on
                filter_img  = adentro_subj > Nimg_thr;
                Nimg_ok     = sum(filter_img);
                % plot subjects values
                plot(zeros(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(mean_dist_img(filter_img,ind_metric),2),'.', 'Color','k','MarkerSize',8)
                for ind_model=1:length(models)
                    filter_img  = models(ind_model).adentro_subj >= Nimg_thr;
                    Nimg_ok     = sum(filter_img);
                    plot(ind_model*ones(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(models(ind_model).mean_dist_img(filter_img,ind_metric),2),...
                        '.','Color', models(ind_model).cols,'MarkerSize',8)
                end
                mean_mm_distance_img = mean_dist_img(:, ind_metric);
                % drop images with no views (4 and 132)
                mean_mm_distance_img_nonzero = mean_mm_distance_img(mean_mm_distance_img~=0);
                % calculate subject stats
                mean_mm_distance = mean(mean_mm_distance_img_nonzero);
                std_mm_distance = std(mean_mm_distance_img_nonzero);
                % debug
                % size(mean_mm_metric_img_nonzero)
                plot([-1 length(models)+1], (mean_mm_distance)*[1 1],'k-')
                plot([-1 length(models)+1], (mean_mm_distance + std_mm_distance)*[1 1],'k--')
                plot([-1 length(models)+1], (mean_mm_distance - std_mm_distance)*[1 1],'k--')

                % plot mean values for models
                for ind_model=1:length(models)
                    plot(ind_model+0.15,nanmean(nanmean(models(ind_model).mean_dist_img(:,ind_metric),2)),'o','Color',models(ind_model).cols)
                end
            hold off
            ylabel(mm_names(ind_metric))
            set(gca,'XLim',[-1 length(models)+1],'XTickLabel',[' ', 'Subjects', short_names, ' '])
          set(gcf,'name', strcat(fig_preffix_2, mm_names{ind_metric}))
        pause(0.1)    
end

%% SCATTERPLOT Mean MultiMatch across 4 dimensions models comparison

figure(420)
    set(gcf,'Renderer', 'painters', 'Position', [100 100 800 450])
    set(gcf,'Color','w')
        % consider only trials with ...(completar)       
        Nimg_thr = 15;
            hold on
                filter_img  = adentro_subj > Nimg_thr;
                Nimg_ok     = sum(filter_img);
                % plot subjects values
                plot(zeros(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(mean_dist_img(filter_img,1:end-1),2),'.', 'Color','k','MarkerSize',8)
                for ind_model=1:length(models)
                    filter_img  = models(ind_model).adentro_subj >= Nimg_thr;
                    Nimg_ok     = sum(filter_img);
                    plot(ind_model*ones(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(models(ind_model).mean_dist_img(filter_img,1:end-1),2),...
                        '.','Color', models(ind_model).cols,'MarkerSize',8)
                end
                mean_mm_distance_img = mean_dist_img(:, 1:end-1);
                % drop images with no views (4 and 132)
                mean_mm_distance_img([4 132],:) = [];
                % mean_mm_distance_img_nonzero = mean_mm_distance_img(mean_mm_distance_img~=0);
                % calculate subject stats
                mean_mm_distance = mean(mean(mean_mm_distance_img,2));
                std_mm_distance = std(mean(mean_mm_distance_img,2));
                % debug
                % size(mean_mm_metric_img_nonzero)
                plot([-1 length(models)+1], (mean_mm_distance)*[1 1],'k-')
                plot([-1 length(models)+1], (mean_mm_distance + std_mm_distance)*[1 1],'k--')
                plot([-1 length(models)+1], (mean_mm_distance - std_mm_distance)*[1 1],'k--')

                % plot mean values for models
                for ind_model=1:length(models)
                    plot(ind_model+0.15,nanmean(nanmean(models(ind_model).mean_dist_img(:,1:end-1),2)),'o','Color',models(ind_model).cols)
                end
            hold off
            ylabel('Average Multimatch')
            set(gca,'XLim',[-1 length(models)+1],'XTickLabel',[' ', 'Subjects', short_names, ' '])
            pause(0.1)
         set(gcf,'name', strcat(fig_preffix_2, 'mean_mm'))


%% Save seachers figs

if save_plots
    disp('guardando figuras...')
    FolderFigsName = './../figs/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    SAVEPATH = strcat(FolderFigsName, path_to_save, 'supp/');
    if ~exist(SAVEPATH,'dir') 
        mkdir(SAVEPATH);
    end
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, strcat(SAVEPATH, FigName, '.fig'));
      %saveas(FigHandle, strcat(SAVEPATH, FigName, '.svg'));
      print(FigHandle,strcat(SAVEPATH, FigName,'.png'),'-dpng','-r400')
    end
    close all
end
 
%% AVERAGE MULTIMATCH Boxplot ALL multimatch models comparison

figure(421); 
    group_mm_average = [];
    for mod = 1:length(models)
        group_mm_average = [group_mm_average nanmean(models(mod).mean_dist_img(:,1:4),2)];
    end

    boxplot(group_mm_average, 'Color','k','Notch','on','Labels', short_names) 
    title('Average Multimatch')

    h = findobj(gca,'Tag','Box');
    colors = {models.cols};
    % findobj gets boxplot reversed so I reverse colors with 5-j
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors{length(models)+1-j},'FaceAlpha',.5)
    end
    
    % add humans data
    hold on
        % filter nan images
        mean_mm_distance_img = mean_dist_img(:, 1:end-1);
        % drop images with no views (4 and 132)
        mean_mm_distance_img([4 132],:) = [];
        % mean_mm_distance_img_nonzero = mean_mm_distance_img(mean_mm_distance_img~=0);
        % calculate subject stats
        betweenh_mean = mean(mean(mean_mm_distance_img,2));
        betweenh_std = std(mean(mean_mm_distance_img,2));

        yline(betweenh_mean,'b--')
        yline(betweenh_mean+betweenh_std,'b-.')
        yline(betweenh_mean-betweenh_std,'b-.')
    hold off
% 
    set(gcf,'name', strcat(fig_preffix_1, 'boxplots'))
%% AVERAGE MULTIMATCH SCATTER ALL multimatch models comparison

for mod = 1:length(models)
    x = group_mm_average(:,mod);
    y = nanmean(mean_dist_img(:,1:4),2);
    figure(430+mod)
    line('Color','k','LineStyle','--')
    hold on
    sf = scatter(x,y,[], models(mod).cols, 'filled', 'MarkerFaceAlpha', 0.65);
%     xlim([0.6 1])
%     ylim([0.6 1])
    xlim([0.7 0.95])
    ylim([0.7 0.95])
    % Correlation and linear model fit
    not_nan = isfinite(x) & isfinite(y); 
    [models(mod).corr, models(mod).pvalue, models(mod).modclb, models(mod).icub] = ...
        corrcoef(x(not_nan),y(not_nan));
    [models(mod).corr_spearman, models(mod).pvalue_spearman] = ...
        corr(x(not_nan),y(not_nan),'Type','Spearman');
    mdl = fitlm(x,y,'Intercept',false);
    models(mod).linearregression = mdl;
    hold on
    h = plot(mdl);
    legend off
    delete(h(1));
%     legend([h sf], {'', '', '','' ,models(mod).name})
    xlabel(models(mod).name)
    set(gcf,'name', strcat(fig_preffix_1, short_names{mod}))
    f = gcf;
    f.Position(3:4) = [500 650];
    pause(0.2)
end

%% Save figs

if save_plots
    disp('guardando figuras...')
    FolderFigsName = './../figs/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    SAVEPATH = strcat(FolderFigsName, path_to_save, 'fig4/');
    if ~exist(SAVEPATH,'dir') 
        mkdir(SAVEPATH);
    end
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, strcat(SAVEPATH, FigName, '.fig'));
      %saveas(FigHandle, strcat(SAVEPATH, FigName, '.svg'));
      print(FigHandle,strcat(SAVEPATH, FigName,'.png'),'-dpng','-r400')
    end
    close all
end