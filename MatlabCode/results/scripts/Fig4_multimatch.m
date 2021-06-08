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

trials_tmp                      = load(strcat(src_path, 'info_all_subj.mat'));
[ids_images, ~, images_order]   = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                    = images_order';
image_size                      = trials_tmp.info_per_subj_final(1).image_size;
delta                           = 32;
grid_size                       = image_size/delta;
trials                          = reduce_scanpaths(trials_tmp.info_per_subj_final, delta, image_size);
clear trials_tmp

%% Load results

addpath('./results_metrics/')
load('mm_bh_reduced.mat')
load('mm_hm_reduced_searchers.mat')

%% A
ind_model = 1;
mm_mean   = nanmean(models(ind_model).mean_dist_img);
mm_std    = nanstd(models(ind_model).mean_dist_img);
gold_standard = mean(mean_dist_img);

%% Plot

metric_ind = 1;
metric_labels = {'','','','',''};
x = models(ind_model).mean_dist_img(:,metric_ind);
y = mean_dist_img(:,metric_ind);
not_nan = isfinite(x) & isfinite(y); 
figure(1)
    scatter(x,y,[], models(ind_model).cols, 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0.75 1])
    ylim([0.75 1])
    refline(1,0)

%% all searchers boxplot

metric_ind = 1;
figure(2)
    for ind_m=1:length(models)
        % var to plot
        var = models(ind_m).mean_dist_img(:,metric_ind);
        % plot options
        set(gcf,'Color','w')
        x0=10;
        y0=10;
        width=700;
        height=500;
        set(gcf,'position',[x0,y0,width,height])
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        
        hold on
            boxplot(var,'Notch','on','Color','k')
%             for ind_model=1:length(models)
%                 [correct, M_f, N_f,Pc] = fun_evaluate_experiment(models(ind_model).Nfix_img_model, 57);
%                 plot(nanmean(Pc'),'.-','Color',models(ind_model).cols,'LineWidth',1.5)
% %                 plot(models(ind_model).P_target_found_model,'.-','Color',models(ind_model).cols,'LineWidth',1.5)
%             end
            
        hold off
    end
        names = {models.name};
        [hleg1, icons] = legend(names, 'Location', 'NorthWest', 'FontSize', 11);
        set(gca,'XTickLabels',models.name)
    %             icon = findobj(icons,'Type','line');
    %             icon = findobj(icon,'Marker','none','-xor');
    %             set(icon,'MarkerSize',30);
    %             pause(0.1)
    %             textobj = findobj(icons, 'type', 'text');
    %             set(textobj, 'fontsize', 12);
        xlabel('Number of max saccades') 

%%

figure(3)
    var = structfun(@(x) x(:,metric_ind),models.mean_dist_img);
