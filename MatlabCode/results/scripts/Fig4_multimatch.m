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

% graph_type = 'priors';
% graph_type = 'ibs';
% graph_type = 'searchers';
graph_type = 'ssim';
save_plots = 0;

switch graph_type
    case 'searchers'
        load('mm_hm_reduced_5searchers.mat')
        short_names    = {'cIBS+DGII','sIBS+DGII','IBS+DGII','Greedy+DGII','SaliencyB.'};
        fig_preffix    = 'Fig4_5searcher_boxplot_';
        fig_preffix_1  = 'Fig4_5searcher_mm_all_';
        fig_preffix_2  = 'Fig4_5searcher_scatter_';
        path_to_save   = 'searchers/';
    case 'priors'
        load('mm_hm_reduced_priors.mat')
        short_names    = {'cIBS+DGII','cIBS+Center','cIBS+Flat','cIBS+Noisy'};
        fig_preffix    = 'Fig4_priors_boxplot_';
        fig_preffix_1  = 'Fig4_prior_mm_all_';
        fig_preffix_2  = 'Fig4_5searcher_scatter_';
        path_to_save   = 'priors/';
    case 'ssim'
        load('mm_hm_reduced_ssim.mat')
        short_names    = {'sIBS+DGII','sIBS+Center','sIBS+Noisy','sIBS+Flat'};
        fig_preffix    = 'Fig4_prior_ssim_boxplot_';
        fig_preffix_1  = 'Fig4_prior_ssim_mm_all_';
        fig_preffix_2  = 'Fig4_prior_ssim_scatter_';
        path_to_save   = 'priors/';
    case 'ibs'
        %
        %
        %
end

%% Calculate subjects statistics - REVISAR

ind_model = 1;
mm_mean   = nanmean(models(ind_model).mean_dist_img);
mm_std    = nanstd(models(ind_model).mean_dist_img);
gold_standard = mean(mean_dist_img);

% mean_mm_distance_img = mean_dist_img(:, ind_metric);
% % drop images with no views (4 and 132)
% mean_mm_distance_img_nonzero = mean_mm_distance_img(mean_mm_distance_img~=0);
% % calculate subject stats
% mean_mm_distance = mean(mean_mm_distance_img_nonzero);
% std_mm_distance = std(mean_mm_distance_img_nonzero);

metric_ind = 1;
x = models(ind_model).mean_dist_img(:,metric_ind);
y = mean_dist_img(:,metric_ind);
not_nan = isfinite(x) & isfinite(y); 
% figure(1)
%     scatter(x,y,[], models(ind_model).cols, 'filled', 'MarkerFaceAlpha',0.6);
%     xlim([0.75 1])
%     ylim([0.75 1])
%     refline(1,0)
    
%% BOXPLOT Plot multimatch models comparison - searchers

%close all
for mm = 1:length(mm_names)
    figure(mm+400); 
%         group = [models(1).(mm_names{mm}),... 
%             models(2).(mm_names{mm}),...
%             models(3).(mm_names{mm}),...
%             models(4).(mm_names{mm})];
%         group = zeros(Nimg,length(models));
        group = [];
        for mod = 1:length(models)
            group = [group models(mod).(mm_names{mm})];
        end
%         boxplot(group,'Notch','on','Labels', short_names,'color', [models.cols])   
        boxplot(group, 'Color','k','Notch','on','Labels', short_names) 
        title(mm_titles{mm})
        ylabel(mm_names{mm})
%        set(gca,'XTickLabel',short_names)
        
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

%% Save seachers figs

if save_plots
    disp('guardando figuras...')
    FolderName = './../figs/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, strcat(FolderName, FigName, '.fig'));
      saveas(FigHandle, strcat(FolderName, FigName, '.svg'));
%       saveas(FigHandle, strcat(FolderName, FigName, '.png'));
      print(FigHandle,strcat(FolderName, FigName,'.png'),'-dpng','-r300')
    end
end
close all

%% SCATTERPLOT Plot multimatch models comparison - searchers

for ind_metric=1:4
    figure(410+ind_metric)
    set(gcf,'Renderer', 'painters', 'Position', [100 100 800 450])
    set(gcf,'Color','w')
        % consider only trials with ...(completar)       
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
    %         set(gca,'YTick',0:0.1:0.4)
            ylabel(mm_names(ind_metric))
            %set(gca,'XLim',[-1 length(models)+1],'XTick',0:length(models),'XTickLabel',['humans' {models.name}])
            %set(gca,'XTickLabelRotation',-45)
            set(gca,'XLim',[-1 length(models)+1],'XTickLabel',[' ', 'Subjects', short_names, ' '])
    %         set(gca,'YLim',[0.6,1])
    %         xlabel('Models')
            %[hleg1, icons] = legend(['humans' {models.name}],'Location', 'Best');
            %[hleg1, icons] = legend({'humans','cIBS+DeepGaze2','cIBS+Center','cIBS+Flat','cIBS+Noisy'},'Location', 'Best');

            %[hleg1, icons] = legend({'humans','cIBS+DeepGaze2','IBS+DeepGaze2','Greedy+DeepGaze2','SaliencyBased'},'Location', 'Best');
    %         [hleg1, icons] = legend(['humans' short_names],'Location', 'Best');

    %         icon = findobj(icons,'Type','line');
    %         icon = findobj(icon,'Marker','none','-xor');
    %         set(icon,'MarkerSize',30);
            pause(0.1)
    %         textobj = findobj(icons, 'type', 'text');
    %         set(textobj, 'fontsize', 12);
          set(gcf,'name', strcat(fig_preffix_2, mm_names{ind_metric}))
%         ax = gca;
%         outerpos = ax.OuterPosition;
%         ti = ax.TightInset; 
%         left = outerpos(1) + ti(1);
%         bottom = outerpos(2) + ti(2);
%         ax_width = outerpos(3) - ti(1) - ti(3);
%         ax_height = outerpos(4) - ti(2) - ti(4);
%         ax.Position = [left bottom ax_width ax_height];
%         pause(0.1)    
    %     textobj = findobj(icons, 'type', 'text');
    %     set(textobj, 'fontsize', 12);
end

%% all mm mean

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
    %         set(gca,'YTick',0:0.1:0.4)
            ylabel('Average Multimatch')
            %set(gca,'XLim',[-1 length(models)+1],'XTick',0:length(models),'XTickLabel',['humans' {models.name}])
            %set(gca,'XTickLabelRotation',-45)
            set(gca,'XLim',[-1 length(models)+1],'XTickLabel',[' ', 'Subjects', short_names, ' '])
    %         set(gca,'YLim',[0.6,1])
    %         xlabel('Models')
            %[hleg1, icons] = legend(['humans' {models.name}],'Location', 'Best');
            %[hleg1, icons] = legend({'humans','cIBS+DeepGaze2','cIBS+Center','cIBS+Flat','cIBS+Noisy'},'Location', 'Best');

            %[hleg1, icons] = legend({'humans','cIBS+DeepGaze2','IBS+DeepGaze2','Greedy+DeepGaze2','SaliencyBased'},'Location', 'Best');
    %         [hleg1, icons] = legend(['humans' short_names],'Location', 'Best');

    %         icon = findobj(icons,'Type','line');
    %         icon = findobj(icon,'Marker','none','-xor');
    %         set(icon,'MarkerSize',30);
            pause(0.1)
    %         textobj = findobj(icons, 'type', 'text');
    %         set(textobj, 'fontsize', 12);
          set(gcf,'name', strcat(fig_preffix_2, 'mean_mm'))
%         ax = gca;
%         outerpos = ax.OuterPosition;
%         ti = ax.TightInset; 
%         left = outerpos(1) + ti(1);
%         bottom = outerpos(2) + ti(2);
%         ax_width = outerpos(3) - ti(1) - ti(3);
%         ax_height = outerpos(4) - ti(2) - ti(4);
%         ax.Position = [left bottom ax_width ax_height];
%         pause(0.1)    
    %     textobj = findobj(icons, 'type', 'text');
    %     set(textobj, 'fontsize', 12);

%% Save seachers figs
save_plots = 0;
if save_plots
    disp('guardando figuras...')
    FolderName =  '/home/gastonb/Imágenes/figpaper/supp1/';   % './../figs/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
%       savefig(FigHandle, strcat(FolderName, FigName, '.fig'));
      saveas(FigHandle, strcat(FolderName, FigName, '.svg'));
      print(FigHandle,strcat(FolderName, FigName,'.png'),'-dpng','-r400')
    end
end
close all
 
%% AVERAGE MULTIMATCH Boxplot and scatter ALL multimatch models comparison - searchers

close all
figure(421); 
    group_mm_average = [];
    for mod = 1:length(models)
        group_mm_average = [group_mm_average nanmean(models(mod).mean_dist_img(:,1:4),2)];
    end

    boxplot(group_mm_average, 'Color','k','Notch','on','Labels', short_names) 
    title('Average Multimatch')
%     ylabel(mm_names{mm})
%     set(gca,'XTickLabel',short_names)

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
%%
close all
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
%     legend(sf, models(mod).name)
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

%%

save_plots = 1;
if save_plots
    disp('guardando figuras...')
    FolderName =  strcat('/home/gastonb/Imágenes/figpaper/fig4/', path_to_save);   % './../figs/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
%       savefig(FigHandle, strcat(FolderName, FigName, '.fig'));
      saveas(FigHandle, strcat(FolderName, FigName, '.svg'));
      print(FigHandle,strcat(FolderName, FigName,'.png'),'-dpng','-r400')
    end
end
close all

%% AVERAGE MULTIMATCH - Corrplot Searchers mean/median multimatch

figure(450); gcf
    set(gcf,'Color','w')   
%     set(gcf,'Position',[565 70 430 970])
    y = nanmedian(mean_dist_img, 2);
    % layout
    l1 = 0.09;
    l2 = 0.55;
    b1 = 0.58;
    b2 = 0.088;
    s = 0.39;
    positions = [l1 b1 s s;...
        l2 b1 s s; ...
        l1 b2 s s; ...
        l2 b2 s s];
    for i=1:length(models)
        ha(i) = axes('Position',positions(i,:)); 
        %subplot(2,2,i)
            line('Color','k','LineStyle','--')
            hold on
            x = nanmedian(models(i).mean_dist_img,2);
            scatter(x,y,[], models(i).cols, 'filled', 'MarkerFaceAlpha',...
                0.6);
            xlim([0 0.40])
            ylim([0 0.3])
            % Correlation and linear model fit
            not_nan = isfinite(x) & isfinite(y); 
            [models(i).corr, models(i).pvalue, models(i).iclb, models(i).icub] = ...
                corrcoef(x(not_nan),y(not_nan));
            [models(i).corr_spearman, models(i).pvalue_spearman] = ...
                corr(x(not_nan),y(not_nan),'Type','Spearman');
            mdl = fitlm(x,y,'Intercept',false);
            models(i).linearregression = mdl;
            hold on
            h = plot(mdl);
            delete(h(1));
            h(2).LineStyle = '-';
            h(2).LineWidth = 0.8;
            h(3).LineStyle = '-.';
            h(3).LineWidth = 0.6;
            h(4).LineStyle = '-.';
            h(4).LineWidth = 0.6;
            legend off;
            xlabel(models(i).name,'FontSize',11)
%             if i==1 || i==3
%                 ylabel('Humans','FontSize',11)
%             end
%             if i==2 || i==4
%                 ylabel('')
%             end
%             title('')
            %refline([1,0],'Color','k')
%             mat = [(1:134); x' ; y';abs(x-y)'];
%             err = sortrows(mat',4, 'descend');
%             % Dropnans
%             err(any(isnan(err), 2), :) = [];
%             % elijo las primeras tres y las ploteo
%             im = err(1:3,1)';
%             xt = err(1:3,2)' + 0.003;
%             yt = err(1:3,3)' + 0.003;
%             %str = {arrayfun(@num2str, im', 'UniformOutput', false)};
%             str = {strcat('\leftarrow',num2str(im(1))), ...
%                     strcat('\leftarrow',num2str(im(2))), ...
%                     strcat('\leftarrow',num2str(im(3)))};
            %text(xt,yt,str)
            hold off
    end
    
%% Load priors

load('mm_hm_reduced_priors.mat')
short_names = {'cIBS+DGII','cIBS+Center','cIBS+Flat','cIBS+Noisy'};
fig_preffix = 'Fig4_prior_mm_';
save_priors = false;

%% Plot multimatch models comparison boxplot - priors

for mm = 1:length(mm_names)
    figure(470+mm)
        group = [models(1).(mm_names{mm}),... 
            models(2).(mm_names{mm}),...
            models(3).(mm_names{mm}),...
            models(4).(mm_names{mm})];

        boxplot(group,'Notch','on','Color','k')   
        set(gca,'XTickLabel',short_names)
        title(mm_titles{mm})
        ylabel(mm_names{mm})
        
        h = findobj(gca,'Tag','Box');
        colors = {models.cols};
        % findobj gets boxplot reversed so I reverse colors with 5-j
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors{5-j},'FaceAlpha',.5)
        end
        
        % add humans data
        hold on
        betweenh_mean = nanmean(multimatch_bh.(mm_names{mm}));
        yline(betweenh_mean,'b--')
        hold off
        
        set(gcf,'name', strcat(fig_preffix,mm_names{mm}))
    pause(0.05);
end

%% Save priors figs

if save_priors
    FolderName = './../figs/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, strcat(FolderName, FigName, '.fig'));
      saveas(FigHandle, strcat(FolderName, FigName, '.svg'));
    end
    % close all
end