clear all
close all
clc

addpath('~/Dropbox/my_functions/')
addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
% load('../matrix/info_per_subj_additional_data.mat'); % loads images, img_order, templates, template_order
% addpath(genpath('../static_models/utils/'));
% 
% load('../matrix/initial_fixations.mat');

load('../matrix/info_per_subj_final.mat');
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);

a=dir('../matrix/images/*.mat');    filenames_img = {a.name}'; Nimg = length(filenames_img);
a=dir('../matrix/subjects/*.mat');  filenames_subj = {a.name};

% models = fun_define_models;
models = fun_define_models_tmp(3);

%% Human data
nntrthr     = 20;

% Average number of fixations needed by the observers to find the target
Nfix_img_mean   = nan(Nimg,1);
Nfix_img_std    = nan(Nimg,1);
Nfix_img_nsuj   = nan(Nimg,1);
for ind_img=1:Nimg % images
    path = char(strcat('../matrix/images/info_per_subj_img_', num2str(ind_img), '.mat')); fprintf('%s\n',path);
    load(path);
    if (ind_img ~= 132)
        founds = [info_per_subj.target_found];
        fixs = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_subj);
        Nfix_img_mean(ind_img)  = mean(fixs(founds));
        Nfix_img_std(ind_img)   = std(fixs(founds));
        Nfix_img_nsuj(ind_img)  = sum(founds);
    end
end
nntrfilt    = (Nfix_img_nsuj>nntrthr);

% Proportion of targets found as function of the number of saccades allowed for the humans
mean_nsac_target_found_img  = nan(Nsubj,4); % The number of actual saccades on the grid (could be less than the allowed)
target_found_img            = nan(Nsubj,4);
target_found_img_ntrials    = nan(Nsubj,4);
for ind_subj=[1:43 45:Nsubj]
    path = char(strcat('../matrix/subjects/info_per_subj_subj_', num2str(ind_subj), '.mat')); fprintf('%s\n',path);
    load(path);
    
    founds              = [info_per_subj.target_found];
%     fixspp              = arrayfun(@(x) length(x.fixations), info_per_subj);
%     nsaccades           = arrayfun(@(x) x.exp_data.nsaccades, info_per_subj);
    nfixs               = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_subj);
    nsaccades_allowed   = arrayfun(@(x) x.exp_data.nsaccades_thr, info_per_subj);
    NSACCADES_allowed   = [2 4 8 12]; %unique(nsaccades_allowed);
    for isacc=1:length(NSACCADES_allowed)
        target_found_img(ind_subj,isacc)            = sum(founds(nsaccades_allowed==NSACCADES_allowed(isacc)));
        target_found_img_ntrials(ind_subj,isacc)    = sum(nsaccades_allowed==NSACCADES_allowed(isacc));
        mean_nsac_target_found_img(ind_subj,isacc)  = mean(nfixs(nsaccades_allowed==NSACCADES_allowed(isacc)));
    end
    if size(target_found_img,2)>4; break; end
end
P_target_found = target_found_img./target_found_img_ntrials;

% Proportion of targets found as function of the number of saccades for the humans
target_found_img_accum            = nan(Nsubj,13);
target_found_img_accum_ntrials    = nan(Nsubj,13);
for ind_subj=[1:43 45:Nsubj]
    path = char(strcat('../matrix/subjects/info_per_subj_subj_', num2str(ind_subj), '.mat')); fprintf('%s\n',path);
    load(path);
    
    founds              = [info_per_subj.target_found];
    nfixs               = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_subj);
    for ifixs=1:13
        target_found_img_accum(ind_subj,ifixs)        = sum(founds(nfixs==ifixs));
        target_found_img_accum_ntrials(ind_subj,ifixs)= sum(nfixs==ifixs);
    end
%     if size(target_found_img_accum,2)>4; break; end
end
P_target_found_accum = target_found_img_accum./target_found_img_accum_ntrials;

figure(100); clf
    set(gcf,'Color','w')
    boxplot(P_target_found_accum,'Notch','on','Color','k')
    

%% NOTA IMPORTANTE: TIRAR SUJETO 44!!! TIENE 0 FIJACIONES EN TODOS LOS TRIALS
% ind_subj=44;
% path = char(strcat('../matrix/subjects/info_per_subj_subj_', num2str(ind_subj), '.mat')); fprintf('%s\n',path);
% load(path);
% arrayfun(@(x) length(x.fixations), info_per_subj)
% 
% ver chequeo_suj44_MVA.m

%% NOTA IMPORTANTE: QUE TENIA LA IMAGEN 132??

%% NOTA IMPORTANTE: DECIDIR SI COMPARAMOS NUMERO DE FIJACIONES TOTALES, O REDUCIDAS A LA GRILLA DEL MISMO TAMAÑO QUE EL MODELO
% Estimar cuanto colapsamos y eventualmente se pueden graficar los dos.

%% Model data
% Number of fixations predicted by the models
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
for ind_model=1:length(models)
    target_found_img_model     = zeros(Nimg,4);
    target_notfound_img_model  = zeros(Nimg,4);
    for ind_img=1:Nimg
        path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                                    models(ind_model).prior,...
                                    models(ind_model).searcher,...
                                    models(ind_model).params,...
                                    ind_img);
        fprintf('%s\n',path);    
        if ind_img ~= 132
            load(path); % scanpath = ( (x,y) , fixation number )
            
            count       = length(scanpath) - 1;     % Numero de sacadas
            countthr    = [2 4 8 12];               % Numero de sacadas permitidas (igual que el experimento).
            
            
            
            for j=1:length(countthr)
                if (count <= countthr(j))   % Si el modelo lo encontro en menos, lo encontro!!
                    target_found_img_model(ind_img,j)     = target_found_img_model(ind_img,j) + 1;
                end
            end            
        end        
    end
    models(ind_model).target_found_img_model    = target_found_img_model;
    models(ind_model).P_target_found_model      = mean(target_found_img_model([1:131 133:134],:));
    models(ind_model).P_target_found_model_std  = std(target_found_img_model([1:131 133:134],:));
end

%% Scanpath similarity between humans
delta       = 32;
min_fix     = 2;
max_fix     = 13;
image_size  = [768 1024];
grid_size   = image_size/delta;
addpath('../compare_models')
addpath('../dynamic_models/utils')
    
mean_dist_img   = nan(Nimg,Nsubj);
std_dist_img    = nan(Nimg,Nsubj);

eliminados_subj = nan(Nimg,1);
adentro_subj    = nan(Nimg,1);
    
for ind_img=1:Nimg
    if ind_img ~= 132
        path = char(strcat('../matrix/images/info_per_subj_img_', num2str(ind_img), '.mat')); fprintf('%s\n',path);
        load(path)
        
        % (2019-09-18) Map fixations into the Grid
        % Update fixations_matrix_reduced, timefix_matrix_reduced, mascara_matrix_reduced and delta
        % fields in the info_per_subj structure.
        info_per_subj = subjMapFixationToMatrix( info_per_subj, path, delta, image_size );
        % Calculo la distancia de todos contra todos para los sujetos que
        % tengan min_fix < #fix < max_fix y que ambos hayan encontrado el
        % target 
        
        % Calculate the distance between subjects for each image
        ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                            (length(x.fixations_matrix_reduced(:,1)) > min_fix && ...
                                            length(x.fixations_matrix_reduced(:,1)) < max_fix) ), info_per_subj);
        distance = nan(Nsubj,Nsubj);
        for subj_i=1:length(info_per_subj)
            if ind_subj_target_found(subj_i)
                for subj_j=(subj_i+1):length(info_per_subj) 
                    if ind_subj_target_found(subj_j)
                        distance(subj_i, subj_j) = scanpathDistance(info_per_subj(subj_i).fixations_matrix_reduced,...
                                                                    info_per_subj(subj_j).fixations_matrix_reduced,...
                                                                    grid_size);
                    end
                end
            end
        end
        
        adentro_subj(ind_img)   = sum(ind_subj_target_found);
        eliminados_subj(ind_img)= sum(~ind_subj_target_found & [info_per_subj.target_found]);
        
        % Calculo la distancia promedio de todos contra todos
        full_distance = [distance distance'];
        
        %guardo la distancia por imagen 
        mean_dist_img(ind_img,:)    = nanmean(full_distance,2)';
        std_dist_img(ind_img,:)     = nanstd(full_distance,0,2)';
    end
end

%% Scanpath similarity between humans and model
delta       = 32;
min_fix     = 2;
max_fix     = 13;
image_size  = [768 1024];
grid_size   = image_size/delta;
addpath('../compare_models')
addpath('../dynamic_models/utils')
    
for ind_model=1:length(models)
    models(ind_model).mean_dist_img   = nan(Nimg,Nsubj);
    models(ind_model).std_dist_img    = nan(Nimg,Nsubj);

    models(ind_model).eliminados_subj = nan(Nimg,1);
    models(ind_model).adentro_subj    = nan(Nimg,1);
end

for ind_img=1:Nimg
    if ind_img ~= 132
        path = char(strcat('../matrix/images/info_per_subj_img_', num2str(ind_img), '.mat')); fprintf('%s\n',path);
        load(path)

        % (2019-09-18) Map fixations into the Grid
        % Update fixations_matrix_reduced, timefix_matrix_reduced, mascara_matrix_reduced and delta
        % fields in the info_per_subj structure.
        info_per_subj = subjMapFixationToMatrix( info_per_subj, path, delta, image_size );
        % Calculo la distancia de todos contra todos para los sujetos que
        % tengan min_fix < #fix < max_fix y que ambos hayan encontrado el
        % target 

        % Calculate the distance between subjects for each image
        ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                            (length(x.fixations_matrix_reduced(:,1)) > min_fix && ...
                                            length(x.fixations_matrix_reduced(:,1)) < max_fix) ), info_per_subj);

        for ind_model=1:length(models)
            path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                            models(ind_model).prior,...
                            models(ind_model).searcher,...
                            models(ind_model).params,...
                            ind_img);
            load(path)

            if (size(scanpath,1) < max_fix)
                distance = nan(Nsubj,1);
                for subj_i=1:length(info_per_subj)
                    if ind_subj_target_found(subj_i)
                        distance(subj_i) = scanpathDistance(info_per_subj(subj_i).fixations_matrix_reduced,...
                                                                    scanpath,...
                                                                    grid_size);
                    end
                end

                models(ind_model).adentro_subj(ind_img)   = sum(ind_subj_target_found);
                models(ind_model).eliminados_subj(ind_img)= sum(~ind_subj_target_found & [info_per_subj.target_found]);

                %guardo la distancia por imagen 
                models(ind_model).mean_dist_img(ind_img,:)    = nanmean(distance)';
                models(ind_model).std_dist_img(ind_img,:)     = nanstd(distance)';
            end
        end
    end
end

%% Figure 4
% Fig 4A: Proportion of targets found for each model (lines) and humans (boxplot)

ha=[];
figure(2); clf
    set(gcf,'Color','w')
    set(gcf,'Position',[565 70 430 970])
%     ha(1)=axes('Position',[0.075 0.650 0.300 0.300]); 
    subplot(3,1,1)
        hold on
            boxplot(P_target_found,'Notch','on','Color','k','Labels',[2,4,8,12])
            for ind_model=1:length(models)
                [~,~,~,Pc] = fun_evaluate_experiment(models(ind_model).Nfix_img_model, 57);
                plot(nanmean(Pc'),'.-','Color',models(ind_model).cols,'LineWidth',1.5)
%                 plot(models(ind_model).P_target_found_model,'.-','Color',models(ind_model).cols,'LineWidth',1.5)
            end
        hold off

        set(gca,'XTickLabels',{'2','4','8','12'})

        xlabel('Number of max saccades')
        ylabel('Proportion of targets found')
        ylim([0 1])
        
% Fig 4B: Number of fixations needed to find the target (models and humans)
%     ha(2)=axes('Position',[0.575 0.650 0.300 0.300]); 
    subplot(3,1,2)
        hold on
            x = 1:15;
            y = hist(Nfix_img_mean,x); y = y/sum(y);
            plot(x,y,'k-','LineWidth',1.5)
            for ind_model=1:length(models)
                [~,~,NFix,~] = fun_evaluate_experiment(models(ind_model).Nfix_img_model, 57);
                y = hist(nanmean(NFix'),x); y = y/nansum(y); 
%                 y = hist(models(ind_model).Nfix_img_model,x); y = y/nansum(y); 
                plot(x,y,'-','Color',models(ind_model).cols,'LineWidth',1.5)
            end
        hold off
        box on
%         set(gca,'YLim',[0 0.5])
        set(gca,'YTick',0:0.1:0.4)
        set(gca,'XLim',[0 13])
        xlabel('Number of saccades (to find the target)')
        ylabel('Frequency')

    subplot(3,1,3)
        Nimg_thr = 25;
        hold on
            filter_img  = adentro_subj>Nimg_thr;
            Nimg_ok     = sum(filter_img);
            plot(zeros(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(mean_dist_img(filter_img),2),'k.')
%             plot(zeros(Nimg,1) + (2*rand(Nimg,1)-1)*0.1,nanmean(mean_dist_img,2),'k.')
            for ind_model=1:length(models)
                filter_img  = models(ind_model).adentro_subj>Nimg_thr;
                Nimg_ok     = sum(filter_img);
                plot(ind_model*ones(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(models(ind_model).mean_dist_img(filter_img),2),'.','Color',models(ind_model).cols)
%                 plot(ind_model*ones(Nimg,1) + (2*rand(Nimg,1)-1)*0.1,nanmean(models(ind_model).mean_dist_img,2),'.','Color',models(ind_model).cols)
            end
            plot([-1 length(models)+1],nanmean(nanmean(mean_dist_img,2))*[1 1],'k-')
            plot([-1 length(models)+1],(nanmean(nanmean(mean_dist_img,2)) + nanstd(nanmean(mean_dist_img,2)))*[1 1],'k--')
            plot([-1 length(models)+1],(nanmean(nanmean(mean_dist_img,2)) - nanstd(nanmean(mean_dist_img,2)))*[1 1],'k--')
            for ind_model=1:length(models)
                plot(ind_model+0.15,nanmean(nanmean(models(ind_model).mean_dist_img,2)),'o','Color',models(ind_model).cols)
            end
        hold off
        set(gca,'YTick',0:0.1:0.4)
        ylabel('Scanpath Similarity')
        set(gca,'XLim',[-1 length(models)+1],'XTick',0:length(models),'XTickLabel',['humans' {models.name}])
        set(gca,'XTickLabelRotation',-45)
        
%% Figure 4
% % Fig 4A: Proportion of targets found for each model (lines) and humans (boxplot)
% ha=[];
% figure(1); clf
%     set(gcf,'Color','w')
%     set(gcf,'Position',[565 70 795 775])
% %     ha(1)=axes('Position',[0.075 0.650 0.300 0.300]); 
%     subplot(2,2,1)
%         hold on
%             boxplot(P_target_found,'Notch','on','Color','k','Labels',[2,4,8,12])
%             for ind_model=1:length(models)
%                 plot(models(ind_model).P_target_found_model,'.-','Color',models(ind_model).cols,'LineWidth',1.5)
%             end
%         hold off
% 
%         set(gca,'XTickLabels',{'2','4','8','12'})
% 
%         xlabel('Number of max saccades')
%         ylabel('Proportion of targets found')
%         ylim([0 1])
%         
% % Fig 4B: Number of fixations needed to find the target (models and humans)
% %     ha(2)=axes('Position',[0.575 0.650 0.300 0.300]); 
%     subplot(2,2,2)
%         hold on
%             [y,x] = hist(Nfix_img_mean,1:15);
%             plot(x,y,'k-','LineWidth',1.5)
%             for ind_model=1:length(models)
%                 [y,x] = hist(models(ind_model).Nfix_img_model,1:15);
%                 plot(x,y,'-','Color',models(ind_model).cols,'LineWidth',1.5)
%             end
%         hold off
%         box on
%         set(gca,'YLim',[0 40])
%         set(gca,'XLim',[0 15])
%         xlabel('Number of saccades (to find the target)')
%         ylabel('Counts')
% 
%     subplot(2,2,4)
%         hold on
%             title('Number of saccades (to find the target)')
%             for ind_model=1:length(models)
%                 plot(Nfix_img_mean, models(ind_model).Nfix_img_model+(2*rand(length(models(ind_model).Nfix_img_model),1)-1)*0.2, ...
%                 '.','Color',models(ind_model).cols,'LineWidth',1.5);
% %                 lsline();
%             end
%             plot([0 15],[12 12],'k-')
%             plot([12 12],[0 15],'k-')
%         hold off
%         box on
%         set(gca,'YLim',[0 15])
%         set(gca,'XLim',[0 15])
%         xlabel('humans')
%         ylabel('models')
%         
%     subplot(2,2,3)
%         hold on
%             plot(zeros(Nimg,1) + (2*rand(Nimg,1)-1)*0.1,nanmean(mean_dist_img,2),'k.')
%             for ind_model=1:length(models)
%                 plot(ind_model*ones(Nimg,1) + (2*rand(Nimg,1)-1)*0.1,nanmean(models(ind_model).mean_dist_img,2),'.','Color',models(ind_model).cols)
%             end
%             plot([-1 length(models)+1],nanmean(nanmean(mean_dist_img,2))*[1 1],'k-')
%             plot([-1 length(models)+1],(nanmean(nanmean(mean_dist_img,2)) + nanstd(nanmean(mean_dist_img,2)))*[1 1],'k--')
%             plot([-1 length(models)+1],(nanmean(nanmean(mean_dist_img,2)) - nanstd(nanmean(mean_dist_img,2)))*[1 1],'k--')
%             for ind_model=1:length(models)
%                 plot(ind_model+0.15,nanmean(nanmean(models(ind_model).mean_dist_img,2)),'o','Color',models(ind_model).cols)
%             end
%         hold off
%         set(gca,'XLim',[-1 length(models)+1],'XTick',0:length(models),'XTickLabel',['humans' {models.name}])
%         set(gca,'XTickLabelRotation',45)