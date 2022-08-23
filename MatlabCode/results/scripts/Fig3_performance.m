% Figure 3 - Model Performance
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

% save metrics
guardar = 0;

trials_tmp                      = load(strcat(src_path, 'info_all_subj.mat'));
[ids_images, ~, images_order]   = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                    = images_order';
image_size                      = trials_tmp.info_per_subj_final(1).image_size;
delta                           = 32;
grid_size                       = image_size/delta;
trials                          = reduce_scanpaths(trials_tmp.info_per_subj_final, delta, image_size);
clear trials_tmp

%% Calculations: Average number of fixations needed by the observers to find the target

nntrthr     = 15;

Nfix_img_mean   = nan(Nimg,1);
Nfix_img_std    = nan(Nimg,1);
Nfix_img_nsuj   = nan(Nimg,1);

for ind_img=1:Nimg % images
    img_id = ids_images{ind_img};
    info_per_img = trials(strcmp({trials.image_name}, img_id));
    if (ind_img ~= 132)
        founds = [info_per_img.target_found];
        fixs = arrayfun(@(x) length(x.x_grid), info_per_img);
        Nfix_img_mean(ind_img)  = mean(fixs(founds));
        Nfix_img_std(ind_img)   = std(fixs(founds));
        Nfix_img_nsuj(ind_img)  = sum(founds);
    end
end

nntrfilt    = (Nfix_img_nsuj>nntrthr);

humans_data_img.fix_found_mean = Nfix_img_mean;
humans_data_img.fix_found_std  = Nfix_img_std;
humans_data_img.fix_found_prop = Nfix_img_nsuj;

%%  Proportion of targets found as function of the number of saccades allowed for the humans

mean_nsac_target_found_img  = nan(Nsubj,4); % The number of actual saccades on the grid (could be less than the allowed)
target_found_img            = nan(Nsubj,4);
target_found_img_ntrials    = nan(Nsubj,4);

for ind_subj=[1:43 45:Nsubj]
    subj_id = ids_subjs{ind_subj};
    info_per_subj = trials(strcmp({trials.subj}, subj_id));
    
    founds              = [info_per_subj.target_found];
    nfixs               = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_subj);
    nsaccades_allowed   = arrayfun(@(x) x.nsaccades_allowed, info_per_subj);
    NSACCADES_allowed   = [2 4 8 12]; %unique(nsaccades_allowed);
    for isacc=1:length(NSACCADES_allowed)
        target_found_img(ind_subj,isacc)            = sum(founds(nsaccades_allowed==NSACCADES_allowed(isacc)));
        target_found_img_ntrials(ind_subj,isacc)    = sum(nsaccades_allowed==NSACCADES_allowed(isacc));
        mean_nsac_target_found_img(ind_subj,isacc)  = mean(nfixs(nsaccades_allowed==NSACCADES_allowed(isacc)));
    end

    if size(target_found_img,2)>4; break; end
end

P_target_found = target_found_img./target_found_img_ntrials;

humans_data_subj.target_found           = target_found_img;
humans_data_subj.target_found_ntrials   = target_found_img_ntrials;
humans_data_subj.target_found_mean      = mean_nsac_target_found_img;   
humans_data_subj.proportion             = P_target_found;

%% Models data

% models = fun_define_models('priors-correlation');
% models = fun_define_models('searchers-deepgaze-ssim');
models = fun_define_models('priors-ssim');
% models = fun_define_models('priors-ibs');
% models = fun_define_models('searchers-deepgaze');

% Number of fixations predicted by the models
for ind_model = 1:length(models)
    Nfix_img_model = nan(Nimg,1);
    for ind_img=1:Nimg
        path_scanpath = sprintf('../results_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                        models(ind_model).prior,...
                        models(ind_model).searcher,...
                        models(ind_model).params,...
                        ind_img);
%         fprintf('%s\n',path);
        if ind_img ~= 132
            load(path_scanpath)
            Nfix_img_model(ind_img)  = size(scanpath,1) - 1;
        end
    end
    models(ind_model).Nfix_img_model = Nfix_img_model;
end

% Proportion of targets found as function of the number of saccades allowed for the models
exps_thr = nan(Nimg, Nsubj);
for ind_img=1:Nimg % images
    path = char(strcat(src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat')); 
%     fprintf('%s\n',path);
    info_per_img_aux = load(path);

    for  ind_subj = 1:length(info_per_img_aux.info_per_img)
        exps_thr(ind_img,ind_subj)   = info_per_img_aux.info_per_img(ind_subj).nsaccades_allowed;
    end
    %clear info_per_img_aux
end

% Calculo de targets encontrados MODELOS
% Gaston 7/4 - Creo que esto no se usa en otro lado
for ind_model=1:length(models)
    target_found_img_model     = zeros(Nimg,4);
%     target_notfound_img_model  = zeros(Nimg,4);
    for ind_img=1:Nimg
        path_scanpath = sprintf('../results_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                        models(ind_model).prior,...
                        models(ind_model).searcher,...
                        models(ind_model).params,...
                        ind_img);
%         fprintf('%s\n',path);    
        if ind_img ~= 132
            load(path_scanpath); % scanpath = ( (x,y) , fixation number )
            
            % No me convence esto... hacerlo bien. Tengo que empezar con las fijaciones de cero!!!            
            count       = length(scanpath);% - 1;     % Numero de sacadas
            tmpcountthr = exps_thr(ind_img,:);      % Numero de sacadas permitidas (igual que el experimento).
            tmpcountthr = tmpcountthr(~isnan(tmpcountthr));% Numero de sacadas permitidas (igual que el experimento).
            countthr    = [2 4 8 12];
            for j=1:length(countthr)
                tmp = (count <= tmpcountthr+ 1);
                target_found_img_model(ind_img,j) = mean( tmp( tmpcountthr==countthr(j) ) );
                % keyboard
            end
        end        
    end
    % ? aca hay casos donde el target_found_img_model([1:131 133:134],:)
    % tiene nans y me esta matando las medias y los desvios
    models(ind_model).target_found_img_model    = target_found_img_model;
    models(ind_model).P_target_found_model      = nanmean(target_found_img_model([1:131 133:134],:));
    models(ind_model).P_target_found_model_std  = nanstd(target_found_img_model([1:131 133:134],:));
    % keyboard
end

for ind_model=1:length(models)
    models(ind_model).target_found_distances_new = mean(((nanmean(humans_data_subj.proportion, 1) - ...
        models(ind_model).P_target_found_model) ./ ...
        nanstd(humans_data_subj.proportion, 0, 1)) .^2);
    models(ind_model).target_found_distances_new_std = std(((nanmean(humans_data_subj.proportion, 1) - ...
        models(ind_model).P_target_found_model) ./ ...
        nanstd(humans_data_subj.proportion, 0, 1)) .^2);
    models(ind_model).target_found_distances_by_subj = nanmean(abs(nanmean(humans_data_subj.proportion - ...
        models(ind_model).P_target_found_model , 1)));
    models(ind_model).target_found_distances_std = std(abs(nanmean(humans_data_subj.proportion, 1) - ...
        models(ind_model).P_target_found_model));
    models(ind_model).target_found_distances_by_subj_std = nanmean(abs(nanstd(humans_data_subj.proportion - ...
        models(ind_model).P_target_found_model , 1)));
end

%% Figure 4
% Fig 4A: Proportion of targets found for each model (lines) and humans (boxplot)

ha=[];

figure(2); clf
    set(gcf,'Color','w')
    %set(gcf,'Position',[565 70 430 970])
    %ha(1)=axes('Position',[0.075 0.550 0.85 0.40]); 
%    subplot(3,1,1)
    %    subplot(2,1,1)
    x0=10;
    y0=10;
    width=450;
    height=400;
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
            boxplot(P_target_found,'Notch','on','Color','k','Labels',[2,4,8,12])
            
            for ind_model=1:length(models)
                [correct, M_f, N_f,Pc] = fun_evaluate_experiment(models(ind_model).Nfix_img_model, 57);
                plot(nanmean(Pc'),'.-','Color',models(ind_model).cols,'LineWidth',1.5)
%                 plot(models(ind_model).P_target_found_model,'.-','Color',models(ind_model).cols,'LineWidth',1.5)
            end
        hold off
        names = {models.name};
        [hleg1, icons] = legend(names, 'Location', 'NorthWest', 'FontSize', 11);
        set(gca,'XTickLabels',{'2','4','8','12'})
%             icon = findobj(icons,'Type','line');
%             icon = findobj(icon,'Marker','none','-xor');
%             set(icon,'MarkerSize',30);
%             pause(0.1)
%             textobj = findobj(icons, 'type', 'text');
%             set(textobj, 'fontsize', 12);
        xlabel('Number of max saccades')
        ylabel('Proportion of targets found')
        ylim([0 1])
%             x0=10;
%             y0=10;
%             width=470;
%             height=470;
%             set(gcf,'position',[x0,y0,width,height])


% %Fig 4B: Number of fixations needed to find the target (models and humans)
%     ha(2)=axes('Position',[0.575 0.650 0.300 0.300]); 
%    subplot(3,1,2)
%        hold on
%            x = 1:15;
%            y = hist(Nfix_img_mean,x); y = y/sum(y);
%            plot(x,y,'k-','LineWidth',1.5)
%            for ind_model=1:length(models)
%                [~,~,NFix,~] = fun_evaluate_experiment(models(ind_model).Nfix_img_model, 57);
%                y = hist(nanmean(NFix'),x); y = y/nansum(y); 
%                 %y = hist(models(ind_model).Nfix_img_model,x); y = y/nansum(y); 
%                plot(x,y,'-','Color',models(ind_model).cols,'LineWidth',1.5)
%            end
%        hold off
%        box on
%         %set(gca,'YLim',[0 0.5])
%        set(gca,'YTick',0:0.1:0.4)
%        set(gca,'XLim',[0 13])
%        xlabel('Number of saccades (to find the target)')
%        ylabel('Frequency')
% 
%    subplot(3,1,3)
%    subplot(2,1,2)
%    ha(2)=axes('Position',[0.117 0.075 0.817 0.40]);

%% Save

if 1
%     FigName    = 'correlation-prior';
    FigName    = 'ssim-prior';
%     FigName    = 'ibs-prior';
%     FigName    = 'searchers';
    FolderName = '/home/gastonb/Imágenes/figpaper/fig3/';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    FigHandle = FigList(1);
    saveas(FigHandle, strcat(FolderName, FigName, '.svg'));
    print(gcf,strcat(FolderName, FigName,'.png'),'-dpng','-r400')
end
