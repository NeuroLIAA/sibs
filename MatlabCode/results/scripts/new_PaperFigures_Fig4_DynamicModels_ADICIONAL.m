%% Copiado de Fig4 GeneralComparison v3 - Mas análisis sobre las distancias
%TODO: Pasar a dos figuras nuevas:
%  new_PaperFigure_Fig4_PRIORS_DynamicModels_Correlation
%  new_PaperFigure_Fig4_SEARCHER_DynamicModels_Correlation
clear all
close all
clc
%%
addpath('./utils/')
addpath('../data_analysis/utils/')
addpath('../dynamic_models/utils/')
addpath('../data_analysis/utils/heatmap_code/')
addpath('../compare_models')

src_path = '../new_data/new_matrix/';

load(strcat(src_path,'info_all_subj.mat'));
[subjects, ~, subj_order] = unique({info_per_subj_final(:).subj});
[images_files, ~, images_order] = unique({info_per_subj_final(:).image_name});
images_order = images_order';
subj_order  = subj_order'; % {info_per_subj_final.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj_final);

a=dir('../matrix/images/*.mat');    filenames_img = {a.name}'; Nimg = length(filenames_img);
a=dir('../matrix/subjects/*.mat');  filenames_subj = {a.name};

% Para variar priors case = 6
% Para variar serchers case = 5
models = fun_define_models_tmp(6);
delta = 32;
image_size = info_per_subj_final(1).image_size;
% para no graficar y solo calcular variables
make_plots = true;
%% 
nntrthr     = 20;

% Average number of fixations needed by the observers to find the target
Nfix_img_mean   = nan(Nimg,1);
Nfix_img_std    = nan(Nimg,1);
Nfix_img_nsuj   = nan(Nimg,1);
%data_founds = zeros(Nimg, Nsubj);
for ind_img=1:Nimg % images
    path = char(strcat(src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat'));
    fprintf('%s\n',path);
    load(path); % Loads variable info_per_img
    % Add reduced image data to struct
    info_per_img = new_subjMapFixationToMatrix(info_per_img, '', delta, image_size);
    if (ind_img ~= 132)
        founds = [info_per_img.target_found];
        %data_founds(ind_img,:) = founds;
        fixs = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_img);
        Nfix_img_mean(ind_img)  = mean(fixs(founds));
        Nfix_img_std(ind_img)   = std(fixs(founds));
        Nfix_img_nsuj(ind_img)  = sum(founds);
    end
end
% ?
nntrfilt    = (Nfix_img_nsuj>nntrthr);

%%  
humans_data_img.mean = Nfix_img_mean;
humans_data_img.std = Nfix_img_std;
humans_data_img.prop = Nfix_img_nsuj;

%%  Proportion of targets found as function of the number of saccades allowed for the humans

mean_nsac_target_found_img  = nan(Nsubj,4); % The number of actual saccades on the grid (could be less than the allowed)
target_found_img            = nan(Nsubj,4);
target_found_img_ntrials    = nan(Nsubj,4);
for ind_subj=[1:43 45:Nsubj]
    path = char(strcat(src_path,'sinfo_subj/info_per_subj_', num2str(ind_subj), '.mat')); fprintf('%s\n',path);
    load(path);
    info_per_subj = new_subjMapFixationToMatrix(info_per_subj, '', delta, image_size);
    founds              = [info_per_subj.target_found];
%     fixspp              = arrayfun(@(x) length(x.fixations), info_per_subj);
%     nsaccades           = arrayfun(@(x) x.exp_data.nsaccades, info_per_subj);
    nfixs               = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_subj);
    nsaccades_allowed   = arrayfun(@(x) x.nsaccades_allowed, info_per_subj);
    %if ind_subj==35; keyboard; end
    NSACCADES_allowed   = [2 4 8 12]; %unique(nsaccades_allowed);
    for isacc=1:length(NSACCADES_allowed)
        target_found_img(ind_subj,isacc)            = sum(founds(nsaccades_allowed==NSACCADES_allowed(isacc)));
        target_found_img_ntrials(ind_subj,isacc)    = sum(nsaccades_allowed==NSACCADES_allowed(isacc));
        mean_nsac_target_found_img(ind_subj,isacc)  = mean(nfixs(nsaccades_allowed==NSACCADES_allowed(isacc)));
    end
    %if ind_subj==35; keyboard; end
    if size(target_found_img,2)>4; break; end
end
P_target_found = target_found_img./target_found_img_ntrials;

%% CONTROL

humans_data_subj.target_found           = target_found_img;
humans_data_subj.target_found_ntrials   = target_found_img_ntrials;
humans_data_subj.target_found_mean      = mean_nsac_target_found_img;   
humans_data_subj.proportion             = P_target_found;
    
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
exps_thr = nan(Nimg,Nsubj);
for ind_img=1:Nimg % images
    path = char(strcat(src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat')); 
    fprintf('%s\n',path);
    info_per_img_aux = load(path);

    for  ind_subj = 1:length(info_per_img_aux.info_per_img)
        exps_thr(ind_img,ind_subj)   = info_per_img_aux.info_per_img(ind_subj).nsaccades_allowed;
    end
    %clear info_per_img_aux
end

%% Calculo de targets encontrados MODELOS
% Gaston 7/4 - Creo que esto no se usa en otro lado

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
%% Gaston - distance between human_mean_proportion y models
for ind_model=1:length(models)
    % Pensar: estan dando lo mismo : MIRAR ERROR CUADRATICO
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

%% Gaston - Calculo model vs human in tot (similar al jaccard index)

% para cada modelo
for ind_mod = 1:length(models)
    trials = load(strcat(src_path,'info_all_subj.mat'));
    humans_found = zeros(Nimg, Nsubj);
    models_found = zeros(Nimg, Nsubj);
    humans_and_model_found = zeros(Nimg, Nsubj);
    jaccard_coef = zeros(Nsubj,1);
    for ind_subj=[1:43 45:Nsubj] 
        path = char(strcat(src_path,'sinfo_subj/info_per_subj_', num2str(ind_subj), '.mat'));
        trial = load(path);
        trials = trial.info_per_subj;
        for ind_img = 1:length(trials) % == Nimg
            if (models(ind_mod).Nfix_img_model(ind_img) <= trials(ind_img).nsaccades_allowed+1)
                models_found(ind_img, ind_subj) = 1;
            end
            if trials(ind_img).target_found
                humans_found(ind_img, ind_subj) = 1;
            end
            if humans_found(ind_img, ind_subj) && models_found(ind_img, ind_subj)
                humans_and_model_found(ind_img, ind_subj) = 1;
            end
        end
        jaccard_coef(ind_subj) = ...
            pdist([humans_found(:, ind_subj)'; models_found(:, ind_subj)'], 'jaccard');
%         a = sum(humans_found(:, ind_subj), 1);
%         b = sum(models_found(:, ind_subj), 1);
%         c = sum(humans_and_model_found(:, ind_subj), 1);
%         jaccard_coef(1, ind_subj) = c/(a+b-c);
    end
    % drop subj 44
    jaccard_coef(44) = [];
    models(ind_mod).jacc = mean(1-jaccard_coef);
    models(ind_mod).jacc_std = std(1-jaccard_coef);
    models(ind_mod).agr_subj = mean(abs(models_found - humans_found));
    models(ind_mod).agr = 1 - mean(models(ind_mod).agr_subj);
    models(ind_mod).agrstd = std(models(ind_mod).agr_subj);
end

%% To compare we compute human mean agreement - auxiliar calculations

nsacc = [2,4,8,12];
bool_subj_trials = nan(Nsubj, Nimg, 4);
for ind_trial=1:length(info_per_subj_final)
    ind_subj = subj_order(ind_trial);
    ind_img = images_order(ind_trial);
    ind_nsacc = find(nsacc==info_per_subj_final(ind_trial).nsaccades_allowed);
    if info_per_subj_final(ind_trial).nsaccades_allowed==16 || ...
            info_per_subj_final(ind_trial).nsaccades_allowed==64
        ind_nsacc = 4;
    end
    if info_per_subj_final(ind_trial).nsaccades_allowed==3
        ind_nsacc = 2;
    end
    if length(ind_nsacc) ~= 1
        keyboard
    end
    if info_per_subj_final(ind_trial).target_found
        bool_subj_trials(ind_subj, ind_img, ind_nsacc) = 1;
    else
        bool_subj_trials(ind_subj, ind_img, ind_nsacc) = 0;
    end
end

%% Human mean agreement calculation 
subj_agreement = nan(Nsubj, Nsubj);
for subj_idx = 1:Nsubj
    for subj_idy = subj_idx+1:Nsubj
        subj1 = squeeze(bool_subj_trials(subj_idx,:,:));
        subj2 = squeeze(bool_subj_trials(subj_idy,:,:));
        diff  = abs(subj1-subj2);
        subj_agreement(subj_idx,subj_idy) = 1-nanmean(diff,'all');
    end
end
notnan_idx = isfinite(subj_agreement);
humans.mean_agreement = mean(subj_agreement(notnan_idx));
humans.std_agreement = std(subj_agreement(notnan_idx));
humans.subj_agr = subj_agreement;

%% Scanpath similarity between humans
delta       = 32;
min_fix     = 2;
max_fix     = 13;
image_size  = [768 1024];
grid_size   = image_size/delta; % edit yo le puse delta
    
mean_dist_img   = nan(Nimg,Nsubj);
std_dist_img    = nan(Nimg,Nsubj);
mean_dist_img_nsac = nan(Nimg,Nsubj,4);
eliminados_subj = nan(Nimg,1);
adentro_subj    = nan(Nimg,1);
    
for ind_img=1:Nimg
    if ind_img ~= 132
        path = char(strcat(src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat')); 
        fprintf('%s\n',path);
        load(path);
        
        % (2019-09-18) Map fixations into the Grid
        % Update fixations_matrix_reduced, timefix_matrix_reduced, mascara_matrix_reduced and delta
        % fields in the info_per_subj structure.
        info_per_img = new_subjMapFixationToMatrix( info_per_img, path, delta, image_size );
        % Calculo la distancia de todos contra todos para los sujetos que
        % tengan min_fix < #fix < max_fix y que ambos hayan encontrado el
        % target 
        
        % Calculate the distance between subjects for each image
        ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                            (length(x.fixations_matrix_reduced(:,1)) > min_fix && ...
                                            length(x.fixations_matrix_reduced(:,1)) <= max_fix) ), info_per_img);
        
        distance = nan(Nsubj,Nsubj);
        for subj_i=1:length(info_per_img)
            if ind_subj_target_found(subj_i)
                %for subj_j=(subj_i+1):length(info_per_img) 
                for subj_j=1:length(info_per_img)
                    if subj_j ~= subj_i
                        if ind_subj_target_found(subj_j)
                            tmp_distance = scanpathDistance(info_per_img(subj_i).fixations_matrix_reduced,...
                                                                        info_per_img(subj_j).fixations_matrix_reduced,...
                                                                        grid_size);
                            if ~isnan(tmp_distance)
                                distance(subj_i, subj_j) = tmp_distance;
                            else
                                keyboard
                            end
                        end
                    end
                end
            end
        end
        
        adentro_subj(ind_img)   = sum(ind_subj_target_found);
        % ?
        eliminados_subj(ind_img)= sum(~ind_subj_target_found & [info_per_img.target_found]);
        
        % Calculo la distancia promedio de todos contra todos
        %full_distance = [distance distance'];
        full_distance  = distance;
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
        path = char(strcat(src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat'));
        fprintf('%s\n',path);
        load(path)

        % (2019-09-18) Map fixations into the Grid
        % Update fixations_matrix_reduced, timefix_matrix_reduced, mascara_matrix_reduced and delta
        % fields in the info_per_subj structure.
        info_per_img = new_subjMapFixationToMatrix( info_per_img, path, delta, image_size );
        % Calculo la distancia de todos contra todos para los sujetos que
        % tengan min_fix < #fix < max_fix y que ambos hayan encontrado el
        % target 

        % Calculate the distance between subjects for each image
        ind_subj_target_found = arrayfun(@(x) (x.target_found && ...
                                            (length(x.fixations_matrix_reduced(:,1)) > min_fix && ...
                                            length(x.fixations_matrix_reduced(:,1)) <= max_fix) ), info_per_img);

        for ind_model=1:length(models)
            path = sprintf('../out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                            models(ind_model).prior,...
                            models(ind_model).searcher,...
                            models(ind_model).params,...
                            ind_img);
            load(path)

            if (size(scanpath,1) < max_fix)
                distance = nan(Nsubj,1);
                for subj_i=1:length(info_per_img)
                    if ind_subj_target_found(subj_i)
                        try
                            % EDIT Gaston 7/4 - vi en el chequeo que estan
                            % al reves de como esta en humans despues de la
                            % correccion de JUAN
                            % rechequeado AGOSTO
                            distance(subj_i) = scanpathDistance(info_per_img(subj_i).fixations_matrix_reduced,...
                                                                    [scanpath(:,2) scanpath(:,1)],...
                                                                    grid_size);
                        catch e
                            fprintf(2,"Identificador de error:\n%s \n",e.identifier);
                            fprintf(1,"Mensaje del error:\n%s \n", e.message);
                            keyboard
                            %break
                        end
                    end
                end

                models(ind_model).adentro_subj(ind_img)   = sum(ind_subj_target_found);
                % los eliminados que lo encontraron pero tienen mas de
                % max_fix y menos de min_fix
                models(ind_model).eliminados_subj(ind_img)= sum(~ind_subj_target_found & [info_per_img.target_found]);

                %guardo la distancia por imagen 
                models(ind_model).mean_dist_img(ind_img,:)    = nanmean(distance)';
                models(ind_model).std_dist_img(ind_img,:)     = nanstd(distance)';
            end
        end
    end
end

%% Next: Plots



%% Figure 4
% Fig 4A: Proportion of targets found for each model (lines) and humans (boxplot)

ha=[];
if make_plots
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
            %[hleg1, icons] = legend({'cIBS+DG2','cIBS+Center','cIBS+Flat','cIBS+Noisy'},'Location', 'NorthWest', 'FontSize', 11);
            [hleg1, icons] = legend({'cIBS+DG2','IBS+DG2','Greedy+DG2','SaliencyBased'},'Location', 'NorthWest', 'FontSize', 11);
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
end

% Fig 4B: Number of fixations needed to find the target (models and humans)
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

%    subplot(3,1,3)
%    subplot(2,1,2)
%    ha(2)=axes('Position',[0.117 0.075 0.817 0.40]);

%%
fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 6 3];
%print('SameAxisLimits','-dpng','-r0')
%print(gcf,'Fig4searcher.eps','-deps','-r600')

%%

if make_plots
    figure(3); gcf

        set(gcf,'Color','w')
        Nimg_thr = 25;
            hold on
                filter_img  = adentro_subj>Nimg_thr;
                Nimg_ok     = sum(filter_img);
                %box on;
                plot(zeros(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(mean_dist_img(filter_img),2),'.', 'Color','k','MarkerSize',8)
    %             plot(zeros(Nimg,1) + (2*rand(Nimg,1)-1)*0.1,nanmean(mean_dist_img,2),'k.')
                for ind_model=1:length(models)
                    filter_img  = models(ind_model).adentro_subj>Nimg_thr;
                    Nimg_ok     = sum(filter_img);
                    plot(ind_model*ones(Nimg_ok,1) + (2*rand(Nimg_ok,1)-1)*0.1,nanmean(models(ind_model).mean_dist_img(filter_img),2),'.','Color',models(ind_model).cols,'MarkerSize',8)
    %                 plot(ind_model*ones(Nimg,1) + (2*rand(Nimg,1)-1)*0.1,nanmean(models(ind_model).mean_dist_img,2),'.','Color',models(ind_model).cols)
                end
                plot([-1 length(models)+1],nanmean(nanmean(mean_dist_img,2))*[1 1],'k-')
                plot([-1 length(models)+1],(nanmean(nanmean(mean_dist_img,2)) + nanstd(nanmean(mean_dist_img,2)))*[1 1],'k--')
                plot([-1 length(models)+1],(nanmean(nanmean(mean_dist_img,2)) - nanstd(nanmean(mean_dist_img,2)))*[1 1],'k--')
                for ind_model=1:length(models)
                    plot(ind_model+0.15,nanmean(nanmean(models(ind_model).mean_dist_img,2)),'o','Color',models(ind_model).cols)
                end
                %box off;
            hold off
            set(gca,'YTick',0:0.1:0.4)
            ylabel('Scanpath dissimilarity')
            %set(gca,'XLim',[-1 length(models)+1],'XTick',0:length(models),'XTickLabel',['humans' {models.name}])
            %set(gca,'XTickLabelRotation',-45)
            set(gca,'XLim',[-1 length(models)+1],'XTickLabel',[])
            xlabel('Models')
            % Legend
            [hleg1, icons] = legend({'Humans','cIBS+DG2','cIBS+Center','cIBS+Flat','cIBS+Noisy'},'Location', 'NorthWest', 'FontSize', 13);
            %[hleg1, icons] = legend({'Humans','cIBS+DG2','IBS+DG2','Greedy+DG2','SaliencyBased'},'Location', 'northwest', 'FontSize', 12);
            %lgd = legend({'humans','cIBS+DeepGaze2','IBS+DeepGaze2','Greedy+DeepGaze2','SaliencyBased'},'Location', 'northwest');
            %lgd.FontSize = 11;
            icon = findobj(icons,'Type','line');
            icon = findobj(icon,'Marker','none','-xor');
            set(icon,'MarkerSize',25);
            %pause(0.1)
            %textobj = findobj(icons, 'type', 'text');
            %set(textobj, 'fontsize', 12);
        set(gca,'FontSize',13)
        x0=10;
        y0=10;
        width=400;
        height=550;
        set(gcf,'position',[x0,y0,width,height])
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        pause(0.1)    
        textobj = findobj(icons, 'type', 'text');
        %set(textobj, 'fontsize', 12);
            %legend(['humans' {models.name}],'Location', 'Best');
            %h=get(gca,'XTickLabel');
            %h=get(gca,'xlabel');
            %set(h, 'FontSize', 100)
            %set(gca,'XTickLabelMode','auto')
end

%% FIGURE 4 - Scatter - Mean Scanpath Distance - Models vs Humans

if make_plots
    figure(5); gcf
        set(gcf,'Color','w')   
    %    set(gcf,'Color','w') meterme con hacer la explicación ajustada al esquema, que todavía lo tengo relativamente fresco desde noviembre. Si querés podemos hacer una llamada de 10min mañana para entender bien cómo te gustaría que quede. Podría a la mañana hasta las 12:30hs y después de las 6pm.
    %    set(gcf,'Position',[565 70 430 970])
        y = nanmedian(mean_dist_img, 2);
        x0=10;
        y0=10;
        width= 450;%800; %950
        height= 450;%900; %800
        set(gcf,'position',[x0,y0,width,height])
        %l1 = 0.08;
        %l2 = 0.54;
        %b1 = 0.57;
        %b2 = 0.078;
        %s = 0.42;
        l1 = 0.13; %0.07;
        l2 = 0.57;%l2 = 0.57;
        b1 = 0.56;
        b2 = 0.08;%0.068;
        h = 0.42;
        w = 0.40;
        positions = [l1 b1 w h;...
            l2 b1 w h; ...
            l1 b2 w h; ...
            l2 b2 w h];
        names={'cIBS+DG2','cIBS+Center','cIBS+Flat','cIBS+Noisy'};
        %names={'cIBS+DG2','IBS+DG2','Greedy+DG2','SaliencyBased'};
        for i=1:4
            %grid on
            ha(i) = axes('Position',positions(i,:)); 
            %subplot(2,2,i)
                line('Color','k','LineStyle','--')
                hold on
                x = nanmedian(models(i).mean_dist_img,2);
                h1 = scatter(x,y,[], models(i).cols, 'filled', 'MarkerFaceAlpha',...
                    0.6);
                xlim([0 0.37])
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
                h2 = plot(mdl);
                delete(h2(1));
                h2(2).LineStyle = '-';
                h2(2).LineWidth = 0.8;
                h2(3).LineStyle = '-.';
                h2(3).LineWidth = 0.6;
                h2(4).LineStyle = '-.';
                h2(4).LineWidth = 0.6;
                %legend off;
                legend(h1,names(i),'Location', 'NorthWest','FontSize',12)
                %legend(models(i).name,'FontSize',11)
                %xlabel(models(i).name,'FontSize',11)
                if i==1 || i==2
                    xlabel('')
                    set(gca,'XTick', [])
                end
                if i==3 || i==4
                    xlabel('hmSD')
                end
                if i==1 || i==3
                    ylabel('bhSD','FontSize',11)
                end
                if i==2 || i==4
                    ylabel('')
                    set(gca,'YTick', [])
                end
                title('')
                set(gca,'FontSize',12)
                %refline([1,0],'Color','k')
                mat = [(1:134); x' ; y';abs(x-y)'];
                err = sortrows(mat',4, 'descend');
                % Dropnans
                err(any(isnan(err), 2), :) = [];
                % elijo las primeras tres y las ploteo
                im = err(1:3,1)';
                xt = err(1:3,2)' + 0.003;
                yt = err(1:3,3)' + 0.003;
                str = {arrayfun(@num2str, im', 'UniformOutput', false)};
                str = {strcat('\leftarrow',num2str(im(1))), ...
                        strcat('\leftarrow',num2str(im(2))), ...
                        strcat('\leftarrow',num2str(im(3)))};
                %text(xt,yt,str)
                hold off
        end
        x0=10;
        y0=10;
        width= 650;%800; %950
        height= 550;%900; %800
        set(gcf,'position',[x0,y0,width,height])
end
%% Save Img
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
%print('SameAxisLimits','-dpng','-r0')
%print(gcf,'Fig4searcher_corr_ajuste.png','-dpng','-r600')
    
%% Adicional - save/load

% exp_data

% humans
humans.mean_dist_img    = mean_dist_img;
humans.std_dist_img     = std_dist_img;
humans.subjects_data    = humans_data_subj;
humans.image_data       = humans_data_img;
humans.adentro_subj     = adentro_subj;
humans.eliminados_subj  = eliminados_subj;
save('humans_data.mat','humans')
%save('./Fig4/humans-summary_delta_32.mat', 'humans')

% models
%save('./Fig4/models-searcher_delta_32.mat', 'models')
%model_loaded = load('./Fig4/models-searcher_delta_32.mat');

%% Análisis de los gráficos - NeurIPS RESPUESTA
%  Setear ind_model

ind_model = 1; % deepgaze
% peores imagenes:
%
total_hum_found = sum(humans_found,2);
%hist(total_hum_found)
bin_num = 6;
[N,edges] = histcounts(total_hum_found, bin_num);
images_bin = discretize(total_hum_found,edges);
statarray = grpstats(total_hum_found, images_bin);
X = images_bin; % pertenecia al bin (images_bin) puede ser en porcentual 
Y = nanmedian(models(ind_model).mean_dist_img,2); % similitud a los humanos

figure(6); gcf
    %scatter(X,Y)
    %xlim([0 bin_num+1])
    boxplot(Y,X)
    set(gca,'XTickLabel', round(statarray))
    xlabel('Mean Target Found')
    title('IBS+DG2 similary per image')

%% by image
% Explicacion, fijada una imagen I
% eje Y: 
%       s_j = similitud media del sujeto j al todos los otros sujetos
%       M_j = similitud media del modelo al sujeto j
% eje X imagenes ordenadas
image_subj_subj_sim = nanmean(humans.mean_dist_img, 2); % similarity between subjs
image_model_sim = models(ind_model).mean_dist_img(:,1);
plot_n = 40; % cantidad de imagenes a dibujar
figure(7); gcf
    % distancia a la identid
    distance_ident = abs(image_subj_subj_sim - image_model_sim);
    [values,index]  = sort(distance_ident, 'descend');
    % nan filter
    mask_nan = isnan(values);
    nan_cantidad = sum(mask_nan);
    values_nan = values(~mask_nan);
    index_nan  = index (~mask_nan);
    % groupby
    mean_distance_model = nanmean(abs(image_subj_subj_sim - image_model_sim));
    std_distance_model = nanstd(abs(image_subj_subj_sim - image_model_sim));
    % barplot with error bars of mean similarity
    bar(values_nan(1:plot_n))
    yline(mean_distance_model,'r')
    yline(mean_distance_model+3*std_distance_model,'-.r')
    yline(mean_distance_model-3*std_distance_model,'-.r')
    numero_imagenes = char(index_nan(1:plot_n));
    set(gca, 'Ylim', [0,0.2])
    % agregar nfix y nsubjs
    nsubjs_found = [];
    nfix_mean = [];
    for i=1:plot_n
        image = images_files(numero_imagenes(i));
        s = info_per_subj_final(strcmp({info_per_subj_final.image_name},image));
        nsubjs_found = [nsubjs_found mean([s.target_found])];
        s1 = s([s.target_found]==1);
        nfix = arrayfun(@(x) length(x.y), s1);
        nfix_mean = [nfix_mean mean(nfix)];
    end
    
    set(gca,'XTickLabel', index_nan(1:plot_n))
    xlabel('Image ID')
    ylabel('Similarity Difference')
    title('Mean Similarity Difference')
    
% Chequeo de outliers para varias opciones
out_lim_1  = mean_distance_model + 1 * std_distance_model;
n_out_1    = sum(values_nan >= out_lim_1);
outliers_1 = n_out_1/length(values_nan);

out_lim_2  = mean_distance_model + 2 * std_distance_model;
n_out_2    = sum(values_nan >= out_lim_2);
outliers_2 = n_out_2/length(values_nan);

out_lim_3  = mean_distance_model + 3 * std_distance_model;
n_out_3    = sum(values_nan >= out_lim_3);
outliers_3 = n_out_3/length(values_nan);


%% CUAL ES EL CORRECTO?
control1 = [];
control2 = [];
% tengos dos maneras de medir la cantidad de targets hallados
for ind_img = 1:134
    for ind_subj = 1:57
        if humans_found(ind_img,ind_subj)
            control1 = [control1 humans.mean_dist_img(ind_img,ind_subj)];
        end
        if ~isnan(humans.mean_dist_img(ind_img,ind_subj))
            control2 = [control2 humans_found(ind_img,ind_subj)];
        end    
    end
end

hallado_distancias_medias = 1-isnan(humans.mean_dist_img);
b = sum(hallado_distancias_medias,1);
a = sum(humans_found,1);
c = 0;
l = [];
for i = 1:57
    if a(i)>=b(i)
        c = c + 1; 
    else
        l = [l i];
    end
end

%% Chequeo sobre los found/not_found
hallado_distancias_medias = 1-isnan(humans.mean_dist_img);
hallado_d = sum(hallado_distancias_medias,2);
hallado_f = sum(humans_found,2);
figure(10); gcf
    %map = brewermap(2,'Set1'); %bajar desde git ColorBrewer
    [counts,bins] = histcounts(hallado_d, 10);
    histogram(hallado_d, bins,'FaceColor','blue', 'FaceAlpha', .5, 'EdgeColor', 'none');
    hold on
    histogram(hallado_f, bins,'FaceColor','red', 'FaceAlpha', .5, 'EdgeColor', 'none');
    box off
    legend('Segun nans en dist', 'Segun los datos', 'location', 'northwest')
    legend boxoff
    hold off
    
%% Scanpaths

%plot_scanpaths(12, 1, []);
%% End