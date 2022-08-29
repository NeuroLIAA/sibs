% Table 1 - Model comparison
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

% plt_type = 'priors';
% plt_type = 'ibs'
plt_type = 'searchers';

%% Load subjects multimatch metrics

addpath('./results_metrics/')
% var_name -> multimatch_bh, mean_dist_img, std_dist_img
subj_data = load('mm_bh_reduced.mat');

%% Subjects Calculations: Average number of fixations needed by the observers to find the target

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

subj_data.fix_found_mean = Nfix_img_mean;
subj_data.fix_found_std  = Nfix_img_std;
subj_data.fix_found_prop = Nfix_img_nsuj;

%% Subjects Calculations: Proportion of targets found as function of the number of saccades allowed for the humans

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

subj_data.target_found           = target_found_img;
subj_data.target_found_ntrials   = target_found_img_ntrials;
subj_data.target_found_mean      = mean_nsac_target_found_img;   
subj_data.proportion             = P_target_found;

%% Load Searchers

% load('mm_hm_reduced_5searchers.mat')
% short_names    = {'cIBS+DGII','sIBS+DGII','IBS+DGII','Greedy+DGII','SaliencyB.'};
load('mm_hm_reduced_all.mat')
% short_names    = {'sIBS+DGII','sIBS+Center','sIBS+Noisy','sIBS+Flat'};

%% Model Calculations: Proportion of targets found as function of the number of saccades allowed for the humans
% Number of fixations predicted by the models
for ind_model = 1:length(models)
    Nfix_img_model = nan(Nimg,1);
    for ind_img=1:Nimg
        path = sprintf('../results_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
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
        path = sprintf('../results_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
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

%% Weighted Distance
for ind_model=1:length(models)

    models(ind_model).target_found_distances_new = mean(...
                ( (nanmean(subj_data.proportion, 1) - models(ind_model).P_target_found_model) ./ ...
                    nanstd(subj_data.proportion, 0, 1)) .^2 ...
                );
    models(ind_model).target_found_distances_new_std = std(...
                ( (nanmean(subj_data.proportion, 1) - models(ind_model).P_target_found_model) ./ ...
                    nanstd(subj_data.proportion, 0, 1)) .^2 ...
                );

    models(ind_model).target_found_distances_by_subj = nanmean(...
                abs(nanmean(subj_data.proportion - models(ind_model).P_target_found_model , 1))...
                );
    
    models(ind_model).target_found_distances_std = std(...
                abs(nanmean(subj_data.proportion, 1) - models(ind_model).P_target_found_model)...
                );
    models(ind_model).target_found_distances_by_subj_std = nanmean(...
                abs(nanstd(subj_data.proportion - models(ind_model).P_target_found_model , 1))...
                );
end

%% Models Mean Agreement and Jaccard Index

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
        % explicit jaccard index
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

%% Human Mean Agreement

% auxiliar calculations
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

% mean agreement
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
subj_data.mean_agreement = mean(subj_agreement(notnan_idx));
subj_data.std_agreement = std(subj_agreement(notnan_idx));
subj_data.subj_agr = subj_agreement;

%% Correlation and Linear Regression

y = nanmedian(subj_data.mean_dist_img, 2);
for i=1:length(models)
    x = nanmedian(models(i).mean_dist_img,2);
    %scatter(x,y,[], models(i).cols, 'filled', 'MarkerFaceAlpha',...
    %    0.6);
    name_aux = {models.name};
    disp(name_aux{i})
    % pearson correlation
    not_nan = isfinite(x) & isfinite(y); 
    [models(i).corr, models(i).pvalue, models(i).iclb, models(i).icub] = ...
        corrcoef(x(not_nan),y(not_nan));
    models(i).corr_value = models(i).corr(1,2);
    
    % linear model
    [models(i).corr_spearman, models(i).pvalue_spearman] = ...
        corr(x(not_nan),y(not_nan),'Type','Spearman');
    mdl = fitlm(x,y,'Intercept',false);
    %models(i).linearregression = mdl;
    models(i).slope       = mdl.Coefficients.Estimate;
    models(i).slope_se    = mdl.Coefficients.SE;
    models(i).slope_tstat = mdl.Coefficients.tStat;
    models(i).slope_pval  = mdl.Coefficients.pValue;
    models(i).lr_rmse     = mdl.RMSE;
end

%% create table

% multimatch_bh = [index, distance];
% multimatch_bh = array2table(multimatch_bh);
% index_names   = {'subj1_id', 'subj2_id', 'img_id'};
% multimatch_bh.Properties.VariableNames = [index_names, multimatch_names];

tabla_metricas_performance = array2table([...
                                            {models.searcher}' {models.prior}' {models.agr}' {models.agrstd}'...
                                            {models.target_found_distances_new}' {models.target_found_distances_new_std}'...
                                            {models.jacc}' {models.jacc_std}' {models.corr_value}'...
                                            ]);
index_names = {'Searcher', 'Prior', 'MAgr', 'MAgrSD', ...
                'WDis', 'WDistSD', 'JaccIn', 'JaccISD', 'Corr'};
tabla_metricas_performance.Properties.VariableNames = index_names;

tabla_mm = array2table(...
                        [nanmean([models.vectorSim])' nanmean([models.directionSim])'...
                         nanmean([models.positionSim])' nanmean([models.lengthSim])']...
                       );

multimatch_names = {'vecSim','dirSim','lenSim','posSim'};
tabla_mm.Properties.VariableNames = multimatch_names;

mm_average = array2table(mean(tabla_mm{:,:},2));
mm_average.Properties.VariableNames = {'MeanMM'};

tabla_resultados = [tabla_metricas_performance mm_average tabla_mm];

%% Optional, basic portage to LaTeX

table2latex(tabla_resultados, './results_metrics/tabla_resultados.tex',2)
