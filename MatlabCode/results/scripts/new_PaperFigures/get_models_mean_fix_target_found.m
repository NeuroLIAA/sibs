function models = get_models_mean_fix_target_found(exp, models)

% Number of fixations predicted by the models
% INPUT:
%     - exp: struct that contains all experiment's variables settings
%     - models: struct that contains all subjects data
% OUTPUT: (for each model)
%     - target_found_img_model
%     - P_target_found_model
%     - P_target_found_model_std
    
    for ind_model = 1:length(models)
        Nfix_img_model = nan(exp.Nimg,1);
        for ind_img=1:exp.Nimg
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
    exps_thr = nan(exp.Nimg, exp.Nsubj);
    for ind_img=1:exp.Nimg % images
        path = char(strcat(exp.src_path,'sinfo_img/info_per_img_', num2str(ind_img), '.mat')); 
        fprintf('%s\n',path);
        info_per_img_aux = load(path);

        for  ind_subj = 1:length(info_per_img_aux.info_per_img)
            exps_thr(ind_img,ind_subj)   = info_per_img_aux.info_per_img(ind_subj).nsaccades_allowed;
        end
        %clear info_per_img_aux
    end

    for ind_model=1:length(models)
        target_found_img_model     = zeros(exp.Nimg,4);
    %     target_notfound_img_model  = zeros(Nimg,4);
        for ind_img=1:exp.Nimg
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
end
    