function fun_genero_new_model_vs_subj_file_dg(dynamic_model,ab_value,delta,min_fix,max_fix,image_size)
    % % Se le aplica la grilla a los sujetos. Si ya fue aplicada, comentar.
    % for img=1:134
    %     if img ~= 132
    %         path = char(strcat('../../matrix/images/info_per_subj_img_', num2str(img), '.mat')); fprintf('%s\n',path);
    %         load(path)
    %         info_per_subj = subjMapFixationToMatrix( info_per_subj, path, delta, image_size );
    %     end
    % end

    a = ab_value(1);
    b = ab_value(2);

    mean_dist_model_corr = [];
    std_dist_model_corr = [];

    eliminados_model = [];
    adentro_model = [];
    mean_length_subj_corr = [];
    length_model = [];
    std_length_subj_corr = [];
    for img=1:134
        if img ~= 132
            aden = 0;
            elim = 0;

            path = char(strcat('../../out_models/deepgaze/',dynamic_model,...
                '/others/a_',num2str(a), '_b_', num2str(b), '_tam_celda_',...
                num2str(delta),'/scanpath/scanpath_', num2str(img), '.mat')); fprintf('%s\n',path);
            load(path)
            
            path = char(strcat('../../new_data/new_matrix/sinfo_img/info_per_img_', num2str(img), '.mat')); fprintf('%s\n',path);
            load(path)
            
            % Calculo fijaciones reducidas
            info_per_img = new_subjMapFixationToMatrix(info_per_img, '', delta, image_size);

            % Calculo la distancia de todos contra todos
            model_distance = [];
            length_subj_img = [];

            if length(scanpath(:,1)) <= max_fix && length(scanpath(:,1)) >= min_fix
                for subj_i=1:length(info_per_img)
                    if info_per_img(subj_i).target_found && ...
                            length(info_per_img(subj_i).fixations_matrix_reduced(:,1)) >= min_fix && ...
                            length(info_per_img(subj_i).fixations_matrix_reduced(:,1)) <= max_fix 
                        model_distance = [model_distance scanpathDistance(info_per_img(subj_i).fixations_matrix_reduced,scanpath)];
                        length_subj_img = [length_subj_img length(info_per_img(subj_i).fixations_matrix_reduced(:,1))];
                    end
                end
            end

            if length(scanpath(:,1)) <= max_fix
                aden = aden + 1;
            else
                elim = elim + 1;
            end

            adentro_model = [adentro_model aden];
            eliminados_model = [eliminados_model elim];

            mean_dist_model_corr = [mean_dist_model_corr mean(model_distance)];
            std_dist_model_corr = [std_dist_model_corr std(model_distance)];

            mean_length_subj_corr = [mean_length_subj_corr mean(length_subj_img)];
            std_length_subj_corr = [std_length_subj_corr std(length_subj_img)];
            length_model = [length_model length(scanpath(:,1))];

        end
    end
    
    save_name = char(strcat('model_vs_subj_',dynamic_model,'_a_',num2str(a), '_b_', num2str(b), '_delta_',num2str(delta),'.mat'));
    save(save_name, 'adentro_model', 'eliminados_model', 'mean_dist_model_corr', 'std_dist_model_corr', 'mean_length_subj_corr', 'length_model', 'std_length_subj_corr')
end