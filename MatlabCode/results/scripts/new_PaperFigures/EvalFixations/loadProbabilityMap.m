function [map, fix, tot_subj] = loadProbabilityMap(image, image_ind, nfix, data, model_path)
    % image: filename of the image
    % image_ind: the index of the image to eval
    % nfix: integer of fixation rank (1 < nfix <= nsacc)
    % model: string with the name of the model
    % data: struct (=info_all_subjs) 
    % NOT TAKING INTO ACCOUNT WHETHER THE SUBJECT FOUND OR NOT THE TARGET
    
    if nargin == 4
        % load human fixations
        sprintf('Loading subjects data...')
        data_image = data(strcmp({data.image_name},image));
        %keyboard
        n_subjs = length(data_image);
        fix = nan(n_subjs,2);
        tot_subj = 0;
        for ind_subj = 1:n_subjs
            a = size(data_image(ind_subj).fixations_matrix_reduced);
            if a(1) >= nfix
                fix(ind_subj,:) = [data_image(ind_subj).fixations_matrix_reduced(nfix,2),...
                                        data_image(ind_subj).fixations_matrix_reduced(nfix,1)];  
                tot_subj = tot_subj + 1;
            else
                sprintf(strcat('WARNING: Subj ',num2str(ind_subj), ' missing ', num2str(nfix), 'th fixation.'))
            end
        end
        % create binary map
        image_size_reduced = data(1).image_size/data(1).delta;
        map =  zeros(image_size_reduced(1), image_size_reduced(2));
        for i_s = 1:n_subjs
            if ~isnan(fix(i_s,:))
                map(fix(i_s,1), fix(i_s,2))=1;
            end
        end
    else
        sprintf('Loading model data...')
        % nfix-1 the probability map is indexed according saccades (first probability map correspond to second fixation)
        path = char(strcat(model_path,'probabilityMap_Image_', num2str(image_ind), '_Saccade_', num2str(nfix-1), '.mat')); 
        fix = 0; tot_subj = 0; % not used for model
        map_struct = load(path);
        map = map_struct.probability_map;
    end
end