function [array1, array2] = createTrialArrayMM(subj1, subj2, image,data,reduced)
% 30/01/21 Gaston: createTrialArrayMM
%   creates arrays as needed by doComparison in MM toolbox between two
%   subjects ids for a specific image. If image name provided it also
%   selects the specific
% TODO: 
%   - No selecciona imagen
%   - CHEQUEAR XY en version reducida
%   - saco la duración para evitar los tamaños CHEQUEAR

% prueba
% data = trail;
% img_name = 'grayscale_11_opensource.jpg';
% subj1_id = 'AFR';
% trial_s1 = data(strcmp({data.subj}, subj1_id));
% trial_s1 = trial_s1(strcmp({trial_s1.image_name}, 'grayscale_11_opensource.jpg'));

    if nargin==4
        reduced = true;
    end

    if reduced
        % select subj1 and slice
        %trial_s1 = data(strcmp({data.subj}, subj1_id));
        trial_s1 = data(subj1);
        array1 = [trial_s1.fixations_matrix_reduced(:,2), trial_s1.fixations_matrix_reduced(:,1), 100*ones(size(trial_s1.fixations_matrix_reduced(:,1)))]; % trial_s1.dur'
        %array1 = [trial_s1.x',trial_s1.y', trial_s1.dur'];

        % select subj2 and slice
        %trial_s2 = data(strcmp({data.subj}, subj2_id));
        trial_s2 = data(subj2);
        array2 = [trial_s2.fixations_matrix_reduced(:,2), trial_s2.fixations_matrix_reduced(:,1), 100*ones(size(trial_s2.fixations_matrix_reduced(:,1)))]; %trial_s2.dur'
        %array2 = [trial_s2.x',trial_s2.y', trial_s2.dur'];
    else
        trial_s1 = data(subj1);
        array1 = [trial_s1.x',trial_s1.y', 100*ones(size(trial_s1.x))'];

        trial_s2 = data(subj2);
        array2 = [trial_s2.x',trial_s2.y', 100*ones(size(trial_s2.x))']; 
    end
end