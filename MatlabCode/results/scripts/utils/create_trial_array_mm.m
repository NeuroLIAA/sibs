function [array1, array2] = create_trial_array_mm(subj1, subj2, data, reduced)
% 30/01/21 Gaston: createTrialArrayMM
%   creates arrays as needed by doComparison in MM toolbox between two
%   subjects ids for a specific image. If image name provided it also
%   selects the specific
% TODO: 
%   - No selecciona imagen
%   - CHEQUEAR XY en version reducida
%   - saco la duración para evitar los tamaños CHEQUEAR
% NOTE: Given current implementation of MultiMatch fixations can not fall
% on screen borders. To fix it, we add to the reduced grid borders of size
% 1.

% TEST:
% data = trial;
% img_name = 'grayscale_11_opensource.jpg';
% subj1_id = 'AFR';
% trial_s1 = data(strcmp({data.subj}, subj1_id));
% trial_s1 = trial_s1(strcmp({trial_s1.image_name}, 'grayscale_11_opensource.jpg'));

    if nargin==3; reduced = true; end

    if reduced
        % select subj1 and slice
        trial_s1 = data(subj1);
        % dummy variable for time
        t1 = 100*ones(length(trial_s1.x_grid),1);
        x1_shif = trial_s1.x_grid +1; y1_shif = trial_s1.y_grid +1;
        array1 = [x1_shif, y1_shif, t1]; 
        
        % select subj2 and slice
        trial_s2 = data(subj2);
        % dummy variable for time
        t2 = 100*ones(length(trial_s2.x_grid),1);
        x2_shif = trial_s2.x_grid +1; y2_shif = trial_s2.y_grid +1;
        array2 = [x2_shif, y2_shif, t2]; 
    else
        trial_s1 = data(subj1);
        array1 = [trial_s1.x', trial_s1.y', 100*ones(size(trial_s1.x))'];

        trial_s2 = data(subj2);
        array2 = [trial_s2.x', trial_s2.y', 100*ones(size(trial_s2.x))']; 
    end
end