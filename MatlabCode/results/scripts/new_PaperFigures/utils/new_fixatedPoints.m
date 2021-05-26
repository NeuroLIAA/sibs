function new_fixatedPoints(radii, minFix, maxFix, src_path, out_path)
% FIXATED POINTS: Creates binary matrix the same size as the original image,
% putting 1 on points that were fixated by at least one subject, 0 otherwise
    % It has the capability of computing 2 different matrices, separating
    % subjects in train and test set (params trainSubjects, testSubjects)
    % For all subjects, put testSubjects = [] 

    % radii = radii around pixel that is considered as fixated (default: 3)
    % [minFix, maxFix] = interval of fixations numbers that are considered
    if nargin == 3
        src_path = '../../new_data/new_matrix/';
        out_path = '../../new_data/new_matrix/';
        old_src_path = '../../matrix/';
        old_file_name = 'info_per_subj_final.mat';
    end
    
    % Maybe I could avoid it if the file exist: if ~isfile()
    num_images = 134;           % can be set from data
    load(strcat(src_path,'info_all_subj.mat')); % it has all subjects
    [subjects, ~, ~] = unique({info_per_subj_final(:).subj});
    dims = info_per_subj_final(:).image_size;   % GASTON - [alto ancho]
    x_lim = dims(2);                            % 1024
    y_lim = dims(1);                            % 768

    %testSubjects = randsample(1:length(subjects), length(subjects) / 4);
    testSubjects = [];
    trainSubjects = setdiff(1:length(subjects), testSubjects);

    fixatedArea      = zeros(dims(1), dims(2), num_images);
    fixatedAreaTrain = zeros(dims(1), dims(2), num_images);
    fixatedAreaTest  = zeros(dims(1), dims(2), num_images);
    for im = 1:num_images
        im
        load(strcat(src_path, 'sinfo_img/info_per_img_', int2str(im), '.mat'));
        filter = info_per_img;
        for su = 1:length(filter)
            fixations = [filter(su).x' filter(su).y'];
            %keyboard
            if size(fixations, 1) < minFix
                continue;
            end
            for f = minFix:min(maxFix, size(fixations, 1))
                fix_x = round(fixations(f,1));
                %keyboard
                fix_y = round(fixations(f,2));
                % remove fixation from non fixated salient
                % Best resolution of fixation: (tan(0.5 degrees) * 60cm * 2) * 1280px / 35cm = 38.3px
                for x = (fix_x-radii):(fix_x+radii)
                    for y = (fix_y-radii):(fix_y+radii)
                        if norm([x y] - [fix_x fix_y]) <= radii && ...
                                0 < x && x <= x_lim && 0 < y && y <= y_lim
                            if ismember(filter(su).subj, subjects(testSubjects))
                                % save matrix
                                fixatedAreaTest(y, x, im) = 1;
                            elseif ismember(filter(su).subj, subjects(trainSubjects))
                                fixatedAreaTrain(y, x, im) = 1;
                            end
                        end
                    end
                end
           end
        end
    end

    if isempty(testSubjects)
        fixatedArea = fixatedAreaTrain;
    end

    if (radii==3)
        save(strcat(out_path, sprintf('fixated_area_fix_%d_to_%d.mat', minFix, maxFix)), ...
           'fixatedArea', 'fixatedAreaTrain', 'fixatedAreaTest', 'trainSubjects', 'testSubjects', 'radii');
    else
        save(strcat(out_path, sprintf('fixated_area_fix_%d_to_%d_radii_%d.mat', minFix, maxFix,radii)), ...
           'fixatedArea', 'fixatedAreaTrain', 'fixatedAreaTest', 'trainSubjects', 'testSubjects', 'radii');
    end
end