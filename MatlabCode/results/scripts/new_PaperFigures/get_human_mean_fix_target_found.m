function [mean_nsac_target_found_img, target_found_img, target_found_img_ntrials, ...
                            P_target_found] = get_human_mean_fix_target_found(exp, data)

% Proportion of targets found as function of the number of saccades allowed for the humans
% INPUT:
%     - exp: struct that contains all experiment's variables settings
%     - data: struct that contains all subjects data
% OUTPUT:
%     - mean_nsac_target_found_img
%     - target_found_img
%     - target_found_img_ntrials
%     - P_target_found

    mean_nsac_target_found_img  = nan(exp.Nsubj,4); % The number of actual saccades on the grid (could be less than the allowed)
    target_found_img            = nan(exp.Nsubj,4);
    target_found_img_ntrials    = nan(exp.Nsubj,4);
    for ind_subj=[1:43 45:exp.Nsubj]
        %path = char(strcat(src_path,'sinfo_subj/info_per_subj_', num2str(ind_subj), '.mat')); fprintf('%s\n',path);
        %load(path);
        info_per_subj_ind = data([exp.subj_order]==ind_subj);
        info_per_subj     = new_subjMapFixationToMatrix(info_per_subj_ind, '', exp.delta, exp.image_size);
        founds            = [info_per_subj.target_found];
        %fixspp            = arrayfun(@(x) length(x.fixations), info_per_subj);
        %nsaccades         = arrayfun(@(x) x.exp_data.nsaccades, info_per_subj);
        nfixs             = arrayfun(@(x) length(x.fixations_matrix_reduced), info_per_subj);
        nsaccades_allowed = arrayfun(@(x) x.nsaccades_allowed, info_per_subj);
        for isacc=1:length(exp.nsacc)
            target_found_img(ind_subj,isacc)            = sum(founds(nsaccades_allowed==exp.nsacc(isacc)));
            target_found_img_ntrials(ind_subj,isacc)    = sum(nsaccades_allowed==exp.nsacc(isacc));
            mean_nsac_target_found_img(ind_subj,isacc)  = mean(nfixs(nsaccades_allowed==exp.nsacc(isacc)));
        end
        %if ind_subj==35; keyboard; end
        if size(target_found_img,2) > 4; keyboard; end
    end
    P_target_found = target_found_img./target_found_img_ntrials;

end