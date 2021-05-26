%% NOTA IMPORTANTE: TIRAR SUJETO 44!!! TIENE 0 FIJACIONES EN TODOS LOS TRIALS
ind_subj=44;
path = char(strcat('../matrix/subjects/info_per_subj_subj_', num2str(ind_subj), '.mat')); fprintf('%s\n',path);
load(path);
arrayfun(@(x) length(x.fixations), info_per_subj)

%%
for ind_img=1:Nimg % images
    path = char(strcat('../matrix/images/info_per_subj_img_', num2str(ind_img), '.mat')); fprintf('%s\n',path);
    load(path);
    if (ind_img ~= 132)
        ind = find(ismember({info_per_subj.subj},'MVA'));
        size(info_per_subj(ind).fixations,1)
    end
end