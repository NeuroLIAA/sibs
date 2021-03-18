load('info_all_subj.mat')
fid = fopen('info_per_subj_final.json','wt');
fprintf(fid, jsonencode(info_per_subj_final););
fclose(fid);
