% Results exploration
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

% save metrics
guardar = 0;

trials_tmp                    = load(strcat(src_path, 'info_all_subj.mat'));
[ids_images, ~, images_order] = unique({trials_tmp.info_per_subj_final(:).image_name});
images_order                  = images_order';
image_size                    = trials_tmp.info_per_subj_final(1).image_size;
delta                         = 32;
grid_size                     = image_size/delta;
trials                        = reduce_scanpaths(trials_tmp.info_per_subj_final, delta, image_size);
multimatch_names              = {'vectorSim','directionSim','lengthSim','positionSim','durationSim'};
clear trials_tmp

%% Count male/female

folder = '../../all_data';
S = dir(folder); % all of the names are in the structure S anyway.
g=0;
for id=1:length(S)
    if contains(S(id).name,'mat')
        aux = load(strcat(folder,'/',S(id).name)); % but if you really want a cell array of names, do this.
        g = g + str2num(aux.exp.subj.gender);
    end
end

%% Wrong trials

correct_trials_ratio = sum([info_per_subj_final.target_found])/length([info_per_subj_final.target_found]);

correct_trials_subj = nan(Nsubj,1);
correct_trials_img = nan(Nimg,1);

for ind_img=1:Nimg
    img_id = ids_images{ind_img};
    info_per_img = trials(strcmp({trials.image_name}, img_id));
    correct_trials_img(ind_img,1) = sum([info_per_img.target_found]);
end

for ind_subj=1:Nsubj
    subj_id = ids_subjs{ind_subj};
    info_per_subj = trials(strcmp({trials.subj}, subj_id));
    correct_trials_subj(ind_subj,1) = sum([info_per_subj.target_found]);
end

% agregar lindos prints
mean(correct_trials_img), std(correct_trials_img)
mean(correct_trials_subj), std(correct_trials_subj)

%% IoR

% load models
models = fun_define_models('all');
for ind_model = 1:length(models)
    count_revisits = 0;
    nscanpath_revisit = 0; 
    total_fixs = 0;
    for ind_img=1:Nimg
        img_id = ids_images{ind_img};
        info_per_img = trials(strcmp({trials.image_name}, img_id));

        % load variable "scanpath"
        path_scanpath = sprintf('../results_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                        models(ind_model).prior,...
                        models(ind_model).searcher,...
                        models(ind_model).params,...
                        ind_img);
        load(path_scanpath)
        total_fixs = total_fixs + length(scanpath);
        dists = pdist2(scanpath, scanpath);
        zeros_mat = dists == 0; nzeros = sum(sum(zeros_mat));
        nzeros = nzeros -length(scanpath);
        if nzeros > 0
            nscanpath_revisit = nscanpath_revisit + 1;
        end
        count_revisits = count_revisits + nzeros;
        models(ind_model).count_revisits = count_revisits;
        models(ind_model).nscanpath_revisit = nscanpath_revisit;
        models(ind_model).total_fixs = total_fixs;
    end
end

models_table = table({models.name}.', [models.count_revisits].', [models.nscanpath_revisit].', [models.total_fixs].', 'VariableNames', {'name', 'count_revisits', 'nscanpath_revisit', 'total_fixs'});

%%
writetable(models_table, 'ior_analysis.csv')

%% IoR subjects

subj_fixations = nan(1,Nsubj);
subj_count_revisits = nan(1,Nsubj);
subj_nscanpath_revisit = nan(1,Nsubj);

for subj_i=1:Nsubj
    subj_id = ids_subjs{subj_i};
    info_per_subj = trials(strcmp({trials.subj}, subj_id));
                                
    count_revisits = 0;
    nscanpath_revisit = 0; 
    total_fixs = 0;
    for i = 1:length(info_per_subj)
        
        scanpath = [info_per_subj(i).y_grid info_per_subj(i).x_grid];
        total_fixs = total_fixs + length(scanpath);
        dists = pdist2(scanpath, scanpath);
        zeros_mat = dists == 0; nzeros = sum(sum(zeros_mat));
        nzeros = nzeros -length(scanpath);
        if nzeros > 0
            nscanpath_revisit = nscanpath_revisit + 1;
        end
        count_revisits = count_revisits + nzeros;
        subj_count_revisits(subj_i)  = count_revisits;
        subj_nscanpath_revisit(subj_i) = nscanpath_revisit;
        subj_fixations(subj_i) = total_fixs;
    end
end

nanmean(subj_count_revisits)
nanstd(subj_count_revisits)

nanmean(subj_nscanpath_revisit)
nanstd(subj_nscanpath_revisit)

nanmean(subj_fixations)
nanstd(subj_fixations)


%% Plot all target boxes

f = figure(1);
    set(f, 'Color','w')
    f.Position = [50 50 1024 768];
    for i=1:134
        load(sprintf('../results_models/deepgaze/structuralsim/a_3_b_4_tam_celda_32/cfg/cfg_%d.mat',i))
        x_1 = cfg.target_center(2);
        y_1 = cfg.target_center(1);
        x_1 = x_1 -36;
        y_1 = y_1 -36;
        h = rectangle('position',[x_1 y_1 72 72]);
        h.EdgeColor = 'r';
        h.LineWidth = 3;
        
        %daspect([1024 768 1])
        daspect([1.33 1 1])
        hold on
    end
    xlim([0 1024])
    ylim([0 768])
    
    
%% MultiMatch Test

screen_size = [1024, 768];
x_1 = [0, 500, 1023];
y_1 = [0, 350, 0];
t_1 = [200, 200, 200];

x_2 = [0, 500, 0];
y_2 = [0, 350, 767];
t_2 = [200, 200, 200];

scanpath_1 = [x_1', y_1', t_1'];
scanpath_2 = [x_2', y_2', t_2'];
format long
tmp_distance = doComparison(scanpath_1, scanpath_2, screen_size, 0)

x_1 = [0, 16, 31];
y_1 = [0, 12, 0];
t_1 = [200, 200, 200];

x_2 = [0, 16, 0];
y_2 = [0, 12, 23];
t_2 = [200, 200, 200];

scanpath_1 = [x_1', y_1', t_1'];
scanpath_2 = [x_2', y_2', t_2'];
format long
tmp_distance = doComparison(scanpath_1, scanpath_2, screen_size, 0)

x_1 = [0, 16, 31];
y_1 = [0, 12, 0];
t_1 = [200, 200, 200];

x_2 = [0, 16, 0];
y_2 = [0, 12, 23];
t_2 = [200, 200, 200];

scanpath_1 = [x_1', y_1', t_1']+1;
scanpath_2 = [x_2', y_2', t_2']+1;
format long
tmp_distance = doComparison(scanpath_1, scanpath_2, screen_size+2, 0)
