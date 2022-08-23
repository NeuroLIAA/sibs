% GB:
% * Agregar un vector con la similaridad entre el scanpath del modelo y cada humano, deberia tener dim 1x57 o 1x[# de sujetos que encontraron el target]
% * Mandar el script por mail.
% 
% JK:
% * Graficar la grilla (aunque despues la ocultemos).
% * Graficar el punto de inicio sobre la grilla.
% * Pintar todas las celdas que incluyen al target.
% * Colorear los scanpath de acuerdo a la similaridad.
% * Agregar colorbar de referencia.

% It's the same as PaperFigures_Fig1_Behavior_Grid.m but with just the last
% figure, and the option of using the Grid or not.
clear all
close all
clc

% addpath('~/Dropbox/my_functions/')
mainpath = '../../'; %GB
src_path = strcat(mainpath,'data_subjects/data_final/'); 
% mainpath = '../'; %JK
addpath('./utils/')
%addpath([mainpath 'data_analysis/utils/'])
%addpath([mainpath 'data_analysis/utils/heatmap_code/'])
addpath([mainpath 'dynamic_models/utils'])
addpath([mainpath 'results/scripts/new_PaperFigures/utils/'])
%addpath([mainpath 'compare_models/'])
load([mainpath 'matrix/initial_fixations.mat']);
load([src_path 'info_all_subj.mat']);


[subjects, ~, subj_order] = unique({info_per_subj_final(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr            = length(info_per_subj_final);

apply_reduction = 0;

%% (JK 2019-11-14) Junta todos los sujetos. Ya no necesito descartar MVA,
% y recalculo el mapeo a la grilla por las dudas, ya no se guarda.
delta       = 32;
image_size  = [768 1024];
grid_size   = image_size/delta;

a = dir(strcat(src_path,'sinfo_subj/*.mat')); filenames = {a.name};
info_per_subj = [];
for i=1:length(filenames) 
    tmp = load([src_path 'sinfo_subj/' filenames{i}]);
    tmp.info_per_subj = new_subjMapFixationToMatrix( tmp.info_per_subj, path, delta, image_size );
    info_per_subj = [info_per_subj tmp.info_per_subj];
end
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
[images, ~, images_order] = unique({info_per_subj(:).image_name});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);
imgnames = unique({info_per_subj.image_name});

%% Choose Image and models

% Malas imagenes (cIBS+DG):[108;2;63;12;41;107] (133 y 23) - la 2 tiene
% solo un trial
% Buenas imagenes(cIBS+DG): 57, 32, 78, 53, 97, 9

models = fun_define_models('searchers-deepgaze'); 
ind_model = 3; 
tr = 78;
Nsubj_plot = 4; % Amount of scanpaths to plot (x2)

fprintf('\nModel: %s\n',models(ind_model).name);
fprintf('Image: %s\n', imgnames{tr});
info = info_per_subj( strcmp({info_per_subj.image_name},imgnames{tr}) );
path = fullfile(sprintf('../results_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                            models(ind_model).prior,...
                            models(ind_model).searcher,...
                            models(ind_model).params,...
                            tr));
load(path); % scanpath = ( (x,y) , fixation number )

path = fullfile(mainpath,sprintf('out_models/%s/%s/%s/cfg/cfg_%d.mat',...
                            models(ind_model).prior,...
                            models(ind_model).searcher,...
                            models(ind_model).params,...
                            tr));
load(path); % scanpath = ( (x,y) , fixation number )

%% 
% figure(); clf;
%     hold on
%         im = imrea
%         d([mainpath 'templates/' info(1).template]);
%         imagesc(im); colormap gray
%     hold off
%     axis tight
%     set(gca,'YDir','reverse')
%     set(gca,'Visible','off')

figure(tr); clf;
    set(gcf,'Color','w')
    set(gcf,'Position',[64 310 1070 665])
    hold on
        c = summer(100);
        clear hc;
        hc(1) = plot(100,100,'.','MarkerSize',30,'Color',c(1,:)); 
        hc(2) = plot(100,100,'.','MarkerSize',30,'Color',c(26,:)); 
        hc(3) = plot(100,100,'.','MarkerSize',30,'Color',c(51,:)); 
        hc(4) = plot(100,100,'.','MarkerSize',30,'Color',c(76,:)); 
        hc(5) = plot(100,100,'.','MarkerSize',30,'Color',c(100,:)); 

        im = imread([mainpath 'data_images/images/' info(1).image_name]);
        imagesc(im); colormap gray

        ini_pos = initial_fixations( strcmp({initial_fixations.image},imgnames{tr}) ).initial_fix;
        IX = info(1).delta*...
                [floor( ini_pos(2)/info(1).delta ), floor( ini_pos(2)/info(1).delta ), ...
                ceil( ini_pos(2)/info(1).delta ),   ceil( ini_pos(2)/info(1).delta ), ...
                floor( ini_pos(2)/info(1).delta )];
        IY = info(1).delta*...
                [floor( ini_pos(1)/info(1).delta ), ceil( ini_pos(1)/info(1).delta ), ...
                ceil( ini_pos(1)/info(1).delta ),   floor( ini_pos(1)/info(1).delta ), ...
                floor( ini_pos(1)/info(1).delta )];
        %patch(IX,IY,[1 0 0],'EdgeColor','None','FaceAlpha',.5)
        
        %ts      = info(1).exp_trial.tarsize(1)/2;
        ts      =  (info(1).target_rect(3) - info(1).target_rect(1))/2;
        TX = info(1).delta*...
                [floor( (info(1).target_center(1)-ts)/info(1).delta ), ...
                    floor( (info(1).target_center(1)-ts)/info(1).delta ), ...
                    ceil( (info(1).target_center(1)+ts)/info(1).delta ), ...
                    ceil( (info(1).target_center(1)+ts)/info(1).delta ), ...
                    floor( (info(1).target_center(1)-ts)/info(1).delta )];
        TY = info(1).delta*...
                [floor( (info(1).target_center(2)-ts)/info(1).delta ), ...
                    ceil( (info(1).target_center(2)+ts)/info(1).delta ), ...
                    ceil( (info(1).target_center(2)+ts)/info(1).delta ), ...
                    floor( (info(1).target_center(2)-ts)/info(1).delta ), ...
                    floor( (info(1).target_center(2)-ts)/info(1).delta )];
        patch(TX,TY,[0 0 1],'EdgeColor','None','FaceAlpha',.5)
        
        % GB edit - esto aparecia mas adelante
        infoc = info( [info.target_found]==1 & arrayfun(@(x) length(x.y), info)>4);
        grid_size = image_size/infoc(1).delta;
        for i=1:grid_size(1)
            plot([0 image_size(2)],i*infoc(1).delta*[1 1],'w-')
        end
        for i=1:grid_size(2)
            plot(i*infoc(2).delta*[1 1],[0 image_size(1)],'w-')
        end
        
        tx      = info(1).target_center(1) + [-ts -ts ts ts -ts];
        ty      = info(1).target_center(2) + [-ts ts ts -ts -ts];
        for i = 1:4
            plot([tx(i) tx(i+1)],[ty(i) ty(i+1)],'b-','LineWidth',2);
        end
        
        Nsp = length(infoc);
        scanpath_similarity = arrayfun( @(x) scanpathDistance(x.fixations_matrix_reduced, ...
                scanpath, grid_size), infoc );
        ssindex = round(2*scanpath_similarity*100); 
        ssindex(ssindex>100) = 100;
        hc_legend = {'0.000','0.125','0.250','0.375','>0.500'};
%         ssindex = round(scanpath_similarity*100); 
%         hc_legend = {'0.00','0.25','0.50','0.75','1.00'};

% Random
%         for s = randi(Nsp,1,5)%1:Nsp % edit de cantidad de scanpath graficados
%             xh = infoc(s).fixations_matrix_reduced(:,1)*delta + delta/2;
%             yh = infoc(s).fixations_matrix_reduced(:,2)*delta + delta/2;
%             plot(xh,yh,'.-','Color',c(ssindex(s),:),'LineWidth',2,'MarkerSize',20)
%         end

        [ssi j] = sort(ssindex);
% % Best Nsubj_plot
%         for s = j((end-Nsubj_plot):end)
%             xh = (infoc(s).fixations_matrix_reduced(:,1)-1)*delta + delta/2;
%             yh = (infoc(s).fixations_matrix_reduced(:,2)-1)*delta + delta/2;
%             plot(xh,yh,'.-','Color',c(ssindex(s),:),'LineWidth',2,'MarkerSize',20)
%         end
% Worst Nsubj_plot
        for s = j(1:Nsubj_plot)
            xh = (infoc(s).fixations_matrix_reduced(:,1)-1)*delta + delta/2;
            yh = (infoc(s).fixations_matrix_reduced(:,2)-1)*delta + delta/2;
            plot(xh,yh,'.-','Color',c(ssindex(s),:),'LineWidth',2,'MarkerSize',20)
        end
        
        plot((scanpath(:,2)-1)*delta + delta/2, (scanpath(:,1)-1)*delta + delta/2,...
                                            '.-','Color','r','LineWidth',2,'MarkerSize',20)
        
        plot(ini_pos(2),ini_pos(1),'ro','LineWidth',2,'MarkerSize',10);
                
    hold off
    xlabel(strcat('Image: ',num2str(tr), ' - Model: ', models(ind_model).name));
    set(gca,'XLim', [0 image_size(2)],'YLim', [0 image_size(1)])
    axis tight
    %set(gca,'XLabel',fprintf('Image: %s\n', tr))
    set(gca,'YDir','reverse')
    set(gca,'Visible','off')
    h=get(gca);
    h.XAxis.Label.Visible='on';
    legend(hc,hc_legend, ...
                'Location','EastOutside'); legend boxoff
    