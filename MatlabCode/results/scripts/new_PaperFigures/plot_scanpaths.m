function []=plot_scanpaths(image_ind, modelo_ind, info)
% info:
% n_scanpaths (int)
% best (bool)
% plot_model (bool)
% plot_grid  (bool) --> a implementar

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
    %clear all
    %close all
    %clc

    % addpath('~/Dropbox/my_functions/')
    %mainpath = '../../'; %GB
    mainpath = '../'; %JK
    addpath([mainpath 'data_analysis/utils/'])
    addpath([mainpath 'data_analysis/utils/heatmap_code/'])
    addpath([mainpath 'dynamic_models/utils'])
    addpath([mainpath 'PaperFigures/'])
    addpath([mainpath 'compare_models/'])
    load([mainpath 'matrix/initial_fixations.mat']);
    load([mainpath 'matrix/info_per_subj_final.mat']);


    [subjects, ~, subj_order] = unique({info_per_subj(:).subj});
    subj_order  = subj_order'; % {info_per_subj.subj}
    Nsubj       = length(unique(subj_order));
    Ntr         = length(info_per_subj);

    apply_reduction = 0;

    %% (JK 2019-11-14) Junta todos los sujetos. Ya no necesito descartar MVA,
    % y recalculo el mapeo a la grilla por las dudas, ya no se guarda.
    delta       = 32;
    image_size  = [768 1024];
    grid_size   = image_size/delta;
    models_tmp = fun_define_models_tmp(3); 
    %ind_model = 1; 
    ind_model = modelo_ind;
    fprintf('%s\n',models_tmp(ind_model).name);

    a = dir([mainpath 'new_data/new_matrix/sinfo_subj/*.mat']); filenames = {a.name};
    info_per_subj = [];
    for i=1:length(filenames) 
        tmp = load([mainpath 'new_data/new_matrix/sinfo_subj/' filenames{i}]);
        tmp.info_per_subj = new_subjMapFixationToMatrix( tmp.info_per_subj, path, delta, image_size );
        info_per_subj = [info_per_subj tmp.info_per_subj];
    end
    [subjects, ~, subj_order] = unique({info_per_subj(:).subj});
    subj_order  = subj_order'; % {info_per_subj.subj}
    Nsubj       = length(unique(subj_order));
    Ntr         = length(info_per_subj);
    imgnames = unique({info_per_subj.image});

    %%    
    %tr = 108;
    tr = image_ind;

    info = info_per_subj( strcmp({info_per_subj.image},imgnames{tr}) );
    new_path = fullfile(mainpath,sprintf('out_models/%s/%s/%s/scanpath/scanpath_%d.mat',...
                                models_tmp(ind_model).prior,...
                                models_tmp(ind_model).searcher,...
                                models_tmp(ind_model).params,...
                                tr));
    load(new_path); % scanpath = ( (x,y) , fixation number )

    new_path = fullfile(mainpath,sprintf('out_models/%s/%s/%s/cfg/cfg_%d.mat',...
                                models_tmp(ind_model).prior,...
                                models_tmp(ind_model).searcher,...
                                models_tmp(ind_model).params,...
                                tr));
    load(new_path); % scanpath = ( (x,y) , fixation number )

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

            im = imread([mainpath 'images/' info(1).image]);
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
            patch(IX,IY,[1 0 0],'EdgeColor','None','FaceAlpha',.5)

            ts      = info(1).exp_trial.tarsize(1)/2;
            TX = info(1).delta*...
                    [floor( (info(1).target_center(1)-ts)/info(1).delta ), floor( (info(1).target_center(1)-ts)/info(1).delta ), ...
                    ceil( (info(1).target_center(1)+ts)/info(1).delta ), ceil( (info(1).target_center(1)+ts)/info(1).delta ), ...
                    floor( (info(1).target_center(1)-ts)/info(1).delta )];
            TY = info(1).delta*...
                    [floor( (info(1).target_center(2)-ts)/info(1).delta ), ceil( (info(1).target_center(2)+ts)/info(1).delta ), ...
                    ceil( (info(1).target_center(2)+ts)/info(1).delta ), floor( (info(1).target_center(2)-ts)/info(1).delta ), ...
                    floor( (info(1).target_center(2)-ts)/info(1).delta )];
            patch(TX,TY,[0 0 1],'EdgeColor','None','FaceAlpha',.5)

            % GB edit - esto aparecia mas adelante
            infoc = info( [info.target_found]==1 & arrayfun(@(x) size(x.fixations,1), info)>4);
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
    % Best 3
            for s = j((end-3):end)
                xh = (infoc(s).fixations_matrix_reduced(:,1)-1)*delta + delta/2;
                yh = (infoc(s).fixations_matrix_reduced(:,2)-1)*delta + delta/2;
                plot(xh,yh,'.-','Color',c(ssindex(s),:),'LineWidth',2,'MarkerSize',20)
            end
    % Worst 3
            for s = j(1:3)
                xh = (infoc(s).fixations_matrix_reduced(:,1)-1)*delta + delta/2;
                yh = (infoc(s).fixations_matrix_reduced(:,2)-1)*delta + delta/2;
                plot(xh,yh,'.-','Color',c(ssindex(s),:),'LineWidth',2,'MarkerSize',20)
            end

            plot((scanpath(:,2)-1)*delta + delta/2, (scanpath(:,1)-1)*delta + delta/2,...
                                                '.-','Color','r','LineWidth',2,'MarkerSize',20)

            plot(ini_pos(2),ini_pos(1),'ro','LineWidth',2,'MarkerSize',10);

        hold off
        set(gca,'XLim', [0 image_size(2)],'YLim', [0 image_size(1)])
        axis tight
        set(gca,'YDir','reverse')
        set(gca,'Visible','off')
        legend(hc,hc_legend, ...
                    'Location','EastOutside'); legend boxoff
end