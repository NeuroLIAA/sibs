% It's the same as PaperFigures_Fig1_Behavior_Grid.m but with just the last
% figure, and the option of using the Grid or not.
clear all
close all
clc

addpath('~/Dropbox/my_functions/')
addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
addpath('../dynamic_models/utils')

load('../matrix/initial_fixations.mat');

load('../matrix/info_per_subj_final.mat');
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

a = dir('../matrix/subjects/*.mat'); filenames = {a.name};
info_per_subj = [];
for i=1:length(filenames)
%     tmp = load(['../matrix/subjects/' filenames{i}]);
%     if ~strcmp(tmp.info_per_subj(1).subj,'MVA')
%         info_per_subj = [info_per_subj tmp.info_per_subj];
%     end
    tmp = load(['../matrix/subjects/' filenames{i}]);
    tmp.info_per_subj = subjMapFixationToMatrix( tmp.info_per_subj, path, delta, image_size );
    info_per_subj = [info_per_subj tmp.info_per_subj];
end
[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);

%% Figure 1. Experimental design and general behavior. 
% Nota (19/06/2018):
%   Se agregaron nuevos campos en los que se colapsan las fijaciones
%   sucesivas dentro de una posicion en la misma grilla que se usa para
%   correr el modelo de Geisler. 
%   Se recalculan los tamaños y direccion de las sacadas.
%   info_per_subj.fixations -> info_per_subj.fixations_matrix_reduced
%   info_per_subj.timeFix -> info_per_subj.timefix_matrix_reduced
%   NO ESTA DEFINIDO EN LA ESTRUCTURA info_per_subj.saccades -> info_per_subj.saccades_matrix_reduced
%   info_per_subj.mascara_matrix_reduced: que fijaciones colapsar

%% A) Timeline of a trial. The search could be stopped by fixating into the 
% target region of after N saccades, where N varied across trials. 

%% B) Proportion of targets found as function of the number of saccades allowed. 
% NOTA: NO USO LA POSICION ACA.
nsacc_allowed   = arrayfun(@(x) x.exp_trial.nsaccades_thr,info_per_subj);
target_found    = [info_per_subj.target_found];

% NSACC_ALLOWED = unique(nsacc_allowed);
NSACC_ALLOWED = [2 4 8 12];
n   = nan(Nsubj, length(NSACC_ALLOWED));
p   = nan(Nsubj, length(NSACC_ALLOWED));
for su = 1:Nsubj
    for nn = 1:length(NSACC_ALLOWED)
        ind = (subj_order==su & nsacc_allowed==NSACC_ALLOWED(nn));
            if any(ind)
                n(su,nn) = sum(ind);
                p(su,nn) = mean(target_found(ind));
            end                
    end
end
    
%% C) Distributions of saccade length for the first (red), the fifth 
% (green) and the tenth (blue) saccade. 
% NOTA: NO USO LA POSICION ACA, SOLO LA CANTIDAD DE SACADAS.
if apply_reduction
    nsacc_used      = arrayfun(@(x) size(x.fixations_matrix_reduced,1),info_per_subj);
else
    nsacc_used      = arrayfun(@(x) size(x.fixations,1),info_per_subj);
end
target_found    = [info_per_subj.target_found]';
[images,~,trialimage_id] = unique({info_per_subj.image});

nsacc = nan(max(trialimage_id),1);
for i = 1:max(trialimage_id)
    nsacc(i) = nanmedian(nsacc_used(target_found==1 & trialimage_id==i));
end
nsacc = nsacc(~isnan(nsacc));

%% D) Mean (± s.e.m.) saccade length as function of saccade rank. 
% For this we're using the fixations instead the saccades itself. Check the
% pipeline to generate the saccade structure.
fixations   = [];
saccades    = [];
for tr = 1:Ntr
    if apply_reduction        
        nfix = size(info_per_subj(tr).fixations_matrix_reduced,1);
    else
        nfix = size(info_per_subj(tr).fixations,1);
    end
    if (nfix > 1)
        if apply_reduction        
            fixations       = [fixations; ...
                [repmat([tr subj_order(tr) nsacc_allowed(tr)],nfix,1), [1:nfix]',  info_per_subj(tr).fixations_matrix_reduced]];
            dx = info_per_subj(tr).fixations_matrix_reduced(2:end,1) - info_per_subj(tr).fixations_matrix_reduced(1:end-1,1);
            dy = info_per_subj(tr).fixations_matrix_reduced(2:end,2) - info_per_subj(tr).fixations_matrix_reduced(1:end-1,2);
%             dx = info_per_subj(tr).fixations_matrix_reduced(2:end,2) - info_per_subj(tr).fixations_matrix_reduced(1:end-1,2);
%             dy = info_per_subj(tr).fixations_matrix_reduced(2:end,1) - info_per_subj(tr).fixations_matrix_reduced(1:end-1,1);
        else
            fixations       = [fixations; ...
                [repmat([tr subj_order(tr) nsacc_allowed(tr)],nfix,1), [1:nfix]',  info_per_subj(tr).fixations]];
            dx = info_per_subj(tr).fixations(2:end,1) - info_per_subj(tr).fixations(1:end-1,1);
            dy = info_per_subj(tr).fixations(2:end,2) - info_per_subj(tr).fixations(1:end-1,2);
%             dx = info_per_subj(tr).fixations(2:end,2) - info_per_subj(tr).fixations(1:end-1,2);
%             dy = info_per_subj(tr).fixations(2:end,1) - info_per_subj(tr).fixations(1:end-1,1);
            
        end
        
        theta = atan(dy./dx); % atan:[-Inf,Inf]->[-pi/2,pi/2], y lo quiero en [0,2*pi] 
            theta(dx>=0 & dy>=0)= theta(dx>=0   & dy>=0);
            theta(dx<0  & dy>=0)= theta(dx<0    & dy>=0)+ pi;
            theta(dx<0  & dy<0) = theta(dx<0    & dy<0) + pi;
            theta(dx>=0 & dy<0) = theta(dx>=0   & dy<0) + 2*pi;

        saccamp = sqrt(dx.^2 + dy.^2);

        saccades        = [saccades; ...
                [repmat([tr subj_order(tr) nsacc_allowed(tr)],nfix-1,1), [1:nfix-1]', saccamp, theta]];
    end
end

% SACCRANK    = 1:12;
% s = nan(Nsubj,length(SACCRANK));
% for su = 1:Nsubj
%     for i=1:length(SACCRANK)
%         s(su,i) = nanmedian(saccades(saccades(:,2)==su & saccades(:,4)==SACCRANK(i),5));
%     end
% end

%% E) Distributions of saccade direction (measured in degrees from the 
% positive horizontal axis) for the first (red), the fifth (green) and the 
% tenth (blue) saccade. 
    
%% F) Spatial distributions of the first, the fifth and the tenth fixation 
% position. Warmed colors indicate larger probabilities.
settings.filterSize = 300; % tama?o del filtro gaussiano en px.
settings.filterSD   = 30;  % SD del filtro gaussiano en px.
settings.max_nhist  = 50;  % controla el colormap. ajustarlo para que no sature el heatmap.
settings.alpha      = 0.8; % entre 0 y 1. cuanto mayor, mas se ven las fijaciones.
settings.grayscale  = 1;
im_file             = '../data_analysis/utils/white.jpg';
im                  = imread(im_file);  
[s1,s2]             = size(im);

FIXRANK    = 1:12;
M   = zeros(length(FIXRANK),s1,s2);

if apply_reduction
    sg1 = max(fixations(:,5));
    sg2 = max(fixations(:,6));
%     sg1 = max(fixations(:,6));
%     sg2 = max(fixations(:,5));
    MG  = zeros(length(FIXRANK),sg1,sg2);
end
for i=1:length(FIXRANK)
    x = fixations(fixations(:,4)==FIXRANK(i),5);
    y = fixations(fixations(:,4)==FIXRANK(i),6);
%     x = fixations(fixations(:,4)==FIXRANK(i),6);
%     y = fixations(fixations(:,4)==FIXRANK(i),5);
    m = fun_heatmap(im_file, x, y, settings);
    M(i,:,:) = m/max(m(:));
    
    if apply_reduction
        mg = zeros(sg1,sg2);
        for ii=1:sg1;
            for jj=1:sg2
                mg(ii,jj) = sum(x==jj & y==ii);
            end
        end
        MG(i,:,:) = mg/max(mg(:));
    end
end
       
%% All-in-one fig...
SACCRANK    = 1:12;
figure(100);clf
    set(gcf,'Color','w')
    set(gcf,'Position',[66 1 1855 1001])
    
    tr = 1;
    axes('Position',[0.050, 0.850, 0.1, 0.1])
        hold on
            im = imread(['../templates/' info_per_subj(1).template]);
            [z1,z2] = size(im);
            tmp     = 128*ones(s1,s2);
            ind1    = (-z1/2:(z1/2-1)) + s1/2; 
            ind2    = (-z2/2:(z2/2-1)) + s2/2; 
            tmp(ind1,ind2) = im;
            
            imagesc(tmp); colormap gray
            plot([1 s2],[1 1],  '-','Color','k')
            plot([1 s2],[s1 s1],'-','Color','k')
            plot([1 1], [1 s1], '-','Color','k')
            plot([s2 s2],[1 s1],'-','Color','k')
        hold off
        axis tight
        set(gca,'YDir','reverse')
        set(gca,'Visible','off')
    axes('Position',[0.110, 0.775, 0.1, 0.1])
        hold on
            tmp     = 128*ones(s1,s2);
            imagesc(tmp); colormap gray
            plot(info_per_subj(tr).exp_trial.pos(2),...
                    info_per_subj(tr).exp_trial.pos(1), ...
                                '.','Color','k','MarkerSize',15)
            plot([1 s2],[1 1],  '-','Color','k')
            plot([1 s2],[s1 s1],'-','Color','k')
            plot([1 1], [1 s1], '-','Color','k')
            plot([s2 s2],[1 s1],'-','Color','k')
        hold off
        axis tight
        set(gca,'Visible','off')
    axes('Position',[0.170, 0.700, 0.1, 0.1])
        hold on
            tr = 1;
            im = imread(['../images/' info_per_subj(tr).image]);
            imagesc(im); colormap gray
            plot([1 s2],[1 1],  '-','Color','k')
            plot([1 s2],[s1 s1],'-','Color','k')
            plot([1 1], [1 s1], '-','Color','k')
            plot([s2 s2],[1 s1],'-','Color','k')
        hold off
        axis tight
        set(gca,'YDir','reverse')
        set(gca,'Visible','off')
    axes('Position',[0.230, 0.625, 0.1, 0.1])
        hold on
            tmp     = 128*ones(s1,s2);
            imagesc(tmp); colormap gray
            plot(info_per_subj(tr).exp_trial.info.rect_target(1),...
                    info_per_subj(tr).exp_trial.info.rect_target(2), ...
                                '.','Color','k','MarkerSize',30)
            plot([1 s2],[1 1],  '-','Color','k')
            plot([1 s2],[s1 s1],'-','Color','k')
            plot([1 1], [1 s1], '-','Color','k')
            plot([s2 s2],[1 s1],'-','Color','k')
        hold off
        axis tight
        set(gca,'Visible','off')
    
    
    axes('Position',[0.40, 0.625, 0.25, 0.325])
        boxplot(p,'Notch','on','Color','k','Labels',NSACC_ALLOWED)
        box on
        set(gca,'XTick',NSACC_ALLOWED)
        xlabel('Number of saccades allowed')
        ylabel('P(target found)')
    
    axes('Position',[0.70, 0.625, 0.25, 0.325])
        hold on
            xh = 1:12;
            [yh]=hist(nsacc,xh); yh=yh/sum(yh);
            h = bar(xh,yh);
        hold on
        box on
        set(h,'FaceColor',[.7 .7 .7])
        xlabel('Number of saccades to find the target')
        ylabel('Prob.')

    if apply_reduction
        xh = 0:2:25;
    else
        xh = 0:20:1000;
    end
    cols        = rainbow_colors(length(SACCRANK));
    for i=1:length(SACCRANK)
        axes('Position',[0.05+0.075*(i-1), 0.40, 0.0725, 0.15])
            yh = hist(saccades(saccades(:,4)==SACCRANK(i),5),xh); yh = yh/sum(yh);
            plot(xh,yh,'Color',cols(i,:),'LineWidth',1.5)
            if apply_reduction            
                set(gca,'YLim', [0 0.305],'XLim',[0 25])
                set(gca,'XTick',[10 20])
            else
                set(gca,'YLim', [0 0.13],'XLim',[0 600])
                set(gca,'XTick',[0 500])
            end
            box off
            grid on
            if (i==1)
                set(gca,'YTick',[0 0.1 0.2 0.3])
%                 xlabel('Saccade amplitude (pxs)')
                xlabel('Saccade amplitude (grid pos)')
                ylabel('Prob.')
            else
                set(gca,'YTickLabel',[],'XTickLabel',[])
            end
    end
    cols    = rainbow_colors(length(SACCRANK));
    for i = 1:length(SACCRANK)
        axes('Position',[0.05+0.075*(i-1), 0.20, 0.0725, 0.15])
            [T,R]=rose(saccades(saccades(:,4)==SACCRANK(i),6)); R = R/length(saccades(saccades(:,4)==SACCRANK(i),6));
            polarplot(T,R,'Color',cols(i,:),'LineWidth',1.5)
            set(gca,'Visible','off')
    end
        axes('Position',[0.05, 0.20, 0.0725, 0.15])
            hold on
                text(s2/2,s1*0.9,'Saccade Direction','HorizontalAlignment','center','FontSize',12); 
                plot([1 s2],[1 1],  '-','Color','w')
                plot([1 s2],[s1 s1],'-','Color','w')
                plot([1 1], [1 s1], '-','Color','w')
                plot([s2 s2],[1 s1],'-','Color','w')
            hold off
            axis tight
            set(gca,'Visible','off')
    cols = rainbow_colors(length(FIXRANK));
    for i = 1:length(FIXRANK)
        axes('Position',[0.05+0.075*(i-1), 0.05, 0.0725, 0.15])
            hold on
                if (i==1); 
                    text(s2/2,s1*1.1,'Fixation position','HorizontalAlignment','center','FontSize',12); 
                end
                if apply_reduction                
                    pcolor(1-squeeze(MG(i,:,:))); colormap gray
                    plot([1 sg2],[1 1],  '-','Color',cols(i,:),'LineWidth',3)
                    plot([1 sg2],[sg1 sg1],'-','Color',cols(i,:),'LineWidth',3)
                    plot([1 1], [1 sg1], '-','Color',cols(i,:),'LineWidth',3)
                    plot([sg2 sg2],[1 sg1],'-','Color',cols(i,:),'LineWidth',3)
                else
                    imagesc(1-squeeze(M(i,:,:))); colormap gray
                    plot([1 s2],[1 1],  '-','Color',cols(i,:),'LineWidth',3)
                    plot([1 s2],[s1 s1],'-','Color',cols(i,:),'LineWidth',3)
                    plot([1 1], [1 s1], '-','Color',cols(i,:),'LineWidth',3)
                    plot([s2 s2],[1 s1],'-','Color',cols(i,:),'LineWidth',3)
                end
            hold off
            axis tight
            if apply_reduction                
                set(gca,'YLim',[0.5 sg1+0.5],'XLim',[0.5 sg2+0.5])
            else
                set(gca,'YLim',[0.5 s1+0.5],'XLim',[0.5 s2+0.5])
            end
            set(gca,'YDir','reverse')
            set(gca,'Visible','off')
    end

