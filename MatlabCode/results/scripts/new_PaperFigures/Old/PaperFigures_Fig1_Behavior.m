clear all
close all
clc

addpath('~/Dropbox/my_functions/')
addpath('../data_analysis/utils/')
load('../matrix/info_per_subj_final.mat');
load('../matrix/initial_fixations.mat');

[subjects, ~, subj_order] = unique({info_per_subj(:).subj});
subj_order  = subj_order'; % {info_per_subj.subj}
Nsubj       = length(unique(subj_order));
Ntr         = length(info_per_subj);

%% Figure 1. Experimental design and general behavior. 

%% A) Timeline of a trial. The search could be stopped by fixating into the 
% target region of after N saccades, where N varied across trials. 

%% B) Proportion of targets found as function of the number of saccades allowed. 
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

cols = rainbow_colors(Nsubj);
figure(2);clf
    set(gcf,'Color','w')
    hold on
%         for i=1:Nsubj
%             plot(NSACC_ALLOWED,p(i,:),'.-','Color',cols(i,:))
%         end
%         plot(NSACC_ALLOWED,nanmean(p,1),'k.-','LineWidth',3)
%         errorbar(NSACC_ALLOWED,nanmean(p,1),nanstd(p)./sqrt(sum(~isnan(n))),'k.-','LineWidth',3)
        boxplot(p,'Notch','on','Color','k','Labels',NSACC_ALLOWED)
    hold on
    box on
    set(gca,'XTick',NSACC_ALLOWED)
    xlabel('Number of saccades allowed')
    ylabel('P(target found)')
%     h=get(gca,'Children')
%     set(h.Children(strcmp({h.Children.Tag},'Outliers')),'Color','k')    
    
%% C) Distributions of saccade length for the first (red), the fifth 
% (green) and the tenth (blue) saccade. 
nsacc_used      = arrayfun(@(x) size(x.fixations,1),info_per_subj);
target_found    = [info_per_subj.target_found]';
[images,~,trialimage_id] = unique({info_per_subj.image});

nsacc = nan(max(trialimage_id),1);
for i = 1:max(trialimage_id)
    nsacc(i) = nanmedian(nsacc_used(target_found==1 & trialimage_id==i));
end
nsacc = nsacc(~isnan(nsacc));

figure(33); clf
    set(gcf,'Color','w')
    hold on
        xh = 1:12;
        [yh]=hist(nsacc,xh); yh=yh/sum(yh);
        h = bar(xh,yh);
    hold on
    box on
    set(h,'FaceColor',[.7 .7 .7])
    xlabel('Number of saccades to find the target)')
    ylabel('Prob.')


%% D) Mean (± s.e.m.) saccade length as function of saccade rank. 
% For this we're using the fixations instead the saccades itself. Check the
% pipeline to generate the saccade structure.
fixations   = [];
saccades    = [];
for tr = 1:Ntr
    nfix = size(info_per_subj(tr).fixations,1);
    if (nfix > 1)
        fixations       = [fixations; ...
                [repmat([tr subj_order(tr) nsacc_allowed(tr)],nfix,1), [1:nfix]',  info_per_subj(tr).fixations]];

        dx = info_per_subj(tr).fixations(2:end,1) - info_per_subj(tr).fixations(1:end-1,1);
        dy = info_per_subj(tr).fixations(2:end,2) - info_per_subj(tr).fixations(1:end-1,2);

        theta = atan(dy./dx); % atan:[-Inf,Inf]->[-pi/2,pi/2], y lo quiero en [0,2*pi] 
        theta(dx>=0 & dy>=0) = theta(dx>=0 & dy>=0);
        theta(dx<0 & dy>=0) = theta(dx<0 & dy>=0) + pi;
        theta(dx<0 & dy<0) = theta(dx<0 & dy<0) + pi;
        theta(dx>=0 & dy<0) = theta(dx>=0 & dy<0) + 2*pi;

        saccamp = sqrt(dx.^2 + dy.^2);

        saccades        = [saccades; ...
                [repmat([tr subj_order(tr) nsacc_allowed(tr)],nfix-1,1), [1:nfix-1]', saccamp, theta]];
    end
end

SACCRANK    = 1:12;
s = nan(Nsubj,length(SACCRANK));
for su = 1:Nsubj
    for i=1:length(SACCRANK)
        s(su,i) = nanmedian(saccades(saccades(:,2)==su & saccades(:,4)==SACCRANK(i),5));
    end
end


% C,D)
% SACCRANK = [1 5 10];
% cols = [1 0 0; 0 1 0; 0 0 1];
SACCRANK    = 1:12;
cols        = rainbow_colors(length(SACCRANK));
figure(3); clf;
    set(gcf,'Color','w')
%     subplot(2,1,1)
        hold on
            xh = 0:20:1000;
            for i=1:length(SACCRANK)
                yh = hist(saccades(saccades(:,4)==SACCRANK(i),5),xh); yh = yh/sum(yh);
                plot(xh,yh,'Color',cols(i,:))
            end
        hold off
        xlabel('Saccade amplitude (pxs)')
        ylabel('Prob.')
%     subplot(2,1,2)
%         hold on
%             boxplot(s,'Notch','on','Color','k','Labels',SACCRANK)
%         hold on
%         box on
%         set(gca,'XTick',SACCRANK)
%         xlabel('Saccade rank')
%         ylabel('Saccade amplitude (pxs)')

%% E) Distributions of saccade direction (measured in degrees from the 
% positive horizontal axis) for the first (red), the fifth (green) and the 
% tenth (blue) saccade. 

% SACCRANK = [1 5 10];
% cols = [1 0 0; 0 1 0; 0 0 1];
% figure(4);clf
%     set(gcf,'Color','w')
%     h = [];
%     for i = 1:length(SACCRANK)
%         subplot(1,3,i);
% %             [T,R]=rose(saccades(saccades(:,4)==SACCRANK(i),6)); R = R/length(saccades(saccades(:,4)==SACCRANK(i),6));
%             [R,T]=hist(saccades(saccades(:,4)==SACCRANK(i),6)); R = R/length(saccades(saccades(:,4)==SACCRANK(i),6));
%             R = [R R(1)];T = [T T(1)];
%             polarplot(T,R,'Color',cols(i,:),'LineWidth',1.5)
%     end

SACCRANK= 1:12;
cols    = rainbow_colors(length(SACCRANK));
figure(44);clf
    set(gcf,'Color','w')
    h = [];
        for i = 1:length(SACCRANK)
            [R,T]=hist(saccades(saccades(:,4)==SACCRANK(i),6)); R = R/length(saccades(saccades(:,4)==SACCRANK(i),6));
            R = [R R(1)];T = [T T(1)];
            polarplot(T,R,'Color',cols(i,:),'LineWidth',1.5)
            hold on
        end
    hold off

%% F) Spatial distributions of the first, the fifth and the tenth fixation 
% position. Warmed colors indicate larger probabilities.
addpath('utils/heatmap_code');
settings.filterSize = 300; % tama?o del filtro gaussiano en px.
settings.filterSD   = 30;  % SD del filtro gaussiano en px.
settings.max_nhist  = 50;  % controla el colormap. ajustarlo para que no sature el heatmap.
settings.alpha      = 0.8; % entre 0 y 1. cuanto mayor, mas se ven las fijaciones.
settings.grayscale  = 1;
im_file             = '../data_analysis/utils/white.jpg';
im                  = imread(im_file);  
[s1,s2]             = size(im);

FIXRANK    = 1:12;
M = zeros(length(FIXRANK),s1,s2);
for i=1:length(FIXRANK)
%     for su = 1:Nsubj
%         x = fixations(fixations(:,2)==su & fixations(:,4)==FIXRANK(i),5);
%         y = fixations(fixations(:,2)==su & fixations(:,4)==FIXRANK(i),6);
%         m = fun_heatmap(im_file, x, y, settings);
%         m = m/sum(m(:));
%     end
%     M(i,:,:) = squeeze(M(i,:,:)) + m;
    x = fixations(fixations(:,4)==FIXRANK(i),5);
    y = fixations(fixations(:,4)==FIXRANK(i),6);
    m = fun_heatmap(im_file, x, y, settings);
    M(i,:,:) = m/max(m(:));
end

cols = rainbow_colors(length(FIXRANK));
figure(5);clf
    set(gcf,'Color','w')
    set(gcf,'Position',[66 860 1852 90])
    for i = 1:length(FIXRANK)
        axes('Position',[0.05+0.075*(i-1),0.05,0.075,0.9])
            hold on
                imagesc(1-squeeze(M(i,:,:))); colormap gray
                plot([1 s2],[1 1],  '-','Color',cols(i,:),'LineWidth',3)
                plot([1 s2],[s1 s1],'-','Color',cols(i,:),'LineWidth',3)
                plot([1 1], [1 s1], '-','Color',cols(i,:),'LineWidth',3)
                plot([s2 s2],[1 s1],'-','Color',cols(i,:),'LineWidth',3)
            hold off
            axis tight
            set(gca,'Visible','off')
    end
       
%% All figs...
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
            
            imagesc(flipud(tmp)); colormap gray
            plot([1 s2],[1 1],  '-','Color','k')
            plot([1 s2],[s1 s1],'-','Color','k')
            plot([1 1], [1 s1], '-','Color','k')
            plot([s2 s2],[1 s1],'-','Color','k')
        hold off
        axis tight
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
            imagesc(flipud(im)); colormap gray
            plot([1 s2],[1 1],  '-','Color','k')
            plot([1 s2],[s1 s1],'-','Color','k')
            plot([1 1], [1 s1], '-','Color','k')
            plot([s2 s2],[1 s1],'-','Color','k')
        hold off
        axis tight
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

    xh = 0:20:1000;
    cols        = rainbow_colors(length(SACCRANK));
    for i=1:length(SACCRANK)
        axes('Position',[0.05+0.075*(i-1), 0.40, 0.0725, 0.15])
            yh = hist(saccades(saccades(:,4)==SACCRANK(i),5),xh); yh = yh/sum(yh);
            plot(xh,yh,'Color',cols(i,:),'LineWidth',1.5)
            set(gca,'YLim', [0 0.13],'XLim',[0 600])
            set(gca,'XTick',[100 500])
            box off
            grid on
            if (i==1)
                set(gca,'YTickLabel',[0 0.05 0.1])
                xlabel('Saccade amplitude (pxs)')
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
                imagesc(1-squeeze(M(i,:,:))); colormap gray
                plot([1 s2],[1 1],  '-','Color',cols(i,:),'LineWidth',3)
                plot([1 s2],[s1 s1],'-','Color',cols(i,:),'LineWidth',3)
                plot([1 1], [1 s1], '-','Color',cols(i,:),'LineWidth',3)
                plot([s2 s2],[1 s1],'-','Color',cols(i,:),'LineWidth',3)
            hold off
            axis tight
            set(gca,'Visible','off')
    end
