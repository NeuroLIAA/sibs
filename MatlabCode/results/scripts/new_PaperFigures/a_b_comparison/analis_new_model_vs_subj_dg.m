clear all
close all
clc

addpath('../../compare_models')
addpath('../../dynamic_models/utils')
        
%% En cada imagen, el modelo correlation calcula la distancia promedio de su scanpath al resto. 
if 0
%     fun_genero_new_model_vs_subj_file(dynamic_model,ab_value,delta,min_fix,max_fix,image_size)
    fun_genero_new_model_vs_subj_file_dg('correlation',[3 4],32,3,4,[768 1024])
end

%% Plot en cada imagen, 
% el modelo correlation calcula la distancia promedio de su scanpath al resto. 
ab_value        = [3 4]; 
nhumans         = 15;
delta           = 32;
load(sprintf('../../TesisFigures/subjVsSubj/subjVsSubj_delta_%d.mat',delta))
a = ab_value(1);
b = ab_value(2);
%%
% dynamic_model   = 'correlation';
%     
%     load_name = char(strcat('model_vs_subj_',dynamic_model,'_a_',num2str(a), '_b_', num2str(b), '_delta_',num2str(delta),'.mat'));
%     load(load_name)
%     
%     ind           = [];
%     zscore        = [];
%     p             = [];
%     length_zscore = [];
% 
%     for i=1:133
%         if sum(~isnan(mean_dist_img(i,:)))>nhumans && ~isnan(mean_dist_model_corr(1,i))
%             ind           = [ind i];
%             zscore        = [zscore ((mean_dist_model_corr(1,i) - nanmean(mean_dist_img(i,:)))/nanstd(mean_dist_img(i,:)))];
%             p             = [p sum(mean_dist_model_corr(1,i) < mean_dist_img(i,~isnan(mean_dist_img(i,:))))/sum(~isnan(mean_dist_img(i,:)))];
%             length_zscore = [length_zscore ((length_model(1,i) - mean_length_subj_corr(1,i))/std_length_subj_corr(1,i))];
% 
%         end
%     end
%     
%     zscore_sorted = sort(zscore);
%     p_sorted      = sort(p);
%     
%     eval(sprintf(char(strcat('Mgreedy.ind_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = ind;'))))
%     eval(sprintf(char(strcat('Mgreedy.zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = zscore;'))))
%     eval(sprintf(char(strcat('Mgreedy.p_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = p;'))))
%     eval(sprintf(char(strcat('Mgreedy.length_zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = length_zscore;'))))
%     eval(sprintf(char(strcat('Mgreedy.zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),'_sorted = zscore_sorted;'))))
%     eval(sprintf(char(strcat('Mgreedy.p_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),'_sorted = p_sorted;'))))
%     eval(sprintf(char(strcat('Mgreedy.mean_dist_img_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = mean_dist_img;'))))
%%
dynamic_model   = 'correlation';
    load_name = char(strcat('model_vs_subj_a_',num2str(a), '_b_', num2str(b), '_delta_',num2str(delta),'.mat'));
    load(load_name)
    
    ind           = [];
    zscore        = [];
    p             = [];
    length_zscore = [];

    for i=1:133
        if sum(~isnan(mean_dist_img(i,:)))>nhumans && ~isnan(mean_dist_model_corr(1,i))
            ind           = [ind i];
            zscore        = [zscore ((mean_dist_model_corr(1,i) - nanmean(mean_dist_img(i,:)))/nanstd(mean_dist_img(i,:)))];
            p             = [p sum(mean_dist_model_corr(1,i) < mean_dist_img(i,~isnan(mean_dist_img(i,:))))/sum(~isnan(mean_dist_img(i,:)))];
            length_zscore = [length_zscore ((length_model(1,i) - mean_length_subj_corr(1,i))/std_length_subj_corr(1,i))];

        end
    end
    
    zscore_sorted = sort(zscore);
    p_sorted      = sort(p);
       
    eval(sprintf(char(strcat('Mcorrel.ind_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = ind;'))))
    eval(sprintf(char(strcat('Mcorrel.zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = zscore;'))))
    eval(sprintf(char(strcat('Mcorrel.p_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = p;'))))
    eval(sprintf(char(strcat('Mcorrel.length_zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = length_zscore;'))))
    eval(sprintf(char(strcat('Mcorrel.zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),'_sorted = zscore_sorted;'))))
    eval(sprintf(char(strcat('Mcorrel.p_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),'_sorted = p_sorted;'))))
    eval(sprintf(char(strcat('Mcorrel.mean_dist_img_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = mean_dist_img;'))))

    
%%    
    
figure(2); clf
    subplot(2,2,1)
        hold on
            title('Zscore de la distancia entre scanpath')
            plot(Mgreedy.zscore_3_4_sorted,'b')
            plot(Mcorrel.zscore_3_4_sorted,'r')
        hold off
        legend({'Greedy','Correlation'});
        set(gca,'YLim',[-2 7],'YTick',-2:7)
        xlabel('imagenes')
        ylabel('zscore')
        grid on
    subplot(2,2,2)
        hold on
            title('Zscore de la distancia entre scanpath')
            xbins = -3:0.5:12;
            [yh,xh]=hist(Mgreedy.zscore_3_4_sorted,xbins); plot(xh,yh,'b')
            [yh,xh]=hist(Mcorrel.zscore_3_4_sorted,xbins); plot(xh,yh,'r')
            set(gca,'XLim',[-2.5 8],'XTick',-2.5:0.5:8)
            xlabel('zscore')
            ylabel('cantidad de imágenes')
        hold off
        legend({'Greedy','Correlation'});
        grid on
    subplot(2,2,3)
        boxX = [abs(Mgreedy.zscore_3_4)'; abs(Mcorrel.zscore_3_4)'];
        boxY = [ones(size(Mgreedy.zscore_3_4')); 2*ones(size(Mcorrel.zscore_3_4'))];
        boxplot(boxX, boxY);
        
        set(gca,'xticklabel',{'Greedy','Correlation'});
        title('Zscore de la distancia entre scanpth en valor absoluto');
        xlabel('modelos')
        ylabel('zscore')
        
    subplot(2,2,4)
        average = [median(abs(Mgreedy.length_zscore_3_4)) median(abs(Mcorrel.length_zscore_3_4))];
        bar(average);
        set(gca,'xticklabel',{'Greedy','Correlation'});
        title('mediana de los Zscore de la longitud');
        xlabel('modelos')
        ylabel('zscore')
%     

figure(3); 
    subplot(1,2,1)
        hold on
            title('porcentaje de sujetos con mayor zscore que los modelos')
            plot(Mgreedy.p_3_4_sorted,'b')
            plot(Mcorrel.p_3_4_sorted,'r')
            %plot(0.05*ones(size(p128_sorted)),'k--')
        hold off
        legend({'Greedy','Correlation'});
        set(gca,'YLim',[0 1])
        xlabel('imagenes')
        ylabel('porcentaje')
        grid on
    subplot(1,2,2)
        hold on
            title('Zscore')
            [yh,xh]=hist(Mgreedy.p_3_4_sorted,15); plot(xh,yh,'b')
            [yh,xh]=hist(Mcorrel.p_3_4_sorted,15); plot(xh,yh,'r')
        hold off
        grid on
        legend({'Greedy','Correlation'});

%%
Mgreedy.count_image_people_found = 0;
for i=1:133
    if sum(~isnan(Mgreedy.mean_dist_img(i,:)))>nhumans
        Mgreedy.count_image_people_found= count_image_people_found+1;
    end
end

Mcorrel.count_image_people_found = 0;
for i=1:133
    if sum(~isnan(Mcorrel.mean_dist_img(i,:)))>nhumans
        Mcorrel.count_image_people_found= count_image_people_found+1;
    end
end

figure(1); 
    %quedarse solo con las imágenes qeu aparecen en todos los modelos
    set(gcf,'Position',[370 400 435 515])
    subplot(2,1,1)
        average = [median(abs(zscore_2_06667)) median(abs(zscore_2_13333)) median(abs(zscore_2_2)) median(abs(zscore_2_26667)) median(abs(zscore_2_32)) median(abs(zscore_2_4))  median(abs(zscore_2_53333));
                   median(abs(zscore_3_1))     median(abs(zscore_3_2))     median(abs(zscore_3_3)) median(abs(zscore_3_4))     median(abs(zscore_3_48)) median(abs(zscore_3_6))  median(abs(zscore_3_8)); 
                   median(abs(zscore_4_13333)) median(abs(zscore_4_26667)) median(abs(zscore_4_4)) median(abs(zscore_4_53333)) median(abs(zscore_4_64)) median(abs(zscore_4_8))  median(abs(zscore_4_106667)); 
                   median(abs(zscore_5_16667)) median(abs(zscore_5_33333)) median(abs(zscore_5_5)) median(abs(zscore_5_66667)) median(abs(zscore_5_8))  median(abs(zscore_5_10)) median(abs(zscore_5_133333))];
        imagesc(average);
        colormap bone
        c = colorbar;
        c.Label.String = 'Zscore';
        xlabel('b/a')
        ylabel('1/a')
        set(gca,'xticklabel',{'0,33','0,66','1','1,33','1.60','2','2,67'},'XTick',1:7);
        set(gca,'yticklabel',{'1/2','1/3','1/4','1/5'},'YTick',1:4);
        title('mediana de los Zscore de la distancia');
    %hacer un gráfico que muestre el % de aciertos del modelo
    %en funcion de cuantas imagenes hay que por lo menos 15
    %personas hayan encontrado el target. 
    subplot(2,1,2)
        siz = [length(zscore_2_06667) length(zscore_2_13333) length(zscore_2_2) length(zscore_2_26667) length(zscore_2_32) length(zscore_2_4)  length(zscore_2_53333);
               length(zscore_3_1)     length(zscore_3_2)     length(zscore_3_3) length(zscore_3_4)     length(zscore_3_48) length(zscore_3_6)  length(zscore_3_8); 
               length(zscore_4_13333) length(zscore_4_26667) length(zscore_4_4) length(zscore_4_53333) length(zscore_4_64) length(zscore_4_8)  length(zscore_4_106667); 
               length(zscore_5_16667) length(zscore_5_33333) length(zscore_5_5) length(zscore_5_66667) length(zscore_5_8)  length(zscore_5_10) length(zscore_5_133333)];

        siz = (siz .* 100) ./ count_image_people_found;
        imagesc(siz);
        colormap bone
        c = colorbar;
        c.Label.String = '% imgs';
        xlabel('b/a')
        ylabel('1/a')
        set(gca,'xticklabel',{'0,33','0,66','1','1,33','1.60','2','2,67'},'XTick',1:7);
        set(gca,'yticklabel',{'1/2','1/3','1/4','1/5'},'YTick',1:4);
        title('porcentaje de objetos encontrados por el modelo con respecto a los humanos')
                