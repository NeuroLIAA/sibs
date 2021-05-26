clear all
close all
clc

%% Notas
% 1 - estamos corriendo el scanpathDistance no el MM
% 1 - algo extraño, la distribucion
% 1 - 
%
%% 1 - Calculo de métricas
% En cada imagen, el modelo correlation calcula la distancia promedio de su scanpath al resto. 
delta = 32;
min_fix = 3;
max_fix = 9;
image_size = [768 1024];

addpath('../../compare_models')
addpath('../../dynamic_models/utils')
addpath('../utils')
dynamic_model = 'correlation';

% % Se le aplica la grilla a los sujetos. Si ya fue aplicada, comentar.
% for img=1:134
%     if img ~= 132
%         path = char(strcat('../../matrix/images/info_per_subj_img_', num2str(img), '.mat')); fprintf('%s\n',path);
%         load(path)
%         info_per_subj = subjMapFixationToMatrix( info_per_subj, path, delta, image_size );
%     end
% end
            
ab = [2 round(2/3,4); 2 round(4/3,4); 2 2; 2 round(8/3,4); 2 round(16/5,4); 2 4; 2 round(16/3,4); 3 1; 3 2; 3 3; 3 4; 3 6; 3 8; 4 4; 4 8; 5 5; 5 8; 5 10; 3 round(24/5,4); 4 round(32/5,4); 
      4 round(4/3,4); 4 round(8/3,4); 4 round(16/3,4); 5 round(5/3,4); 5 round(10/3,4); 
      5 round(20/3,4); 4 round(32/3,4); 5 round(40/3,4)];

for k = 1:length(ab)
    ab_k = ab(k,:);
    a = ab_k(1);
    b = ab_k(2);

    mean_dist_model_corr = [];
    std_dist_model_corr = [];

    eliminados_model = [];
    adentro_model = [];
    mean_length_subj_corr = [];
    length_model = [];
    std_length_subj_corr = [];
    for img=1:134
        if img ~= 132
            aden = 0;
            elim = 0;

            path = char(strcat('../../out_models/deepgaze/',dynamic_model,...
                '/others/a_',num2str(a), '_b_', num2str(b), '_tam_celda_',...
                num2str(delta),'/scanpath/scanpath_', num2str(img), '.mat')); fprintf('%s\n',path);
            load(path)
            
            path = char(strcat('../../new_data/new_matrix/sinfo_img/info_per_img_', num2str(img), '.mat')); fprintf('%s\n',path);
            load(path)
            
            % Calculo fijaciones reducidas
            info_per_img = new_subjMapFixationToMatrix(info_per_img, '', delta, image_size);

            % Calculo la distancia de todos contra todos
            model_distance = [];
            length_subj_img = [];

            if length(scanpath(:,1)) <= max_fix && length(scanpath(:,1)) >= min_fix
                for subj_i=1:length(info_per_img)
                    if info_per_img(subj_i).target_found && ...
                            length(info_per_img(subj_i).fixations_matrix_reduced(:,1)) >= min_fix && ...
                            length(info_per_img(subj_i).fixations_matrix_reduced(:,1)) <= max_fix 
                        model_distance = [model_distance scanpathDistance(info_per_img(subj_i).fixations_matrix_reduced, [scanpath(:,2),scanpath(:,1)])];
                        length_subj_img = [length_subj_img length(info_per_img(subj_i).fixations_matrix_reduced(:,1))];
                    end
                end
            end

            if length(scanpath(:,1)) <= max_fix
                aden = aden + 1;
            else
                elim = elim + 1;
            end

            adentro_model = [adentro_model aden];
            eliminados_model = [eliminados_model elim];

            mean_dist_model_corr = [mean_dist_model_corr nanmean(model_distance)];
            std_dist_model_corr = [std_dist_model_corr nanstd(model_distance)];

            mean_length_subj_corr = [mean_length_subj_corr nanmean(length_subj_img)];
            std_length_subj_corr = [std_length_subj_corr nanstd(length_subj_img)];
            length_model = [length_model length(scanpath(:,1))];

        end
    end

    save_name = char(strcat('model_vs_subj_a_',num2str(a), '_b_', num2str(b), '_delta_',num2str(delta),'.mat'));
    save(save_name, 'adentro_model', 'eliminados_model', 'mean_dist_model_corr', 'std_dist_model_corr', 'mean_length_subj_corr', 'length_model', 'std_length_subj_corr')
    
end

%% 2 - Cálculo de métricas. 

ab = [2 round(2/3,4); 2 round(4/3,4); 2 2; 2 round(8/3,4); 2 round(16/5,4); 2 4; 2 round(16/3,4); 3 1; 3 2; 3 3; 3 4; 3 6; 3 8; 4 4; 4 8; 5 5; 5 8; 5 10; 3 round(24/5,4); 4 round(32/5,4); 
      4 round(4/3,4); 4 round(8/3,4); 4 round(16/3,4); 5 round(5/3,4); 5 round(10/3,4); 
      5 round(20/3,4); 4 round(32/3,4); 5 round(40/3,4)];

% minima cantidad de sujetos que encotraron el target   
nhumans = 15;
load('./subjVsSubj/subjVsSubj_delta_32_max_fix_9.mat')
delta = 32;
models = [];

ab_ind = 1;
for ab_value = ab' 
    a = ab_value(1);
    b = ab_value(2);
    fprintf('\n a: %d, b: %d \n',a,b) 
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
    
    models(ab_ind).a = a;
    models(ab_ind).b = b;
    models(ab_ind).ind = ind;
    models(ab_ind).p = p;
    models(ab_ind).zscore = zscore;
    models(ab_ind).length_zscore = length_zscore;
    
    eval(sprintf(char(strcat('ind_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = ind;'))))
    eval(sprintf(char(strcat('zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = zscore;'))))
    eval(sprintf(char(strcat('p_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = p;'))))
    eval(sprintf(char(strcat('length_zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = length_zscore;'))))
    eval(sprintf(char(strcat('zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),'_sorted = zscore_sorted;'))))
    eval(sprintf(char(strcat('p_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),'_sorted = p_sorted;'))))
    
    ab_ind = ab_ind+1;
end

%% 3 - Visualizar todos los modelos en una grilla

count_image_people_found = 0;
for i=1:133
    if sum(~isnan(mean_dist_img(i,:)))>nhumans
        count_image_people_found= count_image_people_found+1;
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
        title('Mediana de los Zscore de la distancia');
    %gráfico que muestre el total de aciertos del modelo
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
        title('Porcentaje de objetos encontrados por el modelo con respecto a los humanos')

%%

%% 4 - Mejores modelos segun zscore

min_found = 60;
model_min_found_mask = arrayfun(@(x) length(x.ind) >= min_found, models);

model_zscore_median = arrayfun(@(x) median(x.zscore), models);
[model_zscore_median_sorted, index_median_sorted] =  sort(abs(model_zscore_median));
best_ab_median = [];
for i=1:length(models)
    if model_min_found_mask(index_median_sorted(i))
        best_ab_median = [best_ab_median; [models(index_median_sorted(i)).a, models(index_median_sorted(i)).b,...
            index_median_sorted(i), model_zscore_median_sorted(i)]];
    end
end

model_zscore_mean = arrayfun(@(x) mean(x.zscore), models);
[model_zscore_mean_sorted, index_mean_sorted] =  sort(abs(model_zscore_mean));
best_ab_mean = [];
for i=1:length(models)
    if model_min_found_mask(index_mean_sorted(i))
        best_ab_mean = [best_ab_mean; [models(index_mean_sorted(i)).a, models(index_mean_sorted(i)).b,...
            index_mean_sorted(i), model_zscore_mean_sorted(i)]];
    end
end

%% 5 - Gráfico detallado de los mejores modelos

leg = {'a=3 b=8','a=5 b=8','a=3 b=4','a=4 b=8'}; 

figure(2); clf
    subplot(2,2,1)
        hold on
            title('Zscore de la distancia entre scanpath')
            plot(zscore_3_8_sorted,'-.or')
            plot(zscore_5_8_sorted,'-.ob')
            plot(zscore_3_4_sorted,'-.og')
            %plot(zscore_3_6_sorted,'k')
            plot(zscore_4_8_sorted,'-.om')
            %plot(zscore_5_16667_sorted,'m')
        hold off
        legend(leg,'Location','southeast');
        set(gca,'YLim',[-2 7],'YTick',-2:7)
        xlabel('imagenes')
        ylabel('zscore')
        grid on
    subplot(2,2,2)
        hold on
            title('Zscore de la distancia entre scanpath')
            xbins = -3:0.5:12;
            [yh,xh]=hist(zscore_3_8_sorted,xbins); plot(xh,yh,'r')
            [yh,xh]=hist(zscore_5_8_sorted,xbins); plot(xh,yh,'b')
            [yh,xh]=hist(zscore_3_4_sorted,xbins); plot(xh,yh,'g')
            [yh,xh]=hist(zscore_4_8_sorted,xbins); plot(xh,yh,'m')
            %[yh,xh]=hist(zscore_5_16667_sorted,xbins); plot(xh,yh,'r')
            set(gca,'XLim',[-2.5 8],'XTick',-2.5:0.5:8)
            xlabel('zscore')
            ylabel('cantidad de imágenes')
        hold off
        legend(leg);
        grid on
    subplot(2,2,3)
        boxX = [abs(zscore_3_8)'; abs(zscore_5_8)'; abs(zscore_3_4)'; abs(zscore_4_8)'];
        boxY = [ones(size(zscore_3_8')); 2*ones(size(zscore_5_8')); 3*ones(size(zscore_3_4')); 4*ones(size(zscore_4_8'))];
        boxplot(boxX, boxY);
        
        set(gca,'xticklabel',leg);
        title('Zscore de la distancia entre scanpth en valor absoluto');
        xlabel('modelos')
        ylabel('zscore')
        
    subplot(2,2,4)
        average = [median(abs(length_zscore_3_8)) median(abs(length_zscore_5_8)) median(abs(length_zscore_3_4)) median(abs(length_zscore_4_8))];
        bar(average);
        set(gca,'xticklabel',leg);
        title('mediana de los Zscore de la longitud');
        xlabel('modelos')
        ylabel('zscore')
    

figure(3); 
    subplot(1,2,1)
        hold on
            title('porcentaje de sujetos con mayor zscore que los modelos')
            plot(p_3_8_sorted,'-.or')
            plot(p_5_8_sorted,'-.ob')
            plot(p_3_4_sorted,'-.og')
            plot(p_4_8_sorted,'-.om')
            %plot(0.05*ones(size(p128_sorted)),'k--')
        hold off
        legend(leg,'Location','southeast');
        set(gca,'YLim',[0 1])
        xlabel('imagenes')
        ylabel('porcentaje')
        grid on
    subplot(1,2,2)
        hold on
            title('Zscore')
            [yh,xh]=hist(p_3_8_sorted,15); plot(xh,yh,'-.or')
            [yh,xh]=hist(p_5_8_sorted,15); plot(xh,yh,'-.ob')
            [yh,xh]=hist(p_3_4_sorted,15); plot(xh,yh,'-.og')
            [yh,xh]=hist(p_4_8_sorted,15); plot(xh,yh,'-.om')
            %[yh,xh]=hist(p_5_16667_sorted,15); plot(xh,yh,'r')
        hold off
        grid on
        legend(leg);

%% Figures
% figure(2); clf
%     subplot(2,2,1)
%         hold on
%             title('Zscore de la distancia entre scanpath')
%             plot(zscore_3_2_sorted,'m')
%             plot(zscore_3_3_sorted,'b')
%             plot(zscore_3_4_sorted,'g')
%             %plot(zscore_3_6_sorted,'k')
%             plot(zscore_4_26667_sorted,'c')
%             plot(zscore_5_16667_sorted,'r')
%         hold off
%         legend({'a=3 b=2','a=3 b=3','a=3 b=4','a=4 b=8/3','a=5 b=5/3'},'Location','southeast');
%         set(gca,'YLim',[-2 7],'YTick',-2:7)
%         xlabel('imagenes')
%         ylabel('zscore')
%         grid on
%     subplot(2,2,2)
%         hold on
%             title('Zscore de la distancia entre scanpath')
%             xbins = -3:0.5:12;
%             [yh,xh]=hist(zscore_3_2_sorted,xbins); plot(xh,yh,'m')
%             [yh,xh]=hist(zscore_3_3_sorted,xbins); plot(xh,yh,'b')
%             [yh,xh]=hist(zscore_3_4_sorted,xbins); plot(xh,yh,'g')
%             [yh,xh]=hist(zscore_4_26667_sorted,xbins); plot(xh,yh,'c')
%             [yh,xh]=hist(zscore_5_16667_sorted,xbins); plot(xh,yh,'r')
%             set(gca,'XLim',[-2.5 8],'XTick',-2.5:0.5:8)
%             xlabel('zscore')
%             ylabel('cantidad de imágenes')
%         hold off
%         legend({'a=3 b=2','a=3 b=3','a=3 b=4','a=4 b=8/3','a=5 b=5/3'});
%         grid on
%     subplot(2,2,3)
%         boxX = [abs(zscore_3_2)'; abs(zscore_3_3)'; abs(zscore_3_4)'; abs(zscore_4_26667)'; abs(zscore_5_16667)'];
%         boxY = [ones(size(zscore_3_4')); 2*ones(size(zscore_3_2')); 3*ones(size(zscore_3_3')); 4*ones(size(zscore_4_26667')); 5*ones(size(zscore_5_16667'))];
%         boxplot(boxX, boxY);
%         
%         set(gca,'xticklabel',{'a=3 b=2','a=3 b=3','a=3 b=4','a=4 b=8/3','a=5 b=5/3'});
%         title('Zscore de la distancia entre scanpth en valor absoluto');
%         xlabel('modelos')
%         ylabel('zscore')
%         
%     subplot(2,2,4)
%         average = [median(abs(length_zscore_3_2)) median(abs(length_zscore_3_3)) median(abs(length_zscore_3_4)) median(abs(length_zscore_4_26667)) median(abs(length_zscore_5_16667))];
%         bar(average);
%         set(gca,'xticklabel',{'a=3 b=2','a=3 b=3','a=3 b=4','a=4 b=8/3','a=5 b=5/3'});
%         title('mediana de los Zscore de la longitud');
%         xlabel('modelos')
%         ylabel('zscore')
%     
% 
% figure(3); 
%     subplot(1,2,1)
%         hold on
%             title('porcentaje de sujetos con mayor zscore que los modelos')
%             plot(p_3_2_sorted,'m')
%             plot(p_3_3_sorted,'b')
%             plot(p_3_4_sorted,'g')
%             plot(p_4_26667_sorted,'c')
%             plot(p_5_16667_sorted,'r')
%             %plot(0.05*ones(size(p128_sorted)),'k--')
%         hold off
%         legend({'a=3 b=2','a=3 b=3','a=3 b=4','a=4 b=8/3','a=5 b=5/3'},'Location','southeast');
%         set(gca,'YLim',[0 1])
%         xlabel('imagenes')
%         ylabel('porcentaje')
%         grid on
%     subplot(1,2,2)
%         hold on
%             title('Zscore')
%             [yh,xh]=hist(p_3_2_sorted,15); plot(xh,yh,'m')
%             [yh,xh]=hist(p_3_3_sorted,15); plot(xh,yh,'b')
%             [yh,xh]=hist(p_3_4_sorted,15); plot(xh,yh,'g')
%             [yh,xh]=hist(p_4_26667_sorted,15); plot(xh,yh,'c')
%             [yh,xh]=hist(p_5_16667_sorted,15); plot(xh,yh,'r')
%         hold off
%         grid on
%         legend({'a=3 b=2','a=3 b=3','a=3 b=4','a=4 b=8/3','a=5 b=5/3'});        

%% Ejemplo de por que da negativo el zscore
load('/home/usuario/vs_models/matrix/images/info_per_img_img_64.mat')

max = 0;
for i = 1:length(info_per_img)
    if info_per_img(i).target_found && length(info_per_img(i).fixations) < 10 && length(info_per_img(i).fixations) > max
        max = length(info_per_img(i).fixations);
    end
end

figure (4);
imshow('../../images/grayscale_35_oliva.jpg')
hold on
for i = 1:length(info_per_img)
    if info_per_img(i).target_found && length(info_per_img(i).fixations) < 10 && length(info_per_img(i).fixations) > 2 && i ~= 48
        N = max;
        cmap = jet(N);
        for n=1:(length(info_per_img(i).fixations))-1
            plot(info_per_img(i).fixations([n n+1],2),info_per_img(i).fixations([n n+1],1),'color', cmap(n,:))
        end
    end
    
end
load('/home/usuario/vs_models/out_models/mlnet/correlation/a_3_b_3_tam_celda_32/scanpath/scanpath_64.mat')
scanpath = (scanpath .* 32)-31;
plot(scanpath(:,2),scanpath(:,1),'r', 'linewidth',2)
rectangle('Position',[info_per_img(1).target_center(2)-36-10, info_per_img(1).target_center(1)-36-10,72+20,72+20], 'linewidth',2);

% %% Ejemplo de por que da negativo el zscore usando la grilla (queda medio raro
% load('/home/usuario/vs_models/matrix/images/info_per_subj_img_64.mat')
% 
% max = 0;
% for i = 1:length(info_per_subj)
%     if info_per_subj(i).target_found && length(info_per_subj(i).fixations_matrix_reduced(:,1)) < 10 && length(info_per_subj(i).fixations_matrix_reduced(:,1)) > max
%         max = length(info_per_subj(i).fixations_matrix_reduced(:,1));
%     end
% end
% 
% figure (4);
% imshow('../../images/grayscale_35_oliva.jpg')
% hold on
% for i = 1:length(info_per_subj)
%     subj = info_per_subj(i).fixations_matrix_reduced;
%     if info_per_subj(i).target_found && length(subj(:,1)) < 10 && length(subj(:,1)) > 2 && i ~= 48
%         N = max;
%         cmap = jet(N);
%         subj = (subj .* 32) - 16;
%         for n=1:(length(subj(:,1)))-1
%             plot(subj([n n+1],2),subj([n n+1],1),'color', cmap(n,:))
%         end
%     end
%     
% end
% load('/home/usuario/vs_models/out_models/mlnet/correlation/a_3_b_3_tam_celda_32/scanpath/scanpath_64.mat')
% scanpath = (scanpath .* 32)-16;
% plot(scanpath(:,2),scanpath(:,1),'r', 'linewidth',2)
% rectangle('Position',[info_per_subj(1).target_center(2)-36-10, info_per_subj(1).target_center(1)-36-10,72+20,72+20], 'linewidth',2);


%% Ejemplo de por que da positivo el zscore
load('/home/usuario/vs_models/matrix/images/info_per_subj_img_39.mat')

max = 0;
for i = 1:length(info_per_subj)
    if info_per_subj(i).target_found && length(info_per_subj(i).fixations) < 10 && length(info_per_subj(i).fixations) > max
        max = length(info_per_subj(i).fixations);
    end
end

figure (5);
imshow('../../images/grayscale_20_opensource.jpg')
hold on
for i = 1:length(info_per_subj)
    if info_per_subj(i).target_found && length(info_per_subj(i).fixations) < 10 && length(info_per_subj(i).fixations) > 2
        N = max;
        cmap = jet(N);
        for n=1:(length(info_per_subj(i).fixations))-1
            plot(info_per_subj(i).fixations([n n+1],2),info_per_subj(i).fixations([n n+1],1),'color', cmap(n,:))
        end
    end
    
end
load('/home/usuario/vs_models/out_models/mlnet/correlation/a_3_b_3_tam_celda_32/scanpath/scanpath_39.mat')
scanpath = (scanpath .* 32)-30;
plot(scanpath(:,2),scanpath(:,1),'r', 'linewidth',2)
rectangle('Position',[info_per_subj(1).target_center(2)-36-20, info_per_subj(1).target_center(1)-36-10,72+20,72+10], 'linewidth',2);

%% Ejemplo de por que da cerca de 0 el zscore


load('/home/usuario/vs_models/matrix/images/info_per_subj_img_2.mat')

max = 0;
for i = 1:length(info_per_subj)
    if info_per_subj(i).target_found && length(info_per_subj(i).fixations) < 10 && length(info_per_subj(i).fixations) > max
        max = length(info_per_subj(i).fixations);
    end
end

figure (6);  
imshow('../../images/grayscale_10_housebeautiful.jpg')
hold on
for i = 1:length(info_per_subj)
    if info_per_subj(i).target_found && length(info_per_subj(i).fixations) < 10 && length(info_per_subj(i).fixations) > 2 && i ~= 9
        N = max;
        cmap = jet(N);
        for n=1:(length(info_per_subj(i).fixations))-1
            plot(info_per_subj(i).fixations([n n+1],2),info_per_subj(i).fixations([n n+1],1),'color', cmap(n,:))
        end
    end
    
end
load('/home/usuario/vs_models/out_models/mlnet/correlation/a_3_b_4_tam_celda_32/scanpath/scanpath_2.mat')
scanpath = (scanpath .* 32)-31;
plot(scanpath(:,2),scanpath(:,1),'r', 'linewidth',2)
rectangle('Position',[info_per_subj(1).target_center(2)-36-15, info_per_subj(1).target_center(1)-36-15,72+15,72+20], 'linewidth',2);

%% Calculo las imágenes que aparecen en todos los modelos (DIO MAL)
% delta = 32;
% min_fix = 3;
% max_fix = 9;
% image_size = [768 1024];
% 
% addpath('../../compare_models')
% dynamic_model = 'correlation';
% 
% imagenes = zeros(1,134);
% 
% ab = [2 2; 3 1; 3 2; 3 3; 3 4; 3 6; 3 8; 4 4; 4 8; 5 5; 5 8; 5 10; 3 round(24/5,4); 4 round(32/5,4); 
%       4 round(4/3,4); 4 round(8/3,4); 4 round(16/3,4); 5 round(5/3,4); 5 round(10/3,4); 
%       5 round(20/3,4); 4 round(32/3,4); 5 round(40/3,4)];
%   
% for k = 1:length(ab)
%     ab_k = ab(k,:);
%     a = ab_k(1);
%     b = ab_k(2);
% 
%    for img=1:134
%         if img ~= 132
%            
%             path = char(strcat('../../out_models/mlnet/',dynamic_model,'/a_',num2str(a), '_b_', num2str(b), '_tam_celda_',num2str(delta),'/scanpath/scanpath_', num2str(img), '.mat')); fprintf('%s\n',path);
%             load(path);
% 
%             path = char(strcat('../../matrix/images/info_per_subj_img_', num2str(img), '.mat')); fprintf('%s\n',path);
%             load(path);
%             
%             cant_subj = 0;
%             if length(scanpath(:,1)) <= max_fix && length(scanpath(:,1)) >= min_fix
%                 for subj_i=1:length(info_per_subj)
%                     if info_per_subj(subj_i).target_found && ...
%                             length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) >= min_fix && ...
%                             length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) <= max_fix 
%                         cant_subj = cant_subj + 1;
%                     end
%                 end
%             end
% 
%             if cant_subj > 15
%                 imagenes(1, img) = imagenes(1, img) + 1;
%             end
%         end
%     end
% end
% 
% imagenesCompartidasPorLosModelos = [];
% 
% for k= 1:134
%     if imagenes(1,k) == 22
%         imagenesCompartidasPorLosModelos = [imagenesCompartidasPorLosModelos k];
%     end
% end
% 
% 
% % Plot de las imágenes que aparecen en todos los modelos
% delta = 32;
% min_fix = 3;
% max_fix = 9;
% image_size = [768 1024];
% 
% addpath('../../compare_models')
% dynamic_model = 'correlation';
% 
% ab = [2 2; 3 1; 3 2; 3 3; 3 4; 3 6; 3 8; 4 4; 4 8; 5 5; 5 8; 5 10; 3 round(24/5,4); 4 round(32/5,4); 
%       4 round(4/3,4); 4 round(8/3,4); 4 round(16/3,4); 5 round(5/3,4); 5 round(10/3,4); 
%       5 round(20/3,4); 4 round(32/3,4); 5 round(40/3,4)];
%   
% for k = 1:length(ab)
%     ab_k = ab(k,:);
%     a = ab_k(1);
%     b = ab_k(2);
% 
%     mean_dist_model_corr = [];
%     std_dist_model_corr = [];
% 
%     eliminados_model = [];
%     adentro_model = [];
%     mean_length_subj_corr = [];
%     length_model = [];
%     std_length_subj_corr = [];
%     for img_index=1:length(imagenesCompartidasPorLosModelos)
%         img = imagenesCompartidasPorLosModelos(1,img_index);
%         if img ~= 132
%             aden = 0;
%             elim = 0;
% 
%             path = char(strcat('../../out_models/mlnet/',dynamic_model,'/a_',num2str(a), '_b_', num2str(b), '_tam_celda_',num2str(delta),'/scanpath/scanpath_', num2str(img), '.mat')); fprintf('%s\n',path);
%             load(path)
% 
%             path = char(strcat('../../matrix/images/info_per_subj_img_', num2str(img), '.mat')); fprintf('%s\n',path);
%             load(path)
% 
%             % Calculo la distancia de todos contra todos
%             model_distance = [];
%             length_subj_img = [];
% 
%             if length(scanpath(:,1)) <= max_fix && length(scanpath(:,1)) >= min_fix
%                 for subj_i=1:length(info_per_subj)
%                     if info_per_subj(subj_i).target_found && ...
%                             length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) >= min_fix && ...
%                             length(info_per_subj(subj_i).fixations_matrix_reduced(:,1)) <= max_fix 
%                         model_distance = [model_distance scanpathDistance(info_per_subj(subj_i).fixations_matrix_reduced,scanpath)];
%                         length_subj_img = [length_subj_img length(info_per_subj(subj_i).fixations_matrix_reduced(:,1))];
%                     end
%                 end
%             end
% 
%             if length(scanpath(:,1)) <= max_fix
%                 aden = aden + 1;
%             else
%                 elim = elim + 1;
%             end
% 
%             adentro_model = [adentro_model aden];
%             eliminados_model = [eliminados_model elim];
% 
%             mean_dist_model_corr = [mean_dist_model_corr mean(model_distance)];
%             std_dist_model_corr = [std_dist_model_corr std(model_distance)];
% 
%             mean_length_subj_corr = [mean_length_subj_corr mean(length_subj_img)];
%             std_length_subj_corr = [std_length_subj_corr std(length_subj_img)];
%             length_model = [length_model length(scanpath(:,1))];
% 
%         end
%     end
% 
%     save_name = char(strcat('model_vs_subj_a_',num2str(a), '_b_', num2str(b), '_delta_',num2str(delta),'_reducido.mat'));
%     save(save_name, 'adentro_model', 'eliminados_model', 'mean_dist_model_corr', 'std_dist_model_corr', 'mean_length_subj_corr', 'length_model', 'std_length_subj_corr')
% end
% 
% 
% nhumans = 15;
% load('../subjVsSubj/subjVsSubj_delta_32.mat')
% 
% 
% for k = 1:length(ab)
%     ab_k = ab(k,:);
%     a = ab_k(1);
%     b = ab_k(2);
% 
%     
%     load_name = char(strcat('model_vs_subj_a_',num2str(a), '_b_', num2str(b), '_delta_',num2str(delta),'_reducido.mat'));
%     load(load_name)
%     
%     zscore = [];
% 
%     for img_index=1:length(imagenesCompartidasPorLosModelos)
%         i = imagenesCompartidasPorLosModelos(1,img_index);
%     
%         if sum(~isnan(mean_dist_img(i,:)))>nhumans && ~isnan(mean_dist_model_corr(1,img_index))
%             zscore        = [zscore ((mean_dist_model_corr(1,img_index) - nanmean(mean_dist_img(i,:)))/nanstd(mean_dist_img(i,:)))];
%         end
%     end
%     eval(sprintf(char(strcat('zscore_' ,strrep(num2str(a),'.',''),'_' ,strrep(num2str(b),'.',''),' = zscore;'))))
% end
% 
% figure(4); 
%     %quedarse solo con las imágenes qeu aparecen en todos los modelos
%     set(gcf,'Position',[370 400 435 515])
% 
%         average = [median(abs(zscore_2_2))     median(abs(zscore_2_2))     median(abs(zscore_2_2)) median(abs(zscore_2_2))     median(abs(zscore_2_2))  median(abs(zscore_2_2))  median(abs(zscore_2_2));
%                    median(abs(zscore_3_1))     median(abs(zscore_3_2))     median(abs(zscore_3_3)) median(abs(zscore_3_4))     median(abs(zscore_3_48)) median(abs(zscore_3_6))  median(abs(zscore_3_8)); 
%                    median(abs(zscore_4_13333)) median(abs(zscore_4_26667)) median(abs(zscore_4_4)) median(abs(zscore_4_53333)) median(abs(zscore_4_64)) median(abs(zscore_4_8))  median(abs(zscore_4_106667)); 
%                    median(abs(zscore_5_16667)) median(abs(zscore_5_33333)) median(abs(zscore_5_5)) median(abs(zscore_5_66667)) median(abs(zscore_5_8))  median(abs(zscore_5_10)) median(abs(zscore_5_133333))];
%         imagesc(average);
%         colormap bone
%         c = colorbar;
%         c.Label.String = 'Zscore';
%         xlabel('b/a')
%         ylabel('1/a')
%         set(gca,'xticklabel',{'0,33','0,66','1','1,33','1.60','2','2,67'},'XTick',1:7);
%         set(gca,'yticklabel',{'1/2','1/3','1/4','1/5'},'YTick',1:4);
%         title('mediana de los Zscore de la distancia');
%    