%% Cargar mapas de saliencia

clear all
close all
clc

addpath('../data_analysis/utils/')
addpath('../data_analysis/utils/heatmap_code/')
load('../matrix/all_images.mat')
load('../matrix/initial_fixations.mat')
targetLocations = load('../matrix/target_positions_filtered.mat');
load('../data_extraction/image_names.mat');
addpath('../saliency/deepgaze/')

delta = 32;
modo  = 'mean';

dd=dir('../matrix/images/*.mat'); filenames_img = {dd.name}'; Nimg = length(filenames_img);

%% Evaluación con una imagen
fileName   = 'grayscale_70_oliva.jpg';%;'grayscale_100_oliva.jpg';
imagen     = imread(strcat('../images/',fileName));
saliency   = imread(fileName);
%ind = strcmp({targetLocations.target_positions.image},'grayscale_100_oliva.jpg');
%targetInfo = targetLocations.target_positions(ind);
%[imageNrow, imageNcol] = size(s);
%mm = templateMaskReduce(targetInfo, imageNcol, imageNrow, delta);
%imagesc(mm)
indice      = strcmp({targetLocations.target_positions.image}, fileName);
tInfo       = targetLocations.target_positions(indice);
mysearch    = saliencySearcher(saliency, tInfo, delta, modo, 0);
mysearch.scanpath(1:10,:)
disp(mysearch.nfix)

%% Dibujar la solucion

figure(1);
    subplot(2,2,1);
    plot(mysearch.scanpath(1:mysearch.nfix,2), mysearch.scanpath(1:mysearch.nfix,1))
    set(gca, 'YDir','reverse')
    hold on
    scatter(mysearch.scanpath(1:mysearch.nfix,2), mysearch.scanpath(1:mysearch.nfix,1))
    txt = 'Init Fix';
    text(mysearch.scanpath(1,2),mysearch.scanpath(1,1),txt);
    txt = 'Last Fix';
    text(mysearch.scanpath(mysearch.nfix,2),mysearch.scanpath(mysearch.nfix,1),txt);
    [yl, xl] = size(mysearch.rsaliency);
    xlim([0 xl]);
    ylim([0 yl]);
    
    subplot(2,2,2);
    imshow(imagen)
    % le puedo agregar las cosas del target
    
    subplot(2,2,3);
    imagesc(mysearch.rsaliency)
    
    subplot(2,2,4);
    imshow(mysearch.template)

%% Corrida sobre todas las imagenes

nfixNullModel = zeros(Nimg,1);
for i=1:Nimg
    
    fileName         = image_names{i};
    %imagen           = imread(strcat('../images/',fileName));
    saliency         = imread(fileName);
    indice           = strcmp({targetLocations.target_positions.image}, fileName);
    tInfo            = targetLocations.target_positions(indice);
    mysearch         = saliencySearcher(saliency, tInfo, delta, modo, 0);
    nfixNullModel(i,1) = mysearch.nfix;
    
end

%% Evaluacion del searcher basico

figure()
histogram(nfixNullModel)
%xlim([0,15])

%% Prueba con busqueda con IOR

fileName   = 'grayscale_99_oliva.jpg';%;'grayscale_70_oliva.jpg';
imagen     = imread(strcat('../images/',fileName));
saliency   = imread(fileName);
indice      = strcmp({targetLocations.target_positions.image}, fileName);
tInfo       = targetLocations.target_positions(indice);
inhFactor    = 1/3;
inhSize      = 2;
mysearch2    = saliencyIoRSearcher(saliency, tInfo, delta, modo, inhFactor, inhSize, 0);
mysearch2.scanpath(1:10,:)
disp(mysearch2.nfix)

%% Graficamos IOR

figure(1);
    subplot(2,2,1);
    plot(mysearch2.scanpath(1:mysearch2.nfix, 2), mysearch2.scanpath(1:mysearch2.nfix, 1))
    set(gca, 'YDir','reverse')
    hold on
    scatter(mysearch2.scanpath(1:mysearch2.nfix, 2), mysearch2.scanpath(1:mysearch2.nfix, 1))
    txt = 'Init Fix';
    text(mysearch2.scanpath(1,2), mysearch2.scanpath(1,1), txt);
    txt = 'Last Fix';
    text(mysearch2.scanpath(mysearch2.nfix,2), mysearch2.scanpath(mysearch2.nfix,1), txt);
    [yl, xl] = size(mysearch2.rsaliency);
    xlim([0 xl]);
    ylim([0 yl]);
    
    subplot(2,2,2);
    imshow(imagen)
    % le puedo agregar las cosas del target
    
    subplot(2,2,3);
    imagesc(mysearch2.rsaliency)
    
    subplot(2,2,4);
    imshow(mysearch2.template)

%% 

nfixNullModelIOR = zeros(Nimg,1);
inhFactor    = 1/3;
inhSize      = 2;
for i=1:Nimg
    
    fileName         = image_names{i};
    %imagen           = imread(strcat('../images/',fileName));
    saliency         = imread(fileName);
    image_size       = size(saliency);
    if image_size(1) ~= 768 || image_size(2) ~= 1024
        disp(image_size)
    end
    indice           = strcmp({targetLocations.target_positions.image}, fileName);
    tInfo            = targetLocations.target_positions(indice);
    indice_initFix   = strcmp({initial_fixations.image}, fileName);
    initFix          = mapToReducedMatrix(initial_fixations(indice).initial_fix,...
        delta,image_size);
    mysearch         = saliencyIoRSearcher(saliency, tInfo, delta, modo, inhFactor,...
        inhSize, initFix, 0);
    nfixNullModelIOR(i,1) = mysearch.nfix;
    fname = sprintf('../out_models/saliency-based/correlation/a_3_b_4_tam_celda_32/scanpath/scanpath_%d.mat',i);
    scanpath = mysearch.scanpath(1:mysearch.nfix+1, :);
    save(fname, 'scanpath');
end

%%
figure()
histogram(nfixNullModelIOR)
%xlim([0,15])