function metrics = evalSaliencyMap(saliencyMap, imagesSample, images, fixations, baselineMap, otherMap, summary)
    % Función para para calcular AUC, NSS, KLD e InfoGain
    % saliencyMap = Nombre del modelo
    % imagesSample = Lista con los indices de las imagenes a usar
    % images = Nombres de los archivos de las imagenes
    % fixations = Fijaciones cargadas con "fixated area"
    if nargin == 3
        % images, fixations siempre se asumen dadas
        saliencyMap = 'deepgaze';
    end
    if nargin < 4
        imagesSample  = [1:134];
    end
    features = loadFixations(saliencyMap, imagesSample, images);
    Nsamp = length(imagesSample); 
    metrics = nan(Nsamp, 6);
    for img = 1:Nsamp
        metrics(img,1) = AUC_Judd(features(:,:,img), fixations(:,:,imagesSample(img)));
        metrics(img,2) = AUC_shuffled(features(:,:,img), fixations(:,:,imagesSample(img)), otherMap);
        metrics(img,3) = InfoGain(features(:,:,img), fixations(:,:,imagesSample(img)), baselineMap);
        metrics(img,4) = NSS(features(:,:,img), fixations(:,:,imagesSample(img)));
        metrics(img,5) = KLdiv(features(:,:,img), fixations(:,:,imagesSample(img)));
        metrics(img,6) = AUC_Borji(features(:,:,img), fixations(:,:,imagesSample(img)));
    end
    if summary
        metrics = mean(metrics,1);
    end
end