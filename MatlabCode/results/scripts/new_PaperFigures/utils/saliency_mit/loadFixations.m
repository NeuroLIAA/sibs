function features = loadFixations(saliencyMap, imagesSample, images)
    % Función para cargar las fijaciones para las metricas del MIT 
    % saliencyMap = Nombre del modelo
    % imagesSample = Lista con los indices de las imagenes a usar
    % images = Nombres de los archivos de las imagenes
    if nargin == 1
        saliencyMap = 'deepgaze';
    end
    if nargin < 3
        imagesSample  = [1:134];
    end
    modelos_list    = {'mlnet', 'sam_resnet', 'sam_vgg', 'deepgaze', 'icf', 'humans', 'center'};
    %modelos_name    = {'MLNet', 'SAM-ResNet', 'SAM-VGG', 'DeepGaze', 'icf', 'Humans', 'Center'};
    %modelos_col     = {'r', 'b', 'c', 'g', 'm', 'k', [.7 .7 .7]}; 
    %modelo_ind = 4;
    modelo_ind = find(strcmp(saliencyMap,modelos_list));
    features = nan(768,1024,length(imagesSample));
    for im = 1:length(imagesSample)
        if strcmp(modelos_list{modelo_ind}, 'center')
            features(:,:,im) = imread('../saliency/center.jpg');
        elseif strcmp(modelos_list{modelo_ind}, 'humans')
            features(:,:,im) = imread(['../saliency/humans_fix_3_to_3/' images{imagesSample(im)}]);
        else
            try
                features(:,:,im) = imread(['../saliency/' modelos_list{modelo_ind} '/' images{imagesSample(im)}]);
            catch
                keyboard
            end
        end
    end
end