function searchInfo = saliencyIoRSearcher(saliencyMap, targetInfo, delta, ...
                                    mode, inhibitionFactor, inhibitionSize,...
                                    initFix, verbose)
    
    % IoR inhibición de retorno
    [Row, Col]   = size(saliencyMap);
    S            = reduceMatrix(saliencyMap, delta, mode);
    [Nrow, Ncol] = size(S);
    searchMatrix = zeros(Nrow,Ncol);
    notFound     = true;
    Nfix         = 1;
    
    % indices de la matriz reducida
    currentFix   = [0 0];
    
    % inicializo el struct de salida
    searchInfo.scanpath  = zeros(Nrow * Ncol,2);
    searchInfo.found     = false;
    searchInfo.nfix      = 0;
    searchInfo.rsaliency = S;
    
    % creamos el template
    template             = templateMaskReduce(targetInfo, Col, Row, delta);
    
    % harcodeado la distancia mínima
    minDist     = zeros(Nrow*Ncol);
    minDist(2)  = 3;
    minDist(3)  = 2;
    minDist(4)  = 1;
    while notFound
        
        % ordenamos saliencia
        [saliencyValue, saliencyIndex]  = sort(S(:), 'descend');
        [I, J]                          = ind2sub(size(S), saliencyIndex);
        
        % conseguimos la mejor fijacion
        [nextFixLocation, searchMatrix] = getNextFix(Nfix, currentFix, saliencyValue, I, J, searchMatrix, minDist);
        currentFix = nextFixLocation;
        searchInfo.scanpath(Nfix,:) = currentFix;
        Nfix = Nfix + 1;
        
        % aplicamos la inhibicion de retorno
        S = inhibitePatch(S, currentFix, inhibitionFactor, inhibitionSize);
        
        % verbosidad para chequeo
        % TODO: agregar verbosidad con impresiones de los distintos mapas a
        % lo largo de las fijaciones
        if verbose==1
            fprintf('Nfix: ')
            disp(Nfix)
            fprintf('CurrentFix: ')
            disp(currentFix)
        end
        
        % chequeo si lo encontre
        if template(currentFix(1),currentFix(2)) 
            notFound = false;
            % hasta aca fix location solo devuelve la posicion en la matriz
            % reducida, tendria que generarme la ubicacion para mostrar en
            % la imagen original haciendo la transformacion de coordenas
            searchInfo.nfix = Nfix - 1;
            searchInfo.fixposx = currentFix(1);
            searchInfo.fixposy = currentFix(2);
        end       
    end
    % agrego la primer fijacion
    searchInfo.scanpath = [initFix; searchInfo.scanpath]; 
    searchInfo.template = template;
    searchInfo.found = not(notFound);
end

function [nextFixLocation, searchMatrix] = getNextFix(Nfix, currentFix, saliencyValue, I, J, searchMatrix, minDistance)
    % recorrer saliencyValue
    % minDistance pensado como 0, 2, 1, 1, 0 ...
    ind = 1;
    escapeCond = true;
    while (ind < numel(saliencyValue)) && escapeCond
        i = I(ind);
        j = J(ind);
        % si todavia no lo visite
        if not(searchMatrix(i,j))
            % las condiciones las dividi en dos solo para que no me
            % molesten
            % keyboard
            if (max(abs(currentFix-[i,j])) >= minDistance(Nfix))
                % marco que lo estoy visitando y salgo
                % disp(Nfix);
                % keyboard
                searchMatrix(i,j) = 1;
                escapeCond = false;
                nextFixLocation = [i j];
            else
                ind = ind + 1;
            end
        else
            % si no cumple con la condicion de no haberlo visitado y de
            % tener una distancia minima, sigo recorriendo la imagen
            ind = ind + 1;
            %disp(ind)
        end
    end
end

function S = inhibitePatch(S, currentFix, inhFactor, inhSize)
    
    [nrow, ncol] = size(S);
    
    row     = currentFix(1);
    col     = currentFix(2);
    patch = S(max([1,row-inhSize]):min([nrow,row+inhSize]), max([1,col-inhSize]):min([ncol, col+inhSize]));
    patch = patch * inhFactor;
    S(max([1,row-inhSize]):min([nrow,row+inhSize]), max([1,col-inhSize]):min([ncol, col+inhSize])) = patch;
    
end