function reducedMask = templateMaskReduce(targetInfo, imageNcol, imageNrow, delta)
    % - image es el string del nombre del archivo de la imagen
    % - targetInfo struct que contiene la informacion de la ubicacion del 
    % target para la imagen image. Los nombres de los campos deberan ser
    % matched_column y matched_row
    mask = zeros(imageNrow, imageNcol);
    % actualizo los lugares donde esta el target en la mascara
    mask(targetInfo.matched_row:targetInfo.matched_row + targetInfo.template_columns, ...
        targetInfo.matched_column:targetInfo.matched_column + targetInfo.template_side_length) = 1;
    reducedMask = reduceMatrix(mask, delta, 'max');        
end