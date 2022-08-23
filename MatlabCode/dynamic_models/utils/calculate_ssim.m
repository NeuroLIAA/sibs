function ssimval = calculate_ssim(img,template)
    % add padding
    img2 = zeros(size(img,1)+size(template,1),size(img,2)+size(template,2));
    img2 = mean(img(:))*ones(size(img,1)+size(template,1),size(img,2)+size(template,2));
    img2(   size(template,1)/2 + [1:size(img,1)],...
            size(template,2)/2 + [1:size(img,2)] ) = img;
    img2 = uint8(img2);
    
    tic
    ssimval = nan(size(img));
    k = 1; K = size(ssimval,1)*size(ssimval,2);
    for j = 1:size(ssimval,2)
        for i = 1:size(ssimval,1)
            tszx = size(template,1)/2;
            tszy = size(template,2)/2;
            indi = i + tszx + [-tszx:(tszx-1)];
            indj = j + tszy + [-tszy:(tszy-1)];
            % [ssimval, ssimmap] = ssim(A,ref);
            ssimval(i,j) = ssim(img2(indi,indj),template);            
            if mod(k,size(ssimval,2))==0; fprintf('%d/%d: %0.1f secs\n',k,K,toc); end; k=k+1;
        end
    end
end