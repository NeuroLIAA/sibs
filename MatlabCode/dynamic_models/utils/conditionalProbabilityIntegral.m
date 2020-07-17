function result = conditionalProbabilityIntegral(ix, iy, kx, ky, T, p, visibility_map, alpha, norm_cdf_table)
%    m = (-2*log(p(:,:,T) / p(ix,iy,T)) + visibility_map(:,:,kx,ky).^2 + visibility_map(ix,iy,kx,ky).^2) ./ (2 * visibility_map(:,:,kx,ky)); % MAL
%    b =  visibility_map(ix,iy,kx,ky) ./  visibility_map(:,:,kx,ky); % MAL

    b = (-2*log(p(:,:,T) / p(ix,iy,T)) + visibility_map(:,:,kx,ky).^2 + visibility_map(ix,iy,kx,ky).^2) ./ (2 * visibility_map(:,:,kx,ky)); % BIEN
    m =  visibility_map(ix,iy,kx,ky) ./  visibility_map(:,:,kx,ky); % BIEN

    % we ensure that the product is only for i != j (normcdf(1000000) = 1)
    m(ix, iy) = 0;
    b(ix, iy) = 1000000;

    % check limits to integrate (normcdf(-20) = 0 and so will be the product)
    minw = max(max((-20 - b(m>0)) ./ m(m>0)), -20);
    if isempty(minw)
        minw = -20;
    end
    maxw = min(min((-20 - b(m<0)) ./ m(m<0)), 20);
    if isempty(maxw)
        maxw = 20;
    end
        
    if minw >= maxw
%        [ix iy minw maxw] 
        result = 0;
        return;
    end
    
    num = 50;
    wRange = linspace(minw, maxw, num);
    
    tmp = m(:) * wRange + repmat(b(:), 1, length(wRange));
    tmp(isnan(tmp)) = 1;

%    phiW = exp(-0.5 .* wRange .* wRange / sqrt(2 * pi)); % MAL
    phiW = exp(-0.5 .* wRange .* wRange) / sqrt(2 * pi); % BIEN
    if isstruct(norm_cdf_table)
        normcdf_tmp = interp1(norm_cdf_table.x,norm_cdf_table.y,tmp,'nearest','extrap');
    else
        normcdf_tmp = normcdf(tmp);
    end 
        point = phiW .* (prod(alpha * normcdf_tmp, 1) / alpha);

    result = trapz(wRange, point);
end
