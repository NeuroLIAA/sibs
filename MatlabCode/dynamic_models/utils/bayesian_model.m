function [idx_x, idx_y] = bayesian_model(cfg, T, p, visibility_map)
%     global dame1 dame2
    alpha = 1;
    if cfg.parfor
        accum2 = zeros(cfg.size_prior(1)* cfg.size_prior(2), cfg.nsaccades_thr);
        elements = {1:cfg.size_prior(1), 1:cfg.size_prior(2)};
        combs = combinations(elements);    
    
        parfor c = 1:length(combs)
            kx = combs(c,1);ky = combs(c,2);
            tmp = zeros(cfg.size_prior(1), cfg.size_prior(2));
            for ix = 1:cfg.size_prior(1)
                for iy = 1:cfg.size_prior(2)                    
                    tmp(ix, iy) = conditionalProbabilityIntegral(ix, iy, kx, ky, p(:,:,T), visibility_map, alpha, cfg.norm_cdf_table);
                end
            end
           accum2(c,T) = nansum(nansum(p(:,:,T) .* tmp));
        end
        accum = reshape(accum2, cfg.size_prior(1),cfg.size_prior(2),cfg.nsaccades_thr);
    else
        accum = zeros(cfg.size_prior(1), cfg.size_prior(2), cfg.nsaccades_thr);
        tmp = zeros(cfg.size_prior(1), cfg.size_prior(2));

        for kx = 1:cfg.size_prior(1)
            for ky = 1:cfg.size_prior(2)
                for ix = 1:cfg.size_prior(1)
                    for iy = 1:cfg.size_prior(2)                    
                        tmp(ix, iy) = conditionalProbabilityIntegral(ix, iy, kx, ky, p(:,:,T), visibility_map, alpha, cfg.norm_cdf_table);
                    end
                end
                accum(kx,ky,T) = nansum(nansum(p(:,:,T) .* tmp));
            end
        end
        
    end

    
%     dame1 = [dame1 reshape(p(:,:,T),1,[])];
%     dame2 = [dame2 reshape(tmp,1,[])];
    [~, idx_y] = max(max(accum(:,:,T)));
    [~, idx_x] = max(accum(:,idx_y,T));
end
