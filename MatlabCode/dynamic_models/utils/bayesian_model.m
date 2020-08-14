function [idx_x, idx_y] = bayesian_model(cfg, T, p, visibility_map)
    tmp = zeros(cfg.size_prior(1), cfg.size_prior(2));
    accum = zeros(cfg.size_prior(1), cfg.size_prior(2), cfg.nsaccades_thr);
    alpha = 1;
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
    [~, idx_y] = max(max(accum(:,:,T)));
    [~, idx_x] = max(accum(:,idx_y,T));
end
