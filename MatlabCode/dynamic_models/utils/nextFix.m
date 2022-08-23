function [idx_x, idx_y, detectability_map] = nextFix(cfg, T, p, visibility_map)
    % nextFix: Compute the next fix
    
    % Input:
    % cfg            = configuration
    % T              = fix number
    % p              = probability map
    % visibility_map = visibility map

    if strcmp(cfg.dinamic_model, 'greedy') 
        [idx_x, idx_y] = greedy(p,T);
    elseif strcmp(cfg.dinamic_model, 'geisler') || strcmp(cfg.dinamic_model, 'correlation') ...
            || strcmp(cfg.dinamic_model, 'structuralsim')
        [idx_x, idx_y, detectability_map] = bayesian_model(cfg, T, p, visibility_map);
    else
        fprintf('Error: Wrong dinamic model \n')
    end
end
