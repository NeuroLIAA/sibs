function main(incfg)
    % Load data
    addpath(genpath('utils/'));
    load('../matrix/initial_fixations.mat');
    load('../matrix/target_positions_filtered.mat');
    % Gaston edit
    
    % Parameter configuration
    if nargin<1
        cfg.static_model    = 'mlnet';
        cfg.dinamic_model   = 'correlation';
        cfg.delta           = 32;
        cfg.a               = 3;            % integers (?)   
        cfg.b               = 4;            % integers (?)
        cfg.iniimg          = 1;            % integers (?)
    else
        if ~isfield(incfg,'static_model') || isempty(incfg.static_model); 
            cfg.static_model = 'mlnet'; 
        else
            cfg.static_model = incfg.static_model;
        end
        if ~isfield(incfg,'dinamic_model') || isempty(incfg.dinamic_model); 
            cfg.dinamic_model = 'correlation'; 
        else
            cfg.dinamic_model = incfg.dinamic_model;
        end        
        if ~isfield(incfg,'delta') || isempty(incfg.delta); 
            cfg.delta = 32; 
        else
            cfg.delta = incfg.delta;
        end        
        if ~isfield(incfg,'a') || isempty(incfg.a);
            cfg.a = 3;
        else
            cfg.a = incfg.a;
        end
        if ~isfield(incfg,'b') || isempty(incfg.b); 
            cfg.b = 4;
        else
            cfg.b = incfg.b;
        end
        if ~isfield(incfg,'iniimg') || isempty(incfg.iniimg) 
            cfg.iniimg = 1;
        else
            cfg.iniimg = incfg.iniimg;
        end
        
        if ~isfield(incfg,'cache_path') || isempty(incfg.cache_path) 
            cfg.cache_path = '../cache/';
        else
            cfg.cache_path = incfg.cache_path;
        end
        
        if ~exist(cfg.cache_path, 'dir')
            disp('Creating cache folder...')
            mkdir(cfg.cache_path)
        end
        
        if ~isfield(incfg,'norm_cdf_tolerance') || isempty(incfg.norm_cdf_tolerance) 
            cfg.norm_cdf_tolerance = 0.001;       
        else
            cfg.norm_cdf_tolerance = incfg.norm_cdf_tolerance;
        end
        
        if cfg.norm_cdf_tolerance>0
            filename = [cfg.cache_path 'normcdf_table_' sprintf('%g',cfg.norm_cdf_tolerance) '.mat'];
            if exist(filename,'file')
                load(filename,'norm_cdf_table');               
            else
                norm_cdf_table = create_normcdf_lut(cfg.norm_cdf_tolerance);
                save(filename,'norm_cdf_table');
            end
            cfg.norm_cdf_table = norm_cdf_table;
        else
            cfg.norm_cdf_table = struct();
        end
        
        if ~isfield(incfg,'seed') || isempty(incfg.seed)
            cfg.seed = 1234;
        else
            cfg.seed = incfg.seed;
        end
            
        
        
    end
    
    cfg.img_quantity    = 134;
    cfg.nsaccades_thr   = 15; %10; 
    cfg.target_size     = [72 72];
    cfg.image_size      = [768 1024];
    cfg.out_models_path = ['../out_models/' cfg.static_model '/' cfg.dinamic_model '/a_' num2str(cfg.a) '_b_' num2str(cfg.b) '_tam_celda_' num2str(cfg.delta)];

    if ~exist(cfg.out_models_path,'dir') 
        mkdir(cfg.out_models_path);
        mkdir(sprintf('%s/cfg/',cfg.out_models_path));
        mkdir(sprintf('%s/posterior/',cfg.out_models_path));
        mkdir(sprintf('%s/probability_map/',cfg.out_models_path));
        mkdir(sprintf('%s/scanpath/',cfg.out_models_path));
    end
        
    img_time          = [];
    fprintf('\n\na: %d; b: %d; Delta: %d \n', cfg.a, cfg.b, cfg.delta);
    
    % Main loop. Run bayesian model for each image
    for imgnum = cfg.iniimg:cfg.img_quantity
        % Parameter configuration for each image
        fprintf('\nImage: %d  \n', imgnum);
        cfg.initial_fix   = initial_fixations(imgnum).initial_fix;
        cfg.imgnum        = imgnum;
        cfg.imgname       = initial_fixations(imgnum).image;    
        cfg.target_center = [target_positions(imgnum).matched_row target_positions(imgnum).matched_column] + cfg.target_size/2;
        rng(cfg.seed)
        
        % Run bayesian model
        tic;
        bayesianSaliencyModel(cfg); % This function saves ...
        tiempo = toc;
        fprintf('\ntiempo: %d  \n', tiempo);
        img_time = [img_time tiempo];
    end
    save([cfg.out_models_path , '/cfg/time.mat'], 'img_time');
end
