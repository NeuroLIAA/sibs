% grilla tomada de a_b_comparison de Seba
ab = [2 round(2/3,4); 2 round(4/3,4); 2 2; 2 round(8/3,4); 2 round(16/5,4); 2 4; 2 round(16/3,4);
      3 1; 3 2; 3 3; 3 4; 3 6; 3 8; 3 round(24/5,4);
      4 4; 4 8; 4 round(32/5,4); 4 round(4/3,4); 4 round(8/3,4); 4 round(16/3,4); 4 round(32/3,4); 
      5 5; 5 8; 5 10; 5 round(5/3,4); 5 round(10/3,4); 5 round(20/3,4); 5 round(40/3,4)];
  
n_grid = length(ab);
for k = 1:n_grid
    ab_k = ab(k,:);
    a = ab_k(1);
    b = ab_k(2);
    
    fprintf('\nSearch: %d  of %d\n', k, n_grid);
    fprintf('A: %d  \nB: %d \n', a, b);
    
    incfg.dinamic_model   = 'correlation';
    incfg.iniimg 	      = 1;
    incfg.delta           = 32;
    incfg.a               = a;
    incfg.b               = b;  
    incfg.static_model    = 'deepgaze';
    incfg.norm_cdf_tolerance = 0.0001;
    %incfg.norm_cdf_tolerance = 0;
    incfg.parfor = 1;
    
    %main(incfg)
    %clc
    
end