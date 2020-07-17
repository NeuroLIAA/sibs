clc
clear all

% % Noisy
% incfg.dinamic_model   = 'correlation';
% incfg.iniimg 	      = 1;
% incfg.delta           = 32;
% incfg.a               = 3;            % integers (?)   
% incfg.b               = 4;            % integers (?)
% incfg.static_model    = 'noisy';
% main(incfg)
% 
% incfg.dinamic_model   = 'geisler';
% incfg.iniimg 	      = 77;
% main(incfg)
% 
% incfg.dinamic_model   = 'greedy';
% incfg.iniimg 	      = 1;
% main(incfg)

% Flat
incfg.dinamic_model   = 'correlation';
incfg.iniimg 	      = 134;
incfg.delta           = 32;
incfg.a               = 3;            % integers (?)   
incfg.b               = 4;            % integers (?)
incfg.static_model    = 'flat';
incfg.norm_cdf_tolerance = 0.001;
main(incfg)
% 
% incfg.dinamic_model   = 'greedy';
% incfg.iniimg 	      = 1;
% % =======
% incfg.dinamic_model   = 'greedy';
% incfg.delta           = 32;
% incfg.a               = 3;            % integers (?)   
% incfg.b               = 4;            % integers (?)
% incfg.static_model    = 'icf';
% main(incfg)
% 
% incfg.dinamic_model   = 'geisler';
% incfg.delta           = 32;
% incfg.a               = 3;            % integers (?)   
% incfg.b               = 4;            % integers (?)
% incfg.static_model    = 'icf';
% main(incfg)
% 
% incfg.dinamic_model   = 'correlation';
% incfg.delta           = 32;
% incfg.a               = 3;            % integers (?)   
% incfg.b               = 4;            % integers (?)
% incfg.static_model    = 'icf';
% 
% main(incfg)
