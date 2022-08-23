clear all

incfg.dinamic_model   = 'structuralsim';
incfg.iniimg 	      = 1;
incfg.endimg          = 1;
incfg.delta           = 32;
incfg.a               = 3;            % integers (?)   
incfg.b               = 4;            % integers (?)
incfg.static_model    = 'noisy';
incfg.norm_cdf_tolerance = 0.0001;
% incfg.norm_cdf_tolerance = 0;
incfg.parfor = 1;
main(incfg)

% incfg.dinamic_model   = 'structuralsim';
% incfg.iniimg 	      = 1;
% incfg.endimg          = 134;
% incfg.delta           = 32;
% incfg.a               = 3;            % integers (?)   
% incfg.b               = 4;            % integers (?)
% incfg.static_model    = 'flat';
% incfg.norm_cdf_tolerance = 0.0001;
% % incfg.norm_cdf_tolerance = 0;
% incfg.parfor = 1;
% main(incfg)