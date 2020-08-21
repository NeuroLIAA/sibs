% clc
% clear all
% profile on
% 
% % Flat
% incfg.dinamic_model   = 'correlation';
% incfg.iniimg 	      = 134;
% incfg.delta           = 32;
% incfg.a               = 3;            % integers (?)   
% incfg.b               = 4;            % integers (?)
% incfg.static_model    = 'flat';
% incfg.norm_cdf_tolerance = 0;
% incfg.parfor = 0;
% main(incfg)
% 
% p = profile('info');
% save profiles/basic p
% 
% clc
% clear all
% profile on
% 
% % Flat
% incfg.dinamic_model   = 'correlation';
% incfg.iniimg 	      = 134;
% incfg.delta           = 32;
% incfg.a               = 3;            % integers (?)   
% incfg.b               = 4;            % integers (?)
% incfg.static_model    = 'flat';
% incfg.norm_cdf_tolerance = 0.0001;
% % incfg.norm_cdf_tolerance = 0;
% incfg.parfor = 0;
% main(incfg)
% 
% p = profile('info');
% save profiles/interp p


clc
clear all

parpool

profile on

% Flat
incfg.dinamic_model   = 'correlation';
incfg.iniimg 	      = 134;
incfg.delta           = 32;
incfg.a               = 3;            % integers (?)   
incfg.b               = 4;            % integers (?)
incfg.static_model    = 'flat';
% incfg.norm_cdf_tolerance = 0.0001;
incfg.norm_cdf_tolerance = 0;
incfg.parfor = 1;
main(incfg)

p = profile('info');
save profiles/parfor p


clc
clear all
profile on

% Flat
incfg.dinamic_model   = 'correlation';
incfg.iniimg 	      = 134;
incfg.delta           = 32;
incfg.a               = 3;            % integers (?)   
incfg.b               = 4;            % integers (?)
incfg.static_model    = 'flat';
incfg.norm_cdf_tolerance = 0.0001;
% incfg.norm_cdf_tolerance = 0;
incfg.parfor = 1;
main(incfg)

p = profile('info');
save profiles/interp_parfor p


