clc
clear all
%% run all models

% structuralsim
incfg.dinamic_model   = 'structuralsim';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'deepgaze';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'structuralsim';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'flat';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'structuralsim';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'noisy';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'structuralsim';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'center';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

% cIBS
incfg.dinamic_model   = 'correlation';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'deepgaze';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'correlation';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'flat';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'correlation';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'noisy';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'correlation';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'center';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

% IBS
incfg.dinamic_model   = 'geisler';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'deepgaze';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'geisler';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'flat';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'geisler';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'noisy';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'geisler';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'center';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

% Greedy
incfg.dinamic_model   = 'greedy';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'deepgaze';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'greedy';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'flat';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'greedy';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'noisy';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)

incfg.dinamic_model   = 'greedy';
incfg.delta           = 32;
incfg.a               = 3;               
incfg.b               = 4;            
incfg.static_model    = 'center';
incfg.norm_cdf_tolerance = 0.0001;
incfg.parfor = 1;
main(incfg)