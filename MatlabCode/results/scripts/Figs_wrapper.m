% Paper results wrapper
% This script will execute all the scripts needed to compute the figures
% and the tables presented in the paper

clear all

%% Run Figure 1: 
disp('Figure 1')
Fig1_paradigm;

%% Run Figure 2: 
disp('Figure 2')
Fig2_saliency;

%% Run Figure auxiliar calculations: 
disp('Figures auxiliar calculations')
Figs_calculations;

%% Run Figure 3: 
disp('Figure 3')
Fig3_performance;

%% Run Figure 4: 
disp('Figure 4')
Fig4_multimatch;

%% Run Table: 
disp('Table')
Table_metrics;
