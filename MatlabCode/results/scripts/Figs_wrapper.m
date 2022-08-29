% Paper results wrapper
% This script will execute all the scripts needed to compute the figures
% and the tables presented in the paper

clear all
close all

%% Run Figure 1: 
disp('To compute Figure 1 run script manually')
Fig1_paradigm;

%% Run Figure 2:
disp('To compute Figure 2 run script manually')
%Fig2_saliency;

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
