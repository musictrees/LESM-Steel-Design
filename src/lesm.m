%% LESM - Linear Elements Structure Model
%
% This is the main script file of the LESM program.
%
% Run this script to launch the graphical interface
% and start using the program.
%
clc; clearvars; close all;
if ~isdeployed
    addpath(genpath(pwd));
end
launch();