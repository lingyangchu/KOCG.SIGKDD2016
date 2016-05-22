%**************************************************************************
%********  Finding Gangs in War from Signed Networks, KDD2016  ************
%********     Author: Lingyang Chu, Simon Fraser University    ************
%********            Email:  chulingyang@hotmail.com           ************
%********                 All Rights Preserved                 ************
%**************************************************************************
%% This is a simple demo to demonstrate how KOCG works.
%% code as follows.
% clear env
clear;

% load DemoData.mat
load('./DATA/DEMO_DATA/DemoData.mat');

% Run enumKOCG.m
cd KOCG
[X_enumKOCG_cell, time] = enumKOCG(A, B, NG, alpha, beta)
cd ..
