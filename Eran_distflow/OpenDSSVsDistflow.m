clc;
clear all
close all;

%% Load the data .mat file
model = load('IEEE_123.mat');
%% Lossy Case
opt.alpha=0.5;
opt.alpha_method=1;
% Run the distflow model
[BusResults,~]=distflow_multi(model.Bus,model.Branch,opt);


%% Figure, plotting the lateral voltage
figure
lateral_plot(BusResults,model.Branch,model.OpenDSSResults)
figure
lateral_plot(BusResults,model.Branch,model.OpenDSSResults,'diffplot')