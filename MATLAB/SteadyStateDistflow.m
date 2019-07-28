% This code is to do the steady state analysis part 

clc
clear
close all
%  Start the timer
% start_time = tic;

% case_file = 'IEEE_13.mat';
% case_file = 'IEEE_34.mat';
% case_file = 'IEEE_37.mat';
case_file = 'IEEE_123.mat';

load(case_file);
total_simulation = 200; % The number of times i am willing to run the same power flow but initialized everytime, this can be arbitary any number
opt = struct('number_iteration', total_simulation); % This variable ensures that i am not rebuilding the model
% rather running the simulation in a loop
[rbus, rbranch, total_time_power_flow] = distflow_multi(Bus, Branch, opt);
disp(total_time_power_flow)
% For IEEE 123 bus test case the displayed time is 0.0015 seconds


%%
tic 
clear 
case_file = 'IEEE_123.mat';
total_simulation = 200;
load(case_file);
time_array=zeros(1,total_simulation);
opt = struct('number_iteration', 1);
% In this cae i am rebuilding the model every time
for i = 1:total_simulation
    [rbus, rbranch, total_time_power_flow] = distflow_multi(Bus, Branch, opt);
end

toc
% For the IEEE 123 bus test case it took 92.58 seconds