% This code is to do the steady state analysis part 

clc
clear
close all
%  Start the timer
% start_time = tic;


load('IEEE_123.mat');
total_simulation = 200; % The number of times i am willing to run the same power flow but initialized everytime, this can be arbitary any number
time_array=zeros(1,total_simulation);
for i = 1:total_simulation
    [rbus, rbranch, time_array(1,i)] = distflow_multi(Bus, Branch);
end

disp(sum(time_array))
