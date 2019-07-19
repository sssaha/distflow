% This code is to do the steady state analysis part 

clc
clear
close all
%  Start the timer
start_time = tic;
%% Load the components related to OpenDSS
[DSSObj, DSSText, gridpvpath] = DSSStartup;
DSSSolution=DSSObj.ActiveCircuit.Solution;
DSSText.command = 'Clear';
BaseParam=struct();
Baseparam.slack_voltage=1.02;
total_simulation = 200; % The number of times i am willing to run the same power flow but initialized everytime, this can be arbitary any number
time_array=zeros(1,total_simulation); %variable to save all the tic toc values 
%%
for i = 1:total_simulation % I am running 200 steps arbitarily, there is no special reason
    clear DSSObj DSSText gridpvpath DSSSolution % clearing up the variables at every iteration to make sure all the variables initialzed every time 
    [DSSObj, DSSText, gridpvpath] = DSSStartup; 
    DSSSolution=DSSObj.ActiveCircuit.Solution;
    DSSObj.ClearAll(); % clearing the memory completely
%   Choose the feeder we want to test and compile the model

%     DSSText.command = 'Compile C:\feeders\MultiPhase\13Bus\IEEE13Nodeckt.dss';
%     DSSText.command = 'Compile C:\feeders\MultiPhase\34Bus\ieee34Mod1.dss';
%     DSSText.command = 'Compile C:\feeders\MultiPhase\37Bus\ieee37.dss';
    DSSText.command = 'Compile C:\feeders\MultiPhase\123Bus\IEEE123Master.dss';

%     Set the slack bus voltage
    setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
%%
%     Start the power flow time
    power_flow = tic;
    
    DSSSolution.Solve(); % solve the power flow
    if (~DSSSolution.Converged) % check for convergence
        error('Solution Not Converged. Check Model for Convergence');
    end
    t = toc (power_flow); % retrieve the power flow ending time
%% 
%     The following code is to check whether the minimum voltage that i get
%     at each step is the same if i run the model from OpenDSS directly, so
%     far the result has matched for all the cases, but is kept out of the
%     time loop to ensure that we only record the power flow solving time.
    bus = getBusInfo(DSSObj);
    voltage = [bus.phaseVoltagesPU];
    zero_points = find(voltage ==0);
    voltage (zero_points) = [];
    disp(min(voltage))
    %%
    time_array(i) = t;
end
% The total time to run the full code recorded here, it also includes the time for the printing 
total_time = toc (start_time);

sprintf('Displaying Solution Time: %f', sum(time_array))
sprintf('Displaying Total Time: %f', total_time)


%% Time reported to run 200 steps so far on my Desktop running MATLAB 2018

% IEEE 13 Bus 0.316517      3.64594
% IEEE 34 Bus 0.390899      5.99447
% IEEE 37 bus 0.455304      6.461339
% IEEE 123 Bus 0.599445     15.5120


%% matrix multiplication testing for the 13 bus system, considering all the buses are 3 bus, we should have 39 entries
bus_size = 13;
A= rand (bus_size*3,bus_size*3);
x = rand (bus_size*3,1);
B= rand (bus_size*3,bus_size*3);
y = rand (bus_size*3,1);
tic
for i = 1:total_simulation
   c =  A*x + B*y;
end
toc
% Time reported for 13 bus case: 0.001014 second
% Time reported for 123 bus case: 0.017369 second

