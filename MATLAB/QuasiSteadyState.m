% This code is for analyzing the QSTS aprt 
clc
clear
close all
% Start the main time
start_time = tic;
%% Load the components related to OpenDSS

BaseParam=struct();
Baseparam.slack_voltage=1.02;
[DSSObj, DSSText, gridpvpath] = DSSStartup;
DSSSolution=DSSObj.ActiveCircuit.Solution;
DSSObj.ClearAll();

%  Choose the circuit i want to compile
% DSSText.command = 'Compile C:\feeders\MultiPhase\13BusQSTS\IEEE13Nodeckt.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\34Bus\ieee34Mod1.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\37Bus\ieee37.dss';
DSSText.command = 'Compile C:\feeders\MultiPhase\123BusQSTS\IEEE123Master.dss';
setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
% The following number is default to 8760, but if we want to find out the
% minimum voltage at each time step this has to set as 1
DSSSolution.Number= 1; 
% Variable to save the minimum per unit of all the buses
v= zeros(1,8760);
%%
power_flow = tic;
for i = 1:8760
    DSSSolution.Solve();
    if (~DSSSolution.Converged)
        error('Solution Not Converged. Check Model for Convergence');
    end
% Use of the library puts 0 to absent phases to the buses, hence i had to
% filter out the buses with 0 pu voltage
    bus = getBusInfo(DSSObj);
    voltage = [bus.phaseVoltagesPU];
    zero_points = find(voltage ==0);
    voltage (zero_points) = [];
    v (1,i)= min(voltage);
end
t = toc (power_flow);
total_time = toc (start_time);
plot(v) % plotting the minimum voltage 

sprintf('Displaying Power Flow Solution Time: %f', t) % The output value is 710.418 
sprintf('Displaying Total Time: %f', total_time)% The output value is 712.945



%%
% DSSSolution.Solve();


%% Discussion
% Even though playing with many load profiles the power flow part remains
% fast, can we focus on application specific problems rather than focusing
% on load profile specific problems? We probably want to say using the
% matrix vector multiplication speeds up the overall process, rather than
% just looking into the power flow. If you are absolutely sure that we want
% to vary and try thousands of combinations of load profile files, i will
% have to create an automatic script to do that. Shifting load profiles is
% not a preferable thing, because then it becomes hard to justify why we
% did certain things. 


% If we are not interested in getting the minimum voltage at each time
% step, the code will look like this: 
DSSText.command = 'Compile C:\feeders\MultiPhase\123BusQSTS\IEEE123Master.dss';
setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
qsts_start = tic;
DSSSolution.Number= 8760; 
DSSSolution.Solve(); % because the solution number is set as 8760, we no longer can retrieve values at each time step using the COM interface, and hence calcualting the minimum voltage at each time is not possible
t = toc(qsts_start);
sprintf('Displaying Yearly Solution Time: %f', t) % The output value is 2.354529
% To show that it actually ran the 8760 steps i added a monitor assigned to
% a bus which recorded voltage, and then plotted that using the internal
% plot tool of OpenDSS

% As i have shown you two applications, please pick up one, and let me know
% whichever one you want to do. 

% The reason i am talking about application is, there is greater use case
% of using OpenDSS or Gridlabd with external applications, rather than
% running any analysis from the software itself, 

