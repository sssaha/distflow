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
DSSText.command = 'Compile C:\feeders\MultiPhase\13BusQSTS\IEEE13Nodeckt.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\34Bus\ieee34Mod1.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\37Bus\ieee37.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\123BusQSTS\IEEE123Master.dss';
setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
% The following number is default to 8760, but if we want to find out the
% minimum voltage at each time step this has to set as 1
DSSSolution.Number= 1; 
% Variable to save the minimum per unit of all the buses
v= zeros(1,8760);
yy=zeros(8760,18);
%%
power_flow = tic;
for i = 1:8760
    DSSSolution.Solve();
    if (~DSSSolution.Converged)
        error('Solution Not Converged. Check Model for Convergence');
    end
% Use of the library puts 0 to absent phases to the buses, hence i had to
% filter out the buses with 0 pu voltage
%     bus = getBusInfo(DSSObj);
%     voltage = [bus.phaseVoltagesPU];
%     zero_points = find(voltage ==0);
%     voltage (zero_points) = [];
%     v (1,i)= min(voltage);
    DSSObj.ActiveCircuit.SetActiveElement('Load.671');
    yy(i,:) = DSSObj.ActiveCircuit.ActiveCktElement.Yprim;
    
end
t = toc (power_flow);
total_time = toc (start_time);
plot(v) % plotting the minimum voltage 

sprintf('Displaying Power Flow Solution Time: %f', t) % The output value is 710.418 
sprintf('Displaying Total Time: %f', total_time)% The output value is 712.945


% %%
% DSSText.command = 'Compile C:\feeders\MultiPhase\13BusQSTS\IEEE13Nodeckt.dss';
% setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
% qsts_start = tic;
% DSSSolution.Number= 8760; 
% DSSSolution.Solve(); % because the solution number is set as 8760, we no longer can retrieve values at each time step using the COM interface, and hence calcualting the minimum voltage at each time is not possible
% t = toc(qsts_start);
% sprintf('Displaying Yearly Solution Time: %f', t)

