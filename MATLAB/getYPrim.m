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
yy=zeros(8760,18);
%%
power_flow = tic;
for i = 1:8760
    DSSSolution.Solve();
    if (~DSSSolution.Converged)
        error('Solution Not Converged. Check Model for Convergence');
    end
    DSSObj.ActiveCircuit.SetActiveElement('Load.671');
    yy(i,:) = DSSObj.ActiveCircuit.ActiveCktElement.Yprim; % Get and save the Y prim of load 671, the variable is saved as a MAT file, check the odd column 
%     entries for each row as that correspond to the real power portion of the load.
    
end
t = toc (power_flow);
total_time = toc (start_time);


sprintf('Displaying Power Flow Solution Time: %f', t) 
sprintf('Displaying Total Time: %f', total_time)




