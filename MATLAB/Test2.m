clc
clear
close all
start_time = tic;
%% Load the components related to OpenDSS

BaseParam=struct();
Baseparam.slack_voltage=1.02;
[DSSObj, DSSText, gridpvpath] = DSSStartup;
DSSSolution=DSSObj.ActiveCircuit.Solution;
DSSObj.ClearAll();

%     disp(DSSObj.ActiveCircuit)
% DSSText.command = 'Compile C:\feeders\MultiPhase\13BusQSTS\IEEE13Nodeckt.dss';
DSSText.command = 'Compile C:\feeders\MultiPhase\123BusQSTS\IEEE123Master.dss';
setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
DSSSolution.Number= 1; 
v= zeros(1,8760);
%%
power_flow = tic;
DSSSolution.Solve();
for i = 1:8760
    DSSSolution.Solve();
    if (~DSSSolution.Converged)
        error('Solution Not Converged. Check Model for Convergence');
    end
    bus = getBusInfo(DSSObj);
    voltage = [bus.phaseVoltagesPU];
    zero_points = find(voltage ==0);
    voltage (zero_points) = [];
    v (1,i)= min(voltage);

    

end
t = toc (power_flow);
total_time = toc (start_time);
plot(v)

sprintf('Displaying Power Flow Solution Time: %f', t)
sprintf('Displaying Total Time: %f', total_time)

A = randn(39,39);
b = randn (39,1);

tic 
c = A*b;
toc




% 
% The IEEE 34 bus has its minimum voltage to be 0.8259

%     DSSText.command = 'Compile C:\feeders\MultiPhase\34Bus\ieee34Mod1.dss';
%     DSSText.command = 'Compile C:\feeders\MultiPhase\37Bus\ieee37.dss';
%     DSSText.command = 'Compile C:\feeders\MultiPhase\123Bus\IEEE123Master.dss';