clc
clear
close all
%% Load the components related to OpenDSS
[DSSObj, DSSText, gridpvpath] = DSSStartup;
DSSSolution=DSSObj.ActiveCircuit.Solution;
DSSText.command = 'Clear';
BaseParam=struct();
Baseparam.slack_voltage=1.02;
time_array=zeros(1,8760);
%%
% DSSText.command = 'Redirect C:\feeders\MultiPhase\13Bus\IEEE13Nodeckt.dss';
% setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
% disp('Simulating 13 bus')
% tic
for i = 1:8760
    clear DSSObj DSSText gridpvpath DSSSolution
    [DSSObj, DSSText, gridpvpath] = DSSStartup;
    DSSSolution=DSSObj.ActiveCircuit.Solution;
%     DSSObj.ClearAll();
%     disp(DSSObj.ActiveCircuit)
%     DSSText.command = 'Compile C:\feeders\MultiPhase\13Bus\IEEE13Nodeckt.dss';
%     DSSText.command = 'Compile C:\feeders\MultiPhase\34Bus\ieee34Mod1.dss';
    DSSText.command = 'Compile C:\feeders\MultiPhase\37Bus\ieee37.dss';
%     DSSText.command = 'Compile C:\feeders\MultiPhase\123Bus\IEEE123Master.dss';
    setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
    tic 
    DSSSolution.Solve();
    if (~DSSSolution.Converged)
        error('Solution Not Converged. Check Model for Convergence');
    end
    t = toc;
    time_array(i) = t;
end


sprintf('Displaying Solution Time: %f', sum(time_array))

