clc
clear
close all
current_directory=pwd;

savemodel=0;


%% Load the components related to OpenDSS
[DSSObj, DSSText, gridpvpath] = DSSStartup;
% DSSCircuit=DSSObj.ActiveCircuit;
DSSSolution=DSSObj.ActiveCircuit.Solution;
DSSText.command = 'Clear';
% DSSText.command = 'Compile C:\feeders\MultiPhase\34Bus\ieee34Mod1.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\34Bus\ieee34Mod1Test.dss';
DSSText.command = 'Compile C:\feeders\MultiPhase\13Bus\IEEE13Nodeckt.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\37Bus\ieee37.dss';
% DSSText.command = 'Compile C:\feeders\MultiPhase\123Bus\IEEE123Master';

%% Define the base parameters
BaseParam=struct();
Baseparam.slack_voltage=1.02;
BaseParam.VLL=DSSObj.ActiveCircuit.Vsources.BasekV;
Baseparam.sbase=10;
Baseparam.zbase=(BaseParam.VLL)^2/Baseparam.sbase;
Baseparam.delta_index=[0 1 3;1 0 2;3 2 0];
Baseparam.delta2Yconvmatrix=[1 0 -1; -1 1 0; 0 -1 1];

%% Solving the Power Flow
setSourceInfo(DSSObj,{'source'},'pu',Baseparam.slack_voltage);
DSSSolution.Solve();
if (~DSSSolution.Converged)
    error('Solution Not Converged. Check Model for Convergence');
else
    disp('Model Converged.')
end

%% Additional Changes Made
% Loads=getLoadInfo(DSSObj);
% for i = 1:length(Loads)
%     if (Loads(i).isDelta )
%         setLoadInfo(DSSObj,{Loads(i).name},'model',2);
%     end
% end
% DSSSolution.Solve();
% Loads=getLoadInfo(DSSObj);
% for i = 1:length(Loads)
%     if (Loads(i).isDelta )
%         fprintf('%s: %s %d\n',Loads(i).name,Loads(i).model,Loads(i).isDelta);
%     end
% end
% if (~DSSSolution.Converged)
%     error('Solution Not Converged. Check Model for Convergence');
% else
%     disp('Model Converged.')
% end

% Loads=getLoadInfo(DSSObj);
% setLoadInfo(DSSObj,{Loads.name}, 'Model', 2*ones(1,length(Loads)));
% DSSSolution.Solve();
%% OpenDSS simulation changes directory
cd (current_directory);

%% Converting the OpenDSS Model to Distflow Model
[Bus,Branch,BusVoltage]=create_distflow_model(DSSObj,Baseparam);

%% Lossy Case
opt.alpha=0.5;
opt.alpha_method=1;
[BusResults,~]=distflow_multi(Bus,Branch,opt);
% OpenDSSResults=dssvm2bus(Bus,BusVoltage);

%% Figure, plotting the lateral voltage
figure
lateral_plot(BusResults,Branch,Bus)
% figure
% lateral_plot(BusResults,Branch,OpenDSSResults,'diffplot')

% %%
if (savemodel==1)
    save('IEEE_34.mat','Bus','Branch','OpenDSSResults');
end

model = load('IEEE_13.mat');





