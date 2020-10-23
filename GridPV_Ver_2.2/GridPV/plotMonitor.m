%% plotMonitor
% Plots a monitor from the simulation
%
%% Syntax
%  plotMonitor(DSSCircObj,monitorName);
%
%% Description
% Function to plot the simulation results saved in an OpenDSS monitor
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|monitorName|* - string with the name of the OpenDSS monitor
%
%% Outputs
% * *none* - a figure is displayed with the plot
%
%% Copyright 2014
% Georgia Tech Research Corporation, Atlanta, Georgia 30332
% Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
% See the license agreement for full terms and conditions.
%
% Please acknowledge any contributions of the GridPV Toolbox by citing:
% M. J. Reno and K. Coogan, "Grid Integrated Distributed PV (GridPV) Version 2," Sandia National Laboratories SAND2013-20141, 2014.
%
%% Example
% Example of a feeder power monitor plot
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'Set mode=duty number=8760  hour=0  h=1h sec=0';
% DSSText.command = 'Set controlmode = time';
% DSSText.command = 'solve';
% plotMonitor(DSSCircObj,'fdr_05410_Mon_PQ')
%

function plotMonitor(DSSCircObj,monitorName)
%% Parse Inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addRequired('monitorName', @ischar);

p.parse(DSSCircObj, monitorName); %parse inputs

% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
% Define the text interface
DSSText = DSSCircObj.Text;

DSSText.Command = sprintf('export mon %s',monitorName);
    while isempty(DSSText.Result) end %wait for export to finish
monitorFile = DSSText.Result;
MyCSV = importdata(monitorFile);
delete(monitorFile);
Hour = MyCSV.data(:,1); 
Second = MyCSV.data(:,2);
labels = MyCSV.colheaders(3:end);
for ii=1:length(labels)
    data(:,ii) = MyCSV.data(:,2+ii);
end

plot(Hour+Second/3600,data,'LineWidth',2);
grid on;
set(gca,'FontSize',10,'FontWeight','bold')
xlabel('Hour','FontSize',12,'FontWeight','bold')
legend(labels)


end