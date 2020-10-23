%% exampleVoltageAnalysis
% Example analysis of maximum and minimum feeder voltages through time
%
%% Syntax
%  exampleVoltageAnalysis(basecaseFile,solarScenarioFiles);
%  exampleVoltageAnalysis();
%
%% Description
% Example function for analysis of maximum and minimum feeder voltages
% through time. The simulation stops at each time step for MATLAB to
% process the state of the OpenDSS simulation
%
%% Inputs
% * *|basecaseFile|* - optional input with the link to the OpenDSS file
% with the circuit.
% * *|solarScenarioFiles|* - optional input with a cell array of links to
% the OpenDSS files with the solar scenarios to run
%
%% Outputs
% * *none* - generates a plot of maximum and minimum voltage through time
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
% Runs the basecase circuit and the distributed solar case
%%
% exampleVoltageAnalysis('ExampleCircuit\master_ckt24.dss',{'ExampleCircuit\Ckt24_PV_Distributed_7_5.dss'})
%

function exampleVoltageAnalysis(basecaseFile,solarScenarioFiles)

%% Check input arguments and get file paths from user if not in arguments
if nargin~=2
    [FileName,PathName] = uigetfile('*.dss','Basecase: Select the OpenDSS file with the circuit');
    basecaseFile = [PathName,FileName];
    [FileName,PathName] = uigetfile('*.dss','Solar: Select the OpenDSS files with the solar scenarios to run','MultiSelect','on');
    
    if iscell(FileName)
        for ii=1:length(FileName)
            solarScenarioFiles{ii} = [PathName,FileName{ii}];
        end
    else
        solarScenarioFiles = {[PathName,FileName]};
    end
end


%% initiate COM interface (only need to do once when you open MATLAB)
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
% Keep the OpenDSS progress window from popping up during solve
DSSCircObj.Allowforms = false;


%% Run through each scenario to do analysis
caseNames = [basecaseFile,solarScenarioFiles];

for jj=1:length(caseNames)
    
    %% Call OpenDSS compile
    DSSfilename = caseNames{jj};
    location = cd;
    DSSText.command = sprintf('Compile (%s)',caseNames{1}); %run basecase first
    DSSText.command = 'solve';
    cd(location);
    DSSText.command = sprintf('Compile (%s)',DSSfilename); %add solar scenario
    DSSText.command = 'solve';
    cd(location);
    
    
    %% Run simulations every 1-minute and find max/min voltages
    simulationResolution = 60; %in seconds
    simulationSteps = 24*60*7;
    
    DSSText.Command = sprintf('Set mode=duty number=1  hour=0  h=%i sec=0',simulationResolution);
    DSSText.Command = 'Set Controlmode=TIME';
    
    
    minVoltages = zeros(simulationSteps,1);
    minDistances = zeros(simulationSteps,1);
    minDistanceNames = cell(simulationSteps,1);
    maxVoltages = zeros(simulationSteps,1);
    maxDistances = zeros(simulationSteps,1);
    maxDistanceNames = cell(simulationSteps,1);
    
    absoluteMinVoltage = 2;
    absoluteMaxVoltage = 0;
    for ii=1:simulationSteps
        
        %% solve again
        DSSText.Command = 'solve';
        
        %% Get Voltages at buses
        voltages = DSSCircuit.AllBusVmagPU; %PU voltages at all buses per phase
        voltagesInVolts = DSSCircuit.AllBusVmag; %voltages at all buses per phase
        voltageBusNames = DSSCircuit.AllNodeNames; %names of each bus per phase
        distances = DSSCircuit.AllNodeDistances; %distances of each bus per phase
        
        %% remove transmission line voltages (otherwise highest PU voltage will always be the fixed transmission line voltage)
        voltages = voltages(voltagesInVolts<=15000); 
        voltageBusNames = voltageBusNames(voltagesInVolts<=15000);
        distances = distances(voltagesInVolts<=15000);
        
        %% remove zeros
        condition = voltages>0.5; %everything should be above 0.5 pu
        voltages = voltages(condition); 
        voltageBusNames = voltageBusNames(condition);
        distances = distances(condition);
        
        %% find max and min voltage, and store in array
        [minVoltages(ii) minIndex] = min(voltages);
        minDistances(ii) = distances(minIndex); %record distance where min occured
        minDistanceNames{ii} = voltageBusNames{minIndex}; %record name of bus where max occured
        if minVoltages(ii)<absoluteMinVoltage % store max/min voltage location
            minLocation = voltageBusNames(minIndex);
            absoluteMinVoltage = minVoltages(ii);
        end
        
        [maxVoltages(ii) maxIndex] = max(voltages);
        maxDistances(ii) = distances(maxIndex); %record distance where max occured
        maxDistanceNames{ii} = voltageBusNames{maxIndex}; %record name of bus where max occured
        if maxVoltages(ii)>absoluteMaxVoltage % store max/min voltage location
            maxLocation = voltageBusNames(maxIndex);
            absoluteMaxVoltage=maxVoltages(ii);
        end
        
    
    end
        
    %% Plotting
    fileNameNoPath = DSSfilename(find(DSSfilename=='\',1,'last')+1:end-4);
    figure;
    plot(simulationResolution/3600:simulationResolution/3600:simulationSteps*simulationResolution/3600,maxVoltages*120,'LineWidth',2)
    hold all;
    plot(simulationResolution/3600:simulationResolution/3600:simulationSteps*simulationResolution/3600,minVoltages*120,'LineWidth',2)
    set(gca,'FontSize',10,'FontWeight','bold')
    legend('Max Feeder Voltage','Min Feeder Voltage')
    xlabel('Hour','FontSize',12,'FontWeight','bold')
    ylabel('Voltage (120 V Base)','FontSize',12,'FontWeight','bold')
    title([fileNameNoPath,' Voltages'],'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_Voltages.fig'])
    grid on;
    
    figure;
    plot(simulationResolution/3600:simulationResolution/3600:simulationSteps*simulationResolution/3600,maxDistances,'LineWidth',2)
    hold all;
    plot(simulationResolution/3600:simulationResolution/3600:simulationSteps*simulationResolution/3600,minDistances,'LineWidth',2)
    set(gca,'FontSize',10,'FontWeight','bold')
    legend('Distance to Max Voltage Bus','Distance to Min Voltage Bus')
    xlabel('Hour','FontSize',12,'FontWeight','bold')
    ylabel('Distance','FontSize',12,'FontWeight','bold')
    title([fileNameNoPath,' Voltages'],'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_Voltage_Distances.fig'])
    grid on;
    
    %% Save results to MAT file (in case you want to do future anlaysis, so that you don't have to run the simulation again)
    save([DSSfilename(1:end-4),'Voltages.mat'],'minVoltages','maxVoltages','minDistances','maxDistances','absoluteMaxVoltage','absoluteMinVoltage','maxLocation','minLocation','maxDistanceNames','minDistanceNames');
    clear('minVoltages','maxVoltages');

end

