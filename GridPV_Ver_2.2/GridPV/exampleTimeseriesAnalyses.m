%% exampleTimeseriesAnalyses
% Timeseries analysis and plots monitor values from the simulation
%
%% Syntax
%  exampleTimeseriesAnalyses(basecaseFile,solarScenarioFiles);
%  exampleTimeseriesAnalyses();
%
%% Description
% Example function for timeseries analysis and monitor plotting for net 
% feeder power and switching components like LTC and capacitors.  Monitors
% must be setup in the basecaseFile circuit definition.  Place monitors in
% the desired locations, then use the same names in the code in this
% function.
%
%% Inputs
% * *|basecaseFile|* - optional input with the link to the OpenDSS file
% with the circuit.
% * *|solarScenarioFiles|* - optional input with a cell array of links to
% the OpenDSS files with the solar scenarios to run
%
%% Outputs
% * *none* - generates several figures and saves them
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
% exampleTimeseriesAnalyses('ExampleCircuit\master_ckt24.dss',{'ExampleCircuit\Ckt24_PV_Central_7_5.dss'})
%

function exampleTimeseriesAnalyses(basecaseFile,solarScenarioFiles)

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

    
    %% Run OpenDSS simulation for 1-week at 1-minute resolution
    DSSText.command = 'Set mode=duty number=10080  hour=0  h=60 sec=0';
    DSSText.Command = 'Set Controlmode=TIME';
    DSSText.command = 'solve';
    
    
    %% Feeder Power
    fileNameNoPath = DSSfilename(find(DSSfilename=='\',1,'last')+1:end-4);
    plotMonitor(DSSCircObj,'fdr_05410_Mon_PQ');
    ylabel('Power (kW,kVar)','FontSize',12,'FontWeight','bold')
    title([strrep(fileNameNoPath,'_',' '),' Net Feeder 05410 Load'],'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_Net_Power.fig'])
    
    
    %% Feeder Power Factor
    DSSText.Command = 'export mon fdr_05410_Mon_PQ';
    monitorFile = DSSText.Result;
    MyCSV = importdata(monitorFile);
    delete(monitorFile);
    Hour = MyCSV.data(:,1); Second = MyCSV.data(:,2);
    feederPower = MyCSV.data(:,[3,5,7]);
    feederReactivePower = MyCSV.data(:,[4,6,8]);
    
    figure;
    plot(Hour+Second/3600,abs(feederPower(:,1))./hypot(feederPower(:,1),feederReactivePower(:,1)),'LineWidth',2);
    grid on;
    set(gca,'FontSize',10,'FontWeight','bold')
    xlabel('Hour','FontSize',12,'FontWeight','bold')
    ylabel('Power Factor','FontSize',12,'FontWeight','bold')
    title([fileNameNoPath,sprintf(' Feeder Mean PF: %3.2f',mean(abs(feederPower(:,1))./hypot(feederPower(:,1),feederReactivePower(:,1))))],'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_PF.fig'])
    
    
    %% Substation Power
    DSSText.Command = 'export mon OtherFeeder_Mon_PQ';
    monitorFile = DSSText.Result;
    MyCSV = importdata(monitorFile);
    delete(monitorFile);
    otherFeederPower = MyCSV.data(:,[3,5,7]);
    otherFeederReactivePower = MyCSV.data(:,[4,6,8]);
    
    figure;
    plot(Hour+Second/3600,feederPower+otherFeederPower,'LineWidth',2); hold all;
    plot(Hour+Second/3600,feederReactivePower+otherFeederReactivePower,'LineWidth',2);
    legend('P1 (kW)','P2 (kW)','P3 (kW)','Q1 (kvar)','Q2 (kvar)','Q3 (kvar)');
    grid on;
    set(gca,'FontSize',10,'FontWeight','bold')
    xlabel('Hour','FontSize',12,'FontWeight','bold')
    ylabel('Power (kW,kVar)','FontSize',12,'FontWeight','bold')
    title([strrep(fileNameNoPath,'_',' '),' Net Substation Load'],'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_Net_Sub_Power.fig'])
    
    
    %% LTC
    DSSText.Command = 'export mon LTC';
    monitorFile = DSSText.Result;
    MyCSV = importdata(monitorFile);
    delete(monitorFile);
    Hour = MyCSV.data(:,1); Second = MyCSV.data(:,2);
    LTCtap = MyCSV.data(:,3);
    
    figure;
    plot(Hour+Second/3600,LTCtap,'LineWidth',2);
    set(gca,'FontSize',10,'FontWeight','bold')
    xlabel('Hour','FontSize',12,'FontWeight','bold')
    ylabel('LTC Tap','FontSize',12,'FontWeight','bold')
    title(sprintf('LTC Taps: %i Changes during Week',sum(diff(LTCtap)~=0)),'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_LTC.fig'])   
    
    
    %% Substation Voltage
    DSSText.Command = 'export mon subVI';
    monitorFile = DSSText.Result;
    MyCSV = importdata(monitorFile);
    delete(monitorFile);
    Hour = MyCSV.data(:,1); Second = MyCSV.data(:,2);
    subVoltages = MyCSV.data(:,3:2:7);
    
    figure;
    plot(Hour+Second/3600,subVoltages,'LineWidth',2);
    grid on;
    set(gca,'FontSize',10,'FontWeight','bold')
    xlabel('Hour','FontSize',12,'FontWeight','bold')
    ylabel('Voltage','FontSize',12,'FontWeight','bold')
    title([strrep(fileNameNoPath,'_',' '),' Substation Voltages'],'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_Sub_Voltage.fig'])
    
    
    %% Cap1 Voltage
    DSSText.Command = 'export mon Cap1VI';
    monitorFile = DSSText.Result;
    MyCSV = importdata(monitorFile);
    delete(monitorFile);
    Hour = MyCSV.data(:,1); Second = MyCSV.data(:,2);
    cap1voltages = MyCSV.data(:,3:2:size(MyCSV.data,2)/2);
    
    figure;
    plot(Hour+Second/3600,cap1voltages,'LineWidth',2);
    grid on;
    set(gca,'FontSize',10,'FontWeight','bold')
    xlabel('Hour','FontSize',12,'FontWeight','bold')
    ylabel('Voltage','FontSize',12,'FontWeight','bold')
    title([strrep(fileNameNoPath,'_',' '),' Cap1 Voltages'],'FontSize',12,'FontWeight','bold')
    saveas(gcf,[DSSfilename(1:end-4),'_Cap1_Voltage.fig'])
    
    
    
    
    %% Save data (in case you want to do future anlaysis, so that you don't have to run the simulation again)
    save([DSSfilename(1:end-4),'Analyses.mat'],'feederPower','feederReactivePower','LTCtap','subVoltages','cap1voltages');
    clear('feederPower','feederReactivePower','LTCtap','subVoltages','cap1voltages');

end





end
