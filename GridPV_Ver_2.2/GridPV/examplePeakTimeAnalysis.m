%% examplePeakTimeAnalysis
% Runs simulation during peak penetration time and generates plots
%
%% Syntax
%  examplePeakTimeAnalysis(basecaseFile,solarScenarioFiles);
%  examplePeakTimeAnalysis();
%
%% Description
% Function to calculate when the max penetration (PV output / load) time
% occurs.  A snapshot analysis is performed at this peak time, with both a
% voltage contour plot and voltage profile plot being generated.
%
%% Inputs
% * *|basecaseFile|* - optional input with the link to the OpenDSS file
% with the circuit.
% * *|solarScenarioFiles|* - optional input with a cell array of links to
% the OpenDSS files with the solar scenarios to run
%
%% Outputs
% * *none* - generates 2 figures for each analysis scenario and saves them
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
% examplePeakTimeAnalysis('ExampleCircuit\master_ckt24.dss',{'ExampleCircuit\Ckt24_PV_Distributed_7_5.dss'})
%

function examplePeakTimeAnalysis(basecaseFile,solarScenarioFiles)

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


%% Find Peak time in order to do anlaysis
maxTimeIndex = findMaxPenetrationTime();


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
    
    
    %% Run the simulation in static mode for the peak time
    DSSText.command = sprintf('Set mode=duty number=1  hour=%i  h=1.0 sec=%i',floor((maxTimeIndex)/3600),round(mod(maxTimeIndex,3600)));
    DSSText.Command = 'Set Controlmode=Static'; %take control actions immediately without delays
    DSSText.command = 'solve';
    
    
    %% Plot voltage contour
    fileNameNoPath = DSSfilename(find(DSSfilename=='\',1,'last')+1:end-4);
    figure;
    plotCircuitLines(DSSCircObj,'Coloring','voltage120')
    title([fileNameNoPath,' - Voltage Contour'],'FontWeight','bold','FontSize',12);
    saveas(gcf,[DSSfilename,'_Voltage_Contour.fig'])
    

    %% Plot profile (Spider plot)
    figure;
    plotVoltageProfile(DSSCircObj)
    title([DSSfilename,' - Voltage Profile'],'FontWeight','bold','FontSize',12);
    saveas(gcf,[DSSfilename,'_Voltage_Spider.fig'])
    
end  


end