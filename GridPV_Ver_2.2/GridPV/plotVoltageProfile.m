%% plotVoltageProfile
% Plots the voltage profile for the feeder (spider plot)
%
%% Syntax
%  plotVoltageProfile(DSSCircObj);
%  plotVoltageProfile(DSSCircObj, _'PropertyName'_ ,PropertyValue);
%  Handles = plotVoltageProfile(DSSCircObj, _'PropertyName'_ ,PropertyValue);
%
%% Description
% Function to plot the voltage profile for the feeder.  This is the bus
% voltage vs. distance from the substation plot.  Also called a spider
% plot.  Clicking on objects in the figure will display the name of the object, and right clicking 
% will give a menu for viewing properties of the object.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|Properties|* - optional properties as one or more name-value pairs in any order
% * -- *|'SecondarySystem'|* - Property for if the secondary system (<600V)
% should be plotted (if it exists) |{'on'} | 'off'|
% * -- *|'Only3Phase'|* - Property for if only 3-phase power lines
% should be plotted |'on' | {'off'}|
% * -- *|'AveragePhase'|* - Property for if the average voltage should be
% plotted alone or in addition to the phase plots |'on' | {'off'} | 'addition'|
% * -- *|'LineToLine'|* - Property for if the voltage should be line-to-neutral or line-to-line |'on' | {'off'}|
% * -- *|'VoltScale'|* - Property for the y-axis voltage scale |{'120'} | 'pu'|
% * -- *|'DistanceScale'|* - Property for the x-axis distance scale |{'km'} | 'mi' | 'ft'|
% * -- *|'BackgroundShade'|* - Property for if the range of voltage values should be shaded as an area |'on' | {'off'}|
% * -- *|'BusName'|* - Property for the name of the bus (string) that the voltage
% profile should be plotted to.  Only the direct line between the bus and the substation
% will be plotted, unless all buses are selected. |{'all'} | busName|
% * -- *|'Downstream'|* - If a BusName is given, all buses in the electrical 
% path to the substation (upstream) will be plotted, and if this property is on, 
% all buses in the electrical path downstream of BusName will be plotted too |'on' | {'off'}|
% * -- *|'PVMarker'|* - Property for if the PV PCC should be marked |{'on'} | 'off'|
% * -- *|'CapacitorMarker'|* - Property for if capacitors should be marked |'on' | {'off'}|
% * -- *|'Lines'|* - Structure of the circuit lines from getLineInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Transformers'|* - Structure of the circuit transformers from getTransformerInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Capacitors'|* - Structure of the Capacitors from getCapacitorInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'PV'|* - Structure of the PV from getPVInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
%
%% Outputs
% * *|Handles|* - structure of handles for each type of object plotted in the figure
% * a figure is displayed with the plot
%
%% Notes
% For the right-click visualizations, the AllowForms
% field of DSSCircObj must be set to 1, which is the default value.
% Currently, OpenDSS 7.6.3 (the current version as of this writing) does not
% allow for setting the AllowForms field back to 1 after setting it to 0.
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
% Example of a feeder voltage profile plot
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% figure; plotVoltageProfile(DSSCircObj,'BusName','N292743','Downstream','on');
% figure; plotVoltageProfile(DSSCircObj);
% figure; plotVoltageProfile(DSSCircObj,'DistanceScale','ft','VoltScale','pu');
% figure; plotVoltageProfile(DSSCircObj,'SecondarySystem','off','AveragePhase','addition','Only3Phase','on');
% figure; plotVoltageProfile(DSSCircObj,'BackgroundShade','on','SecondarySystem','off');
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_PV_Distributed_7_5.dss"'];
% DSSText.command = 'Set mode=duty number=1  hour=12  h=1 sec=0';
% DSSText.command = 'Set controlmode=static';
% DSSText.command = 'solve';
% figure; plotVoltageProfile(DSSCircObj);
%

function Handles = plotVoltageProfile(DSSCircObj, varargin)

%% Parse input
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addParamValue('SecondarySystem', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('Only3Phase', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('AveragePhase', 'off', @(x)any(strcmp(x,{'on','off','addition'})));
p.addParamValue('LineToLine', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('VoltScale', '120', @(x)any(strcmp(x,{'120','pu'})));
p.addParamValue('DistanceScale', 'km', @(x)any(strcmp(x,{'km','mi','ft'})));
p.addParamValue('BackgroundShade', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('BusName', 'all', @ischar);
p.addParamValue('Downstream', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('PVMarker', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('CapacitorMarker', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('Lines', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('Transformers', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('Capacitors', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('PV', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));

p.parse(DSSCircObj, varargin{:}); %parse inputs


allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end


%% get circuit information
if ischar(Lines) && strcmp(Lines,'noInput')
    Lines = getLineInfo(DSSCircObj);
end
LinesToPlot = Lines;
if ischar(Transformers) && strcmp(Transformers,'noInput')
    Transformers = getTransformerInfo(DSSCircObj);
end
TransformersToPlot = Transformers;
DSSCircuit = DSSCircObj.ActiveCircuit;


%% Setup plotting scale
if strcmp(VoltScale, '120')
    voltPUscale = 120;
    voltScaleString = '120 V Base';
elseif strcmp(VoltScale, 'pu')
    voltPUscale = 1;
    voltScaleString = 'pu';
end

if strcmp(DistanceScale, 'km')
    distanceScale = 1;
elseif strcmp(DistanceScale, 'mi')
    distanceScale = 0.621371;
elseif strcmp(DistanceScale, 'ft')
    distanceScale = 3280.84;
end

%% Remove all 1 and 2 phase lines if Property Only3Phase is on
if strcmp(Only3Phase,'on')
    LinesToPlot = LinesToPlot([LinesToPlot.numPhases]==3);
    TransformersToPlot = TransformersToPlot([TransformersToPlot.numPhases]==3);
end


%% Remove all lines except between BusName and substation (if BusName is set)
if ~strcmp(BusName,'all')
    selectedBuses = findUpstreamBuses(DSSCircObj,BusName,'Lines',Lines,'Transformers',Transformers);
    if strcmp(Downstream,'on') %If Downstream property is on, also plot all buses downstream of BusName
        downstreamBuses = findDownstreamBuses(DSSCircObj,BusName,'Lines',Lines,'Transformers',Transformers);
        selectedBuses = [selectedBuses,downstreamBuses'];
    end
    linesBus1 = regexprep({LinesToPlot.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    linesBus2 = regexprep({LinesToPlot.bus2},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    LinesToPlot = LinesToPlot(ismember(linesBus1,selectedBuses) & ismember(linesBus2,selectedBuses));
    linesBus1 = regexprep({TransformersToPlot.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    linesBus2 = regexprep({TransformersToPlot.bus2},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    TransformersToPlot = TransformersToPlot(ismember(linesBus1,selectedBuses) & ismember(linesBus2,selectedBuses));
end


%% Plot main power lines
legendHandles = []; %empty variabes to store legend information
legendText = {};

% remove service lines
LinesNoSecondary = LinesToPlot([LinesToPlot.bus1Voltage]>600);

if strcmp(AveragePhase,'on') || strcmp(AveragePhase,'addition') %plot average voltage if requested in properties
    if strcmp(AveragePhase,'on')
        LineStyle = '-k';
    elseif strcmp(AveragePhase,'addition')
        LineStyle = '-g';
    end;
    h0 = plot([distanceScale*[LinesNoSecondary.bus1Distance]',distanceScale*[LinesNoSecondary.bus2Distance]']',[voltPUscale*[LinesNoSecondary.bus1VoltagePU]',voltPUscale*[LinesNoSecondary.bus2VoltagePU]']',LineStyle,'LineWidth',2);
    set(h0,{'DisplayName'},strcat('line.',{LinesNoSecondary.name})');
    hold all;
    
    % get step transformers and their voltages
    if ~isempty(TransformersToPlot) && ~strcmp(TransformersToPlot(1).name,'NONE')
        StepTransformers = TransformersToPlot([TransformersToPlot.bus2Voltage]>600);
        if ~isempty(StepTransformers)
            bus1Voltages = voltPUscale*[StepTransformers.bus1VoltagePU];
            bus2Voltages = voltPUscale*[StepTransformers.bus2VoltagePU];

            % plot
            condition = bus1Voltages>0 & bus2Voltages>0;
            h02 = plot([distanceScale*[StepTransformers(condition).bus1Distance]',distanceScale*[StepTransformers(condition).bus2Distance]']',[bus1Voltages(condition)',bus2Voltages(condition)']',LineStyle,'LineWidth',2);
            set(h02,{'DisplayName'},strcat('transformer.',{StepTransformers(condition).name})');
            h0 = [h0;h02];
        end
    else
        StepTransformers = [];
    end
    
    averageGroup = hggroup; set(h0,'Parent',averageGroup); %group together for one item in the legend
    legendHandles = [legendHandles,averageGroup];
    legendText = [legendText,'Average'];
    FigHandles.average_primary = h0;
end

if strcmp(AveragePhase,'off') || strcmp(AveragePhase,'addition') %plot phase voltages
    % find the phases on each line
    phaseA = ~cellfun(@isempty,regexp({LinesNoSecondary.bus1},'.^?(\.1)','Tokens'));
    phaseB = ~cellfun(@isempty,regexp({LinesNoSecondary.bus1},'.^?(\.2)','Tokens'));
    phaseC = ~cellfun(@isempty,regexp({LinesNoSecondary.bus1},'.^?(\.3)','Tokens'));
    isThreePhase = sum([phaseA;phaseB;phaseC])==0;
    phaseA(isThreePhase) = 1;
    phaseB(isThreePhase) = 1;
    phaseC(isThreePhase) = 1;
    
    %check for line to line
    if strcmp(LineToLine, 'off')
        bus1PhaseVoltages = voltPUscale*reshape([LinesNoSecondary.bus1PhaseVoltagesPU],3,[])';
        bus2PhaseVoltages = voltPUscale*reshape([LinesNoSecondary.bus2PhaseVoltagesPU],3,[])';
    else
        phaseA = phaseA & phaseB;
        phaseB = phaseB & phaseC;
        phaseC = phaseC & phaseA;
        bus1PhaseVoltages = voltPUscale*reshape([LinesNoSecondary.bus1PhaseVoltagesLLPU],3,[])';
        bus2PhaseVoltages = voltPUscale*reshape([LinesNoSecondary.bus2PhaseVoltagesLLPU],3,[])';
    end
    
    % plot
    h1 = plot([distanceScale*[LinesNoSecondary(phaseA).bus1Distance]',distanceScale*[LinesNoSecondary(phaseA).bus2Distance]']',[bus1PhaseVoltages(phaseA,1),bus2PhaseVoltages(phaseA,1)]','-k','LineWidth',2);
    set(h1,{'DisplayName'},strcat('line.',{LinesNoSecondary(phaseA).name})');
    hold all;
    h2 = plot([distanceScale*[LinesNoSecondary(phaseB).bus1Distance]',distanceScale*[LinesNoSecondary(phaseB).bus2Distance]']',[bus1PhaseVoltages(phaseB,2),bus2PhaseVoltages(phaseB,2)]','-r','LineWidth',2);
    set(h2,{'DisplayName'},strcat('line.',{LinesNoSecondary(phaseB).name})');
    h3 = plot([distanceScale*[LinesNoSecondary(phaseC).bus1Distance]',distanceScale*[LinesNoSecondary(phaseC).bus2Distance]']',[bus1PhaseVoltages(phaseC,3),bus2PhaseVoltages(phaseC,3)]','-b','LineWidth',2);
    set(h3,{'DisplayName'},strcat('line.',{LinesNoSecondary(phaseC).name})');
    
    % get step transformers and their voltages
    if ~isempty(TransformersToPlot) && ~strcmp(TransformersToPlot(1).name,'NONE')
        StepTransformers = TransformersToPlot([TransformersToPlot.bus2Voltage]>600);
        if ~isempty(StepTransformers)
            % find the phases on each transformer
            phaseA = cellfun(@(x)ismember(1,x),{StepTransformers.bus1NodeOrder});
            phaseB = cellfun(@(x)ismember(2,x),{StepTransformers.bus1NodeOrder});
            phaseC = cellfun(@(x)ismember(3,x),{StepTransformers.bus1NodeOrder});
            
            %check for line to line           
            if strcmp(LineToLine, 'off')
                bus1PhaseVoltages = voltPUscale*reshape([StepTransformers.bus1PhaseVoltagesPU],3,[])';
                bus2PhaseVoltages = voltPUscale*reshape([StepTransformers.bus2PhaseVoltagesPU],3,[])';
            else
                phaseA = phaseA & phaseB;
                phaseB = phaseB & phaseC;
                phaseC = phaseC & phaseA;
                bus1PhaseVoltages = voltPUscale*reshape([StepTransformers.bus1PhaseVoltagesLLPU],3,[])';
                bus2PhaseVoltages = voltPUscale*reshape([StepTransformers.bus2PhaseVoltagesLLPU],3,[])';
            end
            
            % plot
            h12 = plot([distanceScale*[StepTransformers(phaseA).bus1Distance]',distanceScale*[StepTransformers(phaseA).bus2Distance]']',[bus1PhaseVoltages(phaseA,1),bus2PhaseVoltages(phaseA,1)]','k','LineWidth',2);
            set(h12,{'DisplayName'},strcat('transformer.',{StepTransformers(phaseA).name})');
            h1 = [h1;h12];
            h22 = plot([distanceScale*[StepTransformers(phaseB).bus1Distance]',distanceScale*[StepTransformers(phaseB).bus2Distance]']',[bus1PhaseVoltages(phaseB,2),bus2PhaseVoltages(phaseB,2)]','r','LineWidth',2);
            set(h22,{'DisplayName'},strcat('transformer.',{StepTransformers(phaseB).name})');
            h2 = [h2;h22];
            h32 = plot([distanceScale*[StepTransformers(phaseC).bus1Distance]',distanceScale*[StepTransformers(phaseC).bus2Distance]']',[bus1PhaseVoltages(phaseC,3),bus2PhaseVoltages(phaseC,3)]','b','LineWidth',2);
            set(h32,{'DisplayName'},strcat('transformer.',{StepTransformers(phaseC).name})');
            h3 = [h3;h32];
        end
    else
        StepTransformers = [];
    end
    
    phaseAGroup = hggroup; set(h1,'Parent',phaseAGroup); %group together for one item in the legend
    phaseBGroup = hggroup; set(h2,'Parent',phaseBGroup); %group together for one item in the legend
    phaseCGroup = hggroup; set(h3,'Parent',phaseCGroup); %group together for one item in the legend
    FigHandles.phaseA_primary = h1;
    FigHandles.phaseB_primary = h2;
    FigHandles.phaseC_primary = h3;
    if strcmp(LineToLine, 'off')
        legendHandles = [legendHandles,phaseAGroup];
        legendText = [legendText,'PhaseA'];
        legendHandles = [legendHandles,phaseBGroup];
        legendText = [legendText,'PhaseB'];
        legendHandles = [legendHandles,phaseCGroup];
        legendText = [legendText,'PhaseC'];
    else
        legendHandles = [legendHandles,phaseAGroup];
        legendText = [legendText,'AB (L-L)'];
        legendHandles = [legendHandles,phaseBGroup];
        legendText = [legendText,'BC (L-L)'];
        legendHandles = [legendHandles,phaseCGroup];
        legendText = [legendText,'CA (L-L)'];
    end
end

%% Plot secondary system
secHandles = [];
if strcmp(SecondarySystem,'on')
    % service lines
    SecondaryLines = LinesToPlot([LinesToPlot.bus1Voltage]<600);
    
    if strcmp(AveragePhase,'on') || strcmp(AveragePhase,'addition') %plot average voltage if requested in properties
        h0sec = [];
        
        if strcmp(AveragePhase,'on')
            LineStyle = '--k';
        elseif strcmp(AveragePhase,'addition')
            LineStyle = '--g';
        end;
        if strcmp(LineToLine, 'on')
            h0sec = plot([distanceScale*[SecondaryLines.bus1Distance]',distanceScale*[SecondaryLines.bus2Distance]']',[voltPUscale*[SecondaryLines.bus1VoltagePU]',voltPUscale*[SecondaryLines.bus2VoltagePU]']',LineStyle,'LineWidth',2);
            set(h0sec,{'DisplayName'},strcat('line.',{SecondaryLines.name})');
        else
            h0sec = plot([distanceScale*[SecondaryLines.bus1Distance]',distanceScale*[SecondaryLines.bus2Distance]']',[voltPUscale*[SecondaryLines.bus1VoltageLLPU]',voltPUscale*[SecondaryLines.bus2VoltageLLPU]']',LineStyle,'LineWidth',2);
            set(h0sec,{'DisplayName'},strcat('line.',{SecondaryLines.name})');
        end
        hold all;
        
        % get secondary transformers and their voltages
        if ~isempty(TransformersToPlot) && ~strcmp(TransformersToPlot(1).name,'NONE')
            SecondaryTransformers = TransformersToPlot([TransformersToPlot.bus2Voltage]<600);
            if ~isempty(SecondaryTransformers)
                %check for line to line
                if strcmp(LineToLine, 'off')
                    bus1Voltages = voltPUscale*[SecondaryTransformers.bus1VoltagePU];
                    bus2Voltages = voltPUscale*[SecondaryTransformers.bus2VoltagePU];
                else
                    bus1Voltages = voltPUscale*[SecondaryTransformers.bus1VoltageLLPU];
                    bus2Voltages = voltPUscale*[SecondaryTransformers.bus2VoltageLLPU];
                end
                
                % plot
                condition = bus1Voltages>0 & bus2Voltages>0;
                h0sec2 = plot([distanceScale*[SecondaryTransformers(condition).bus1Distance]',distanceScale*[SecondaryTransformers(condition).bus2Distance]']',[bus1Voltages(condition)',bus2Voltages(condition)']',LineStyle,'LineWidth',2);
                set(h0sec2,{'DisplayName'},strcat('transformer.',{SecondaryTransformers(condition).name})');
                h0sec = [h0sec;h0sec2];
            end
        else
            SecondaryTransformers = [];
        end
        averageGroupSec = hggroup; set(h0sec,'Parent',averageGroupSec); %group together for one item in the legend
        FigHandles.average_secondary = h0sec;
        secHandles = [secHandles; h0sec];
    end

    if strcmp(AveragePhase,'off') || strcmp(AveragePhase,'addition') %plot phase voltages
        h1sec = [];
        h2sec = [];
        h3sec = [];
        
        % find the phases on each line
        phaseA = ~cellfun(@isempty,regexp({SecondaryLines.bus1},'.^?(\.1)','Tokens'));
        phaseB = ~cellfun(@isempty,regexp({SecondaryLines.bus1},'.^?(\.2)','Tokens'));
        phaseC = ~cellfun(@isempty,regexp({SecondaryLines.bus1},'.^?(\.3)','Tokens'));
        isThreePhase = sum([phaseA;phaseB;phaseC])==0;
        phaseA(isThreePhase) = 1;
        phaseB(isThreePhase) = 1;
        phaseC(isThreePhase) = 1;
        
        %check for line to line
        if strcmp(LineToLine, 'off')
            bus1PhaseVoltages = voltPUscale*reshape([SecondaryLines.bus1PhaseVoltagesPU],3,[])';
            bus2PhaseVoltages = voltPUscale*reshape([SecondaryLines.bus2PhaseVoltagesPU],3,[])';
        else
            phaseA = phaseA & phaseB;
            phaseB = phaseB & phaseC;
            phaseC = phaseC & phaseA;
            bus1PhaseVoltages = voltPUscale*reshape([SecondaryLines.bus1PhaseVoltagesLLPU],3,[])';
            bus2PhaseVoltages = voltPUscale*reshape([SecondaryLines.bus2PhaseVoltagesLLPU],3,[])';
        end
                
        % plot
        if ~isempty(SecondaryLines)
            h1sec = plot([distanceScale*[SecondaryLines(phaseA).bus1Distance]',distanceScale*[SecondaryLines(phaseA).bus2Distance]']',[bus1PhaseVoltages(phaseA,1),bus2PhaseVoltages(phaseA,1)]','--k','LineWidth',2);
            set(h1sec,{'DisplayName'},strcat('line.',{SecondaryLines(phaseA).name})');
            h2sec = plot([distanceScale*[SecondaryLines(phaseB).bus1Distance]',distanceScale*[SecondaryLines(phaseB).bus2Distance]']',[bus1PhaseVoltages(phaseB,2),bus2PhaseVoltages(phaseB,2)]','--r','LineWidth',2);
            set(h2sec,{'DisplayName'},strcat('line.',{SecondaryLines(phaseB).name})');
            h3sec = plot([distanceScale*[SecondaryLines(phaseC).bus1Distance]',distanceScale*[SecondaryLines(phaseC).bus2Distance]']',[bus1PhaseVoltages(phaseC,3),bus2PhaseVoltages(phaseC,3)]','--b','LineWidth',2);
            set(h3sec,{'DisplayName'},strcat('line.',{SecondaryLines(phaseC).name})');
        end
        
        % get secondary transformers and their voltages
        if ~isempty(TransformersToPlot) && ~strcmp(TransformersToPlot(1).name,'NONE')
            SecondaryTransformers = TransformersToPlot([TransformersToPlot.bus2Voltage]<600);
            if ~isempty(SecondaryTransformers)
                % find the phases on each transformer
                phaseA = ~cellfun(@isempty,regexp({SecondaryTransformers.bus1},'.^?(\.1)','Tokens'));
                phaseB = ~cellfun(@isempty,regexp({SecondaryTransformers.bus1},'.^?(\.2)','Tokens'));
                phaseC = ~cellfun(@isempty,regexp({SecondaryTransformers.bus1},'.^?(\.3)','Tokens'));
                isThreePhase = sum([phaseA;phaseB;phaseC])==0;
                phaseA(isThreePhase) = 1;
                phaseB(isThreePhase) = 1;
                phaseC(isThreePhase) = 1;
                
                %check for line to line
                if strcmp(LineToLine, 'off')
                    phaseA2 = ~cellfun(@isempty,regexp({SecondaryTransformers.bus2},'.^?(\.1)','Tokens'));
                    phaseB2 = ~cellfun(@isempty,regexp({SecondaryTransformers.bus2},'.^?(\.2)','Tokens'));
                    phaseC2 = ~cellfun(@isempty,regexp({SecondaryTransformers.bus2},'.^?(\.3)','Tokens'));
                    isThreePhase = sum([phaseA2;phaseB2;phaseC2])==0;
                    phaseA2(isThreePhase) = 1;
                    phaseB2(isThreePhase) = 1;
                    phaseC2(isThreePhase) = 1;
                    
                    phaseA = phaseA & phaseA2;
                    phaseB = phaseB & phaseB2;
                    phaseC = phaseC & phaseC2;
                    
                    bus1PhaseVoltages = voltPUscale*reshape([SecondaryTransformers.bus1PhaseVoltagesPU],3,[])';
                    bus2PhaseVoltages = voltPUscale*reshape([SecondaryTransformers.bus2PhaseVoltagesPU],3,[])';
                else
                    phaseA = phaseA & phaseB;
                    phaseB = phaseB & phaseC;
                    phaseC = phaseC & phaseA;
                    bus1PhaseVoltages = voltPUscale*reshape([SecondaryTransformers.bus1PhaseVoltagesLLPU],3,[])';
                    bus2PhaseVoltages = voltPUscale*reshape([SecondaryTransformers.bus2PhaseVoltagesLLPU],3,[])';
                end
                
                % plot
                h1sec2 = plot([distanceScale*[SecondaryTransformers(phaseA).bus1Distance]',distanceScale*[SecondaryTransformers(phaseA).bus2Distance]']',[bus1PhaseVoltages(phaseA,1),bus2PhaseVoltages(phaseA,1)]','--k','LineWidth',2);
                set(h1sec2,{'DisplayName'},strcat('transformer.',{SecondaryTransformers(phaseA).name})');
                h1sec = [h1sec; h1sec2];
                h2sec2 = plot([distanceScale*[SecondaryTransformers(phaseB).bus1Distance]',distanceScale*[SecondaryTransformers(phaseB).bus2Distance]']',[bus1PhaseVoltages(phaseB,2),bus2PhaseVoltages(phaseB,2)]','--r','LineWidth',2);
                set(h2sec2,{'DisplayName'},strcat('transformer.',{SecondaryTransformers(phaseB).name})');
                h2sec = [h2sec; h2sec2];
                h3sec2 = plot([distanceScale*[SecondaryTransformers(phaseC).bus1Distance]',distanceScale*[SecondaryTransformers(phaseC).bus2Distance]']',[bus1PhaseVoltages(phaseC,3),bus2PhaseVoltages(phaseC,3)]','--b','LineWidth',2);
                set(h3sec2,{'DisplayName'},strcat('transformer.',{SecondaryTransformers(phaseC).name})');
                h3sec = [h3sec; h3sec2];
            end
        else
            SecondaryTransformers = [];
        end
        if ~isempty(SecondaryLines) & ~isempty(h1sec)
            phaseAGroupSec = hggroup; set(h1sec,'Parent',phaseAGroupSec); %group together for one item in the legend
            FigHandles.phaseA_secondary = h1sec;
        end
        if ~isempty(SecondaryLines) & ~isempty(h1sec)
            phaseBGroupSec = hggroup; set(h2sec,'Parent',phaseBGroupSec); %group together for one item in the legend
            FigHandles.phaseB_secondary = h2sec;
        end
        if ~isempty(SecondaryLines) & ~isempty(h1sec)
            phaseCGroupSec = hggroup; set(h3sec,'Parent',phaseCGroupSec); %group together for one item in the legend
            FigHandles.phaseC_secondary = h3sec;
        end
        secHandles = [secHandles; h1sec; h2sec; h3sec];
    end
end


%% Place the background shading on the figure
if strcmp(BackgroundShade,'on')
    
    %Filter transmission lines
    LinesMV = LinesToPlot([LinesToPlot.bus1Voltage]<21000);
    
    %Divide bins by the total length of the feeder
    maxDistance = max([[LinesMV.bus1Distance],[LinesMV.bus2Distance]]);
    distanceSpacing = maxDistance/50;
    
    %Get bus voltage and distances
    if strcmp(LineToLine, 'off')
        busPhaseVoltages = [vertcat(LinesMV.bus1PhaseVoltagesPU);vertcat(LinesMV.bus2PhaseVoltagesPU)];
    else
        busPhaseVoltages = [vertcat(LinesMV.bus1PhaseVoltagesLLPU);vertcat(LinesMV.bus2PhaseVoltagesLLPU)];
    end
    busDistances = [vertcat(LinesMV.bus1Distance);vertcat(LinesMV.bus2Distance)];
    
    %Add buses for long lines
    LongLines = LinesMV(abs([LinesMV.bus2Distance]-[LinesMV.bus1Distance])>distanceSpacing);
    lineLength = abs([LongLines.bus2Distance]-[LongLines.bus1Distance]);
    lineStart = min([[LongLines.bus2Distance];[LongLines.bus1Distance]]);
    lineEnd = max([[LongLines.bus2Distance];[LongLines.bus1Distance]]);
    for ii=1:length(LongLines)
        for jj=lineStart(ii):distanceSpacing:lineEnd(ii)
            if strcmp(LineToLine, 'off')
                busPhaseVoltages = [busPhaseVoltages; LongLines(ii).bus1PhaseVoltagesPU+(LongLines(ii).bus2PhaseVoltagesPU-LongLines(ii).bus1PhaseVoltagesPU)*abs(jj-LongLines(ii).bus1Distance)/lineLength(ii)];
            else
                busPhaseVoltages = [busPhaseVoltages; LongLines(ii).bus1PhaseVoltagesLLPU+(LongLines(ii).bus2PhaseVoltagesLLPU-LongLines(ii).bus1PhaseVoltagesLLPU)*abs(jj-LongLines(ii).bus1Distance)/lineLength(ii)];
            end
            busDistances = [busDistances; jj];
        end
    end
    
    busPhaseVoltages(busPhaseVoltages==0) = NaN;
    
    %Split into bins and find the max and min in each bin
    index = 1;
    for ii=0:distanceSpacing:maxDistance
        phaseVoltages = busPhaseVoltages(busDistances>ii-distanceSpacing/2 & busDistances<ii+distanceSpacing/2,:);
        tempDistances = busDistances(busDistances>ii-distanceSpacing/2 & busDistances<ii+distanceSpacing/2);
        tempMin = min(phaseVoltages,[],2);
        if ~isempty(tempMin)
            [minVoltage(index) minIndex] = min(tempMin);
            minDistances(index) = tempDistances(minIndex);
        end
        tempMax = max(phaseVoltages,[],2);
        if ~isempty(tempMax)
            [maxVoltage(index) maxIndex] = max(tempMax);
            maxDistances(index) = tempDistances(maxIndex);
            index = index+1;
        end
    end
    hold all;
    fillhandle=fill([maxDistances,fliplr(minDistances)],[maxVoltage,fliplr(minVoltage)]*120,[0.992 0.918 0.796],'EdgeColor',[0.5 0.5 0.5]);
    uistack(fillhandle,'bottom');
    
    legendHandles = [legendHandles,fillhandle];
    legendText = [legendText,'Service Range'];
    FigHandles.BackgroundShade = fillhandle;
end


%% Place PV marker in the figure
if strcmp(PVMarker,'on')
    if ischar(PV) && strcmp(PV,'noInput')
        PV = getPVInfo(DSSCircObj);
    end
    
    %only plot PV markers that are attached to lines being plotted
    if ~isempty(PV) && ~strcmp(PV(1).name,'NONE')
        busesToConsider = [regexprep({LinesNoSecondary.bus1},'(\.[0-9]+)',''),regexprep({LinesNoSecondary.bus2},'(\.[0-9]+)','')];
        if ~isempty(StepTransformers)
            busesToConsider = [busesToConsider, regexprep({StepTransformers.bus1},'(\.[0-9]+)',''),regexprep({StepTransformers.bus2},'(\.[0-9]+)','')];
        end
        if strcmp(SecondarySystem,'on')
            busesToConsider = [busesToConsider, regexprep({SecondaryLines.bus1},'(\.[0-9]+)',''),regexprep({SecondaryLines.bus2},'(\.[0-9]+)','')];
            if ~isempty(SecondaryTransformers)
                busesToConsider = [busesToConsider, regexprep({SecondaryTransformers.bus1},'(\.[0-9]+)',''),regexprep({SecondaryTransformers.bus2},'(\.[0-9]+)','')];
            end
        end
        PV = PV(ismember(regexprep({PV.busName},'(\.[0-9]+)',''),busesToConsider)); 
    end
    
    if ~isempty(PV) && ~strcmp(PV(1).name,'NONE')
        pvHandles = [];
        if strcmp(AveragePhase,'on') || strcmp(AveragePhase,'addition') %plot average voltage if requested in properties
            if strcmp(LineToLine, 'on')
                pvHandles = plot([distanceScale*[PV.distance]',distanceScale*[PV.distance]']',[voltPUscale*[PV.voltageLLPU]',voltPUscale*[PV.voltageLLPU]']','-kh','MarkerSize',11,'MarkerFaceColor','y','LineStyle','none');
                set(pvHandles,{'DisplayName'},strcat('pvsystem.',{PV.name})');
            else
                pvHandles = plot([distanceScale*[PV.distance]',distanceScale*[PV.distance]']',[voltPUscale*[PV.voltagePU]',voltPUscale*[PV.voltagePU]']','-kh','MarkerSize',11,'MarkerFaceColor','y','LineStyle','none');
                set(pvHandles,{'DisplayName'},strcat('pvsystem.',{PV.name})');
            end
        end
        
        if strcmp(AveragePhase,'off') || strcmp(AveragePhase,'addition') %plot phase voltages
            % find the phases on each PV
            phaseA = ~cellfun(@isempty,regexp({PV.busName},'.^?(\.1)','Tokens'));
            phaseB = ~cellfun(@isempty,regexp({PV.busName},'.^?(\.2)','Tokens'));
            phaseC = ~cellfun(@isempty,regexp({PV.busName},'.^?(\.3)','Tokens'));
            isThreePhase = sum([phaseA;phaseB;phaseC])==0;
            phaseA(isThreePhase) = 1;
            phaseB(isThreePhase) = 1;
            phaseC(isThreePhase) = 1;
            
            %check for line to line
            if strcmp(LineToLine, 'off')
                phaseVoltages = voltPUscale*reshape([PV.phaseVoltagesPU],3,[])';
            else
                phaseA = phaseA & phaseB;
                phaseB = phaseB & phaseC;
                phaseC = phaseC & phaseA;
                phaseVoltages = voltPUscale*reshape([PV.phaseVoltagesLLPU],3,[])';
            end
            
            pvHandles1 = plot([distanceScale*[PV(phaseA).distance]',distanceScale*[PV(phaseA).distance]']',[phaseVoltages(phaseA,1),phaseVoltages(phaseA,1)]','-kh','MarkerSize',11,'MarkerFaceColor','y','LineStyle','none');
            set(pvHandles1,{'DisplayName'},strcat('pvsystem.',{PV(phaseA).name})');
            pvHandles2 = plot([distanceScale*[PV(phaseB).distance]',distanceScale*[PV(phaseB).distance]']',[phaseVoltages(phaseB,2),phaseVoltages(phaseB,2)]','-kh','MarkerSize',11,'MarkerFaceColor','y','LineStyle','none');
            set(pvHandles2,{'DisplayName'},strcat('pvsystem.',{PV(phaseB).name})');
            pvHandles3 = plot([distanceScale*[PV(phaseC).distance]',distanceScale*[PV(phaseC).distance]']',[phaseVoltages(phaseC,3),phaseVoltages(phaseC,3)]','-kh','MarkerSize',11,'MarkerFaceColor','y','LineStyle','none');
            set(pvHandles3,{'DisplayName'},strcat('pvsystem.',{PV(phaseC).name})');
            pvHandles = [pvHandles1; pvHandles2; pvHandles3];
        end
        pvGroup = hggroup; set(pvHandles,'Parent',pvGroup); %group together for one item in the legend
        legendHandles = [legendHandles,pvGroup];
        legendText = [legendText,'PV PCC'];
        FigHandles.PV = pvHandles;
    end
end


%% Place Capacitor marker in the figure
if strcmp(CapacitorMarker,'on')
    if ischar(Capacitors) && strcmp(Capacitors,'noInput')
        Capacitors = getCapacitorInfo(DSSCircObj);
    end
    
    %only plot Capacitor markers that are attached to lines being plotted
    if ~isempty(Capacitors) && ~strcmp(Capacitors(1).name,'NONE')
        busesToConsider = [regexprep({LinesNoSecondary.bus1},'(\.[0-9]+)',''),regexprep({LinesNoSecondary.bus2},'(\.[0-9]+)','')];
        if ~isempty(StepTransformers)
            busesToConsider = [busesToConsider, regexprep({StepTransformers.bus1},'(\.[0-9]+)',''),regexprep({StepTransformers.bus2},'(\.[0-9]+)','')];
        end
        if strcmp(SecondarySystem,'on')
            busesToConsider = [busesToConsider, regexprep({SecondaryLines.bus1},'(\.[0-9]+)',''),regexprep({SecondaryLines.bus2},'(\.[0-9]+)','')];
            if ~isempty(SecondaryTransformers)
                busesToConsider = [busesToConsider, regexprep({SecondaryTransformers.bus1},'(\.[0-9]+)',''),regexprep({SecondaryTransformers.bus2},'(\.[0-9]+)','')];
            end
        end
        Capacitors = Capacitors(ismember(regexprep({Capacitors.busName},'(\.[0-9]+)',''),busesToConsider)); 
    end
    
    if ~isempty(Capacitors) && ~strcmp(Capacitors(1).name,'NONE')
        capacitorHandles = [];
        if strcmp(AveragePhase,'on') || strcmp(AveragePhase,'addition') %plot average voltage if requested in properties
            if strcmp(LineToLine, 'on')
                capacitorHandles = plot([distanceScale*[Capacitors.distance]',distanceScale*[Capacitors.distance]']',[voltPUscale*[Capacitors.voltageLLPU]',voltPUscale*[Capacitors.voltageLLPU]']','-ks','MarkerSize',10,'MarkerFaceColor','g','LineStyle','none');
                set(capacitorHandles,{'DisplayName'},strcat('capacitor.',{Capacitors.name})');
            else
                capacitorHandles = plot([distanceScale*[Capacitors.distance]',distanceScale*[Capacitors.distance]']',[voltPUscale*[Capacitors.voltagePU]',voltPUscale*[Capacitors.voltagePU]']','-ks','MarkerSize',10,'MarkerFaceColor','g','LineStyle','none');
                set(capacitorHandles,{'DisplayName'},strcat('capacitor.',{Capacitors.name})');
            end
        end
        
        if strcmp(AveragePhase,'off') || strcmp(AveragePhase,'addition') %plot phase voltages
            % find the phases on each Capacitors
            phaseA = ~cellfun(@isempty,regexp({Capacitors.busName},'.^?(\.1)','Tokens'));
            phaseB = ~cellfun(@isempty,regexp({Capacitors.busName},'.^?(\.2)','Tokens'));
            phaseC = ~cellfun(@isempty,regexp({Capacitors.busName},'.^?(\.3)','Tokens'));
            isThreePhase = sum([phaseA;phaseB;phaseC])==0;
            phaseA(isThreePhase) = 1;
            phaseB(isThreePhase) = 1;
            phaseC(isThreePhase) = 1;
            
            %check for line to line
            if strcmp(LineToLine, 'off')
                phaseVoltages = voltPUscale*reshape([Capacitors.phaseVoltagesPU],3,[])';
            else
                phaseA = phaseA & phaseB;
                phaseB = phaseB & phaseC;
                phaseC = phaseC & phaseA;
                phaseVoltages = voltPUscale*reshape([Capacitors.phaseVoltagesLLPU],3,[])';
            end
            
            capacitorHandles1 = plot([distanceScale*[Capacitors(phaseA).distance]',distanceScale*[Capacitors(phaseA).distance]']',[phaseVoltages(phaseA,1),phaseVoltages(phaseA,1)]','-ks','MarkerSize',10,'MarkerFaceColor','g','LineStyle','none');
            set(capacitorHandles1,{'DisplayName'},strcat('capacitor.',{Capacitors(phaseA).name})');
            capacitorHandles2 = plot([distanceScale*[Capacitors(phaseB).distance]',distanceScale*[Capacitors(phaseB).distance]']',[phaseVoltages(phaseB,2),phaseVoltages(phaseB,2)]','-ks','MarkerSize',10,'MarkerFaceColor','g','LineStyle','none');
            set(capacitorHandles2,{'DisplayName'},strcat('capacitor.',{Capacitors(phaseB).name})');
            capacitorHandles3 = plot([distanceScale*[Capacitors(phaseC).distance]',distanceScale*[Capacitors(phaseC).distance]']',[phaseVoltages(phaseC,3),phaseVoltages(phaseC,3)]','-ks','MarkerSize',10,'MarkerFaceColor','g','LineStyle','none');
            set(capacitorHandles3,{'DisplayName'},strcat('capacitor.',{Capacitors(phaseC).name})');
            capacitorHandles = [capacitorHandles1;capacitorHandles2;capacitorHandles3];
        end
        CapacitorsGroup = hggroup; set(capacitorHandles,'Parent',CapacitorsGroup); %group together for one item in the legend
        legendHandles = [legendHandles,CapacitorsGroup];
        legendText = [legendText,'Capacitor'];
        FigHandles.Capacitors = capacitorHandles;
    end
end


%% Plot Edits
grid on;
set(gca,'FontSize',10,'FontWeight','bold')
xlabel(gca,['Distance from Substation (', DistanceScale ,')'],'FontSize',12,'FontWeight','bold')
ylabel(gca,['Bus Voltage (' voltScaleString ')'],'FontSize',12,'FontWeight','bold')
if strcmp(LineToLine, 'off')
    ylabel(gca,['Bus Voltage (' voltScaleString ')'],'FontSize',12,'FontWeight','bold')
    title('Feeder Voltage Profile','FontWeight','bold','FontSize',12);
else
    ylabel(gca,['Bus Voltage L-L (' voltScaleString ')'],'FontSize',12,'FontWeight','bold')
    title('Feeder Voltage Profile (Line to Line)','FontWeight','bold','FontSize',12);
end
lgh = legend(legendHandles,legendText);
FigHandles.legend = lgh;
FigHandles.legendHandles = legendHandles;
FigHandles.legendText = legendText;

%drawnow;pause(0.1);


%% Add button to toolbar
if isempty(findall(gcf,'tag','NodeViewToggle'))
    tbh = findall(gcf,'Type','uitoolbar');
    buttonIcon = ones(15,15); buttonIcon(6:9,2:5)=0; buttonIcon(6:9,11:14)=0; buttonIcon(7:8,5:11)=0; buttonIcon = repmat(buttonIcon,[1,1,3]);
    if ~isempty(tbh)
        tth = uitoggletool(tbh,'CData',buttonIcon,'Separator','on','Tag','NodeViewToggle','TooltipString','Turn On Node View');
    end
end
if isempty(findall(gcf,'tag','SecondariesToggle'))
    set(gca,'YLimMode','manual')
    tbh = findall(gcf,'Type','uitoolbar');
    buttonIconSec = ones(15,15); buttonIconSec(1:12,7:9)=0; buttonIconSec(12:14,7:15)=0; buttonIconSec = repmat(buttonIconSec,[1,1,3]);
    if strcmp(SecondarySystem,'on')
        commandString = 'if strcmp(get(findall(gcf,''tag'',''SecondariesToggle''),''State''),''off''), set(secHandles,''Visible'',''off'');, else, set(secHandles,''Visible'',''on'');, end';
        buttonState = 'on';
    else
        commandString = 'if strcmp(get(findall(gcf,''tag'',''SecondariesToggle''),''State''),''on''), errordlg(''Cannnot turn on the plotting of the secondary system without the option initially turned on in the command prompt. Please use the property ''''SecondarySystem'''', ''''on'''' in plotVoltageProfile() to use this functionality.'');, end';
        buttonState = 'off';
    end
    if ~isempty(tbh)
        tth2 = uitoggletool(tbh,'CData',buttonIconSec,'Separator','on','Tag','SecondariesToggle','TooltipString','Turn On/Off Secondaries','State',buttonState,'ClickedCallback',['temp = get(gca,''UserData''); secHandles=temp{4}; ' commandString]);
    end
end


%% make button down functions

%remove previous tooltip strings
removeString = 'if ~isempty(findobj(gca,''tag'',''Tooltip'')) delete(findobj(gca,''tag'',''Tooltip'')); end; set(findobj(gca,''Selected'',''on''),''Selected'',''off''); ';
%add in either one text box for the object name, or if in node view, add two text boxes for each bus name
addText = 'if strcmp(get(findall(gcf,''tag'',''NodeViewToggle''),''State''),''off''), text(max(xCoords),max(yCoords),['' '',get(gco,''DisplayName'')],''Interpreter'',''none'',''Tag'',''Tooltip'',''FontWeight'',''bold'',''BackgroundColor'',[1 1 1]);, ';
addText = [addText 'else, temp = get(gca,''UserData''); DSSCircObj=temp{1}; if DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName''))>0, bus1=DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1}; bus2=DSSCircObj.ActiveCircuit.ActiveElement.BusNames{end}; '];
addText = [addText 'text(xCoords(1),yCoords(1),['' '',bus1],''Interpreter'',''none'',''Tag'',''Tooltip'',''FontWeight'',''bold'',''BackgroundColor'',[1 1 1]); textH=(findobj(gca,''tag'',''Tooltip'')); '];
addText = [addText 'text(xCoords(2),yCoords(2),['' '',bus2],''Interpreter'',''none'',''Tag'',''Tooltip'',''FontWeight'',''bold'',''BackgroundColor'',[1 1 1]); textH=(findobj(gca,''tag'',''Tooltip'')); extent=get(textH(1),''Extent''); set(textH(2),''Position'',[xCoords(1) yCoords(1)+extent(4)/1.8]); set(textH(1),''Position'',[xCoords(2) yCoords(2)-extent(4)/1.8]);, end, end'];

%make buttondownfunction of all lines and markers display the name of the object when clicked
set(get(gca,'Children'),'buttondownfcn',[removeString,'set(gco,''Selected'',''on''); xCoords = get(gco,''XData''); yCoords = get(gco,''YData'');' addText]);
%make buttondownfuction for all items that were grouped together
groupChildren =get(findobj(gca,'Type','hggroup'),'Children');
if iscell(groupChildren) groupChildren=vertcat(groupChildren{:}); end 
set(groupChildren,'buttondownfcn',[removeString,'set(gco,''Selected'',''on''); xCoords = get(gco,''XData''); yCoords = get(gco,''YData'');' addText]);
%if clicking not on an object, remove any text objects and selections
set(gca,'buttondownfcn',removeString)


%% Right Click button
% Set AllowForms to true to allow visualizations (for future versions of OpenDSS)
DSSCircObj.AllowForms = 1;

set(gca,'UserData',{DSSCircObj,Lines,Transformers,secHandles});
hcmenu = uicontextmenu;
hcb1 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; if DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName''))>0 DSSCircObj.Text.Command = ''formedit''; else errordlg(sprintf(''Did not find element %s in the circuit.'',get(gco,''DisplayName''))); end'];
hcb2 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.Text.Command = [''visualize voltages element='', get(gco,''DisplayName'')];'];
hcb3 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.Text.Command = [''visualize currents element='', get(gco,''DisplayName'')];'];
hcb4 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.Text.Command = [''visualize powers element='', get(gco,''DisplayName'')];'];
hcb5 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; Lines=temp{2}; Transformers=temp{3}; disp(''Plotting circuit lines with marker...''); objectName=get(gco,''DisplayName''); DSSCircObj.ActiveCircuit.SetActiveElement(objectName); figure; plotCircuitLines(DSSCircObj,''CustomMarker'',DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1},''CustomLegend'',objectName,''Lines'',Lines,''Transformers'',Transformers);'];
hcb6 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; Lines=temp{2}; disp(''Plotting amp profile...''); DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName'')); figure; plotAmpProfile(DSSCircObj,DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1},''Lines'',Lines);'];
hcb7 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; Lines=temp{2}; disp(''Plotting voltage profile...''); DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName'')); figure; plotVoltageProfile(DSSCircObj,''BusName'',DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1},''Lines'',Lines);'];
hcb8 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; Lines=temp{2}; disp(''Plotting kW profile...''); DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName'')); figure; plotKWProfile(DSSCircObj,''BusName'',DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1},''Lines'',Lines);'];
hcb9 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; Lines=temp{2}; disp(''Plotting kVAr profile...''); DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName'')); figure; plotKVARProfile(DSSCircObj,''BusName'',DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1},''Lines'',Lines);'];
hcb10 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName'')); DSSCircObj.Text.Command = sprintf(''Show busflow %s kva elem'',regexprep(DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1},''(\.[0-9]+)'',''''));'];
hcb11 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName'')); DSSCircObj.Text.Command = sprintf(''Show busflow %s kva elem'',regexprep(DSSCircObj.ActiveCircuit.ActiveElement.BusNames{end},''(\.[0-9]+)'',''''));'];

item1 = uimenu(hcmenu, 'Label', 'View Properties', 'Callback', hcb1);
item2 = uimenu(hcmenu, 'Label', 'View Voltages', 'Callback', hcb2);
item3 = uimenu(hcmenu, 'Label', 'View Currents', 'Callback', hcb3);
item4 = uimenu(hcmenu, 'Label', 'View Powers', 'Callback', hcb4);
item5 = uimenu(hcmenu, 'Label', 'View Circuit Plot with Element Marked', 'Callback', hcb5);
item6 = uimenu(hcmenu, 'Label', 'View Amp Profile to Element', 'Callback', hcb6);
item7 = uimenu(hcmenu, 'Label', 'View Voltage Profile to Element', 'Callback', hcb7);
item8 = uimenu(hcmenu, 'Label', 'View kW Profile to Element', 'Callback', hcb8);
item9 = uimenu(hcmenu, 'Label', 'View kVAr Profile to Element', 'Callback', hcb9);
item10 = uimenu(hcmenu, 'Label', 'View Bus Flow for Bus1', 'Callback', hcb10);
item11 = uimenu(hcmenu, 'Label', 'View Bus Flow for Bus2', 'Callback', hcb11);
set(get(gca,'Children'),'UIContextMenu',hcmenu)
set(groupChildren,'UIContextMenu',hcmenu)


%% Setup function outputs
if nargout>0
    Handles = FigHandles;
end


end