%% plotCircuitLines
% Plots the feeder circuit diagram
%
%% Syntax
%  plotCircuitLines(DSSCircObj);
%  plotCircuitLines(DSSCircObj, _'PropertyName'_ ,PropertyValue);
%  Handles = plotCircuitLines(DSSCircObj, _'PropertyName'_ ,PropertyValue);
%
%% Description
% Function to plot the feeder circuit diagram.  The coloring and line
% thickness plotting styles can be customized by the user through the function
% property inputs.  If no properties are selected, the plotCircuitLinesOptions 
% GUI window is displayed to assist the user is selecting plotting options.  Clicking on objects in
% the figure will display the name of the object, and right clicking will give a menu for viewing
% properties of the object.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|Properties|* - optional properties as one or more name-value pairs in any order
% * -- *|'Coloring'|* - Defines how the circuit lines are colored in the figure.
% * ------ |ColorSpec| - three-element RGB vector specifying the line color
% * ------ |'numPhases'| - black for 3-phase lines and a light gray for 1 or 2 phase lines. This is the default
% * ------ |'perPhase'| - colors each phase (or combination of phases) a different color in the figure
% * ------ |'energyMeter'| - colors each energy meter zone a different color in the figure
% * ------ |'voltage120'| - contours the line colors according to the voltage on a 120V base
% * ------ |'voltagePU'| - contours the line colors according to the per unit voltage
% * ------ |'voltage'| - contours the line colors according to the voltage (kV)
% * ------ |'voltage120LL'| - contours the line colors according to the line-to-line voltage on a 120V base
% * ------ |'voltagePULL'| - contours the line colors according to the line-to-line per unit voltage
% * ------ |'voltageLL'| - contours the line colors according to the line-to-line voltage (kV)
% * ------ |'lineLoading'| - contours the line colors according to the line loading (current/line rating)
% * ------ |'realLosses'| - contours the line colors according to the real power line losses (kW/km)
% * ------ |'reactiveLosses'| - contours the line colors according to the reactive power line losses (kVAR/km)
% * ------ |'distance'| - contours the line colors according to the distance from the substation
% * ------ |'unbalance'| - contours the line colors according to the power (kVA) unbalance between phases
% * ------ |'voltageAngle'| - contours the line colors according to the angle of the bus voltage phasor
% * ------ |'powerFactor'| - contours the line colors according to the power factor of the power flow
% * ------ |'powerFlowDirection'| - contours the line colors according to the line power flow (kW) with separate colors for upstream and downstream flow
% * ------ |'impedance'| - contours the line colors according to the positive-sequence short-circuit impedance magnitude
% * ------ |'resistance'| - contours the line colors according to the positive-sequence short-circuit resistance
% * ------ |'reactance'| - contours the line colors according to the positive-sequence short-circuit reactance
% * ------ |'faultCurrent3P'| - contours the line colors according to the fault current for a 3-phase fault
% * ------ |'faultCurrent1P'| - contours the line colors according to the fault current for a 1-phase fault
% * ------ |'faultCurrentLL'| - contours the line colors according to the fault current for a Line-to-Line fault
% * -- *|'ContourScale'|* - Defines the minimum and maximum value for contouring or auto scaling |{'auto'} | [0 5]|
% * -- *|'Thickness'|* - Defines how the thickness of the circuit lines is displayed.  |{'numPhases'} | 'current' | 'lineRating' | 0 - 10|
% * ------ |0 - 10| - numeric value for the fixed line width
% * ------ |'numPhases'| - thicker lines for 3-phase power lines
% * ------ |'current'| - thickness is linearly related to the current flowing through the lines relative to the maximum current in any line
% * ------ |'lineRating'| - thickness is linearly related to the current rating of the line relative to the maximum line rating
% * -- *|'SubstationMarker'|* - Property for if the substation should be marked |{'on'} | 'off'|
% * -- *|'SubEquipmentMarker'|* - Property for if equipment (such as loads, transformers, etc.) in the substation (using distance) whose marker is turned on should should be marked |'on' | {'off'}|
% * -- *|'PVMarker'|* - Property for if the PV PCC should be marked (if it exists) |{'on'} | 'off'|
% * -- *|'GeneratorMarker'|* - Property for if generators should be marked (if it exists) |{'on'} | 'off'|
% * -- *|'LoadMarker'|* - Property for if loads should be marked |{'on'} | 'off'|
% * -- *|'RegulatorMarker'|* - Property for if controlled transfomers such as regulators (LTC and VREG) should be marked |{'on'} | 'off'|
% * -- *|'MVTransformerMarker'|* - Property for if medium-voltage transfomers (>1000V) should be marked |'on' | {'off'}|
% * -- *|'BoosterMarker'|* - Property for if boosters transformers (uncontrolled NLTC) should be marked |{'on'} | 'off'|
% * -- *|'ServiceTransformerMarker'|* - Property for if service transfomers (<1000V) should be marked |'on' | {'off'}|
% * -- *|'CapacitorMarker'|* - Property for if capacitors should be marked |{'on'} | 'off'|
% * -- *|'CapacitorLabel'|* - Property for if capacitors should be labeled with textbox/arrow for capacitor size |'on' | {'off'}|
% * -- *|'EndOfFeederMarker'|* - Property for if the end of the feeder by distance (3-phase section) should be marked |'on' | {'off'}|
% * -- *|'EndOfFeederLabel'|* - Property for if the end of the feeder should be labeled with textbox/arrow for capacitor size |'on' | {'off'}|
% * -- *|'CustomMarker'|* - Property for marking a custom bus by the user specifying a bus name or a cell array of bus names |{'off'} | busNameString|
% * -- *|'CustomLegend'|* - Text to place in the legend describing the custom bus specified in CustomMarker
% * -- *|'EnergyMeter'|* - Name or cell array of names of the energy meter zones to plot |{'all'} | energyMeterName|
% * -- *|'NumPhases'|* - Property for if only lines with the specified number of phases should be plotted |[1,2,3] | 1 | [2,3] | [1,2]|
% * -- *|'PhasesToPlot'|* - Property for which phases to plot (A,B,C). True/False values for each phase |[1,1,1] | [1,0,0]|
% * -- *|'BusName'|* - Property for the name of the bus (string) that the circuit should be plotted to.  Only the direct line between the bus and the substation will be plotted, unless all buses are selected. |{'all'} | busName|
% * -- *|'Downstream'|* - If a BusName is given, all buses in the electrical path to the substation (upstream) will be plotted, and if this property is on, all buses in the electrical path downstream of BusName will be plotted too |'on' | {'off'}|
% * -- *|'MappingBackground'|* - Property for if the satellite image should be displayed in the background. Note, this only works if the coordinates are in latitude/longitude values or if initCoordConversion was performed. |{'none'} | 'hybrid' | 'satellite' | 'roadmap' | 'terrain'|
% * -- *|'Lines'|* - Structure of the circuit lines from getLineInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'PV'|* - Structure of the PV from getPVInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Generators'|* - Structure of the Generators from getGeneratorInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Transformers'|* - Structure of the circuit transformers from getTransformerInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Capacitors'|* - Structure of the circuit capacitors from getCapacitorInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Loads'|* - Structure of the circuit loads from getLoadInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
%
%% Outputs
% * *|Handles|* - structure of handles for each type of object plotted in the figure
% * A figure of the circuit is displayed in the current axes based on the option inputs
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
% Examples of several different circuit plots that can be created
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% figure; plotCircuitLines(DSSCircObj,'CapacitorMarker','on','LoadMarker','on')
% figure; plotCircuitLines(DSSCircObj,'Coloring','perPhase','Thickness',3,'MappingBackground','hybrid')
% figure; plotCircuitLines(DSSCircObj,'Coloring','voltagePU','EndOfFeederMarker','on')
% figure; plotCircuitLines(DSSCircObj,'Coloring','resistance')
% figure; plotCircuitLines(DSSCircObj,'Coloring','faultCurrent1P')
% figure; plotCircuitLines(DSSCircObj,'Coloring','lineLoading')
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_PV_Central_7_5.dss"'];
% DSSText.command = 'solve';
% figure; plotCircuitLines(DSSCircObj,'Coloring','voltage120')
%

function Handles = plotCircuitLines(DSSCircObj, varargin)
global Coloring ContourScale Thickness SubstationMarker SubEquipmentMarker PVMarker GeneratorMarker LoadMarker CapacitorMarker CapacitorLabel RegulatorMarker MVTransformerMarker BoosterMarker ServiceTransformerMarker EndOfFeederMarker EndOfFeederLabel CustomMarker CustomLegend EnergyMeter NumPhases PhasesToPlot MappingBackground BusName Downstream

%% Parse Inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);

p.addParamValue('Coloring', 'numPhases', @(x)any(strcmp(x,{'numPhases','perPhase','energyMeter','voltage120','voltage','voltagePU','voltage120LL','voltageLL','voltagePULL','lineLoading','realLosses','reactiveLosses','distance','unbalance','voltageAngle','powerFactor','powerFlowDirection','impedance','resistance','reactance','faultCurrent3P','faultCurrent1P','faultCurrentLL'}))|| (isnumeric(x) && length(x)==3 && all(x<=1) && all(x>=0)));
p.addParamValue('ContourScale', 'auto', @(x)strcmp(x,'auto')|| (isnumeric(x) && length(x)==2 && x(2)>x(1)));
p.addParamValue('Thickness', 'numPhases', @(x)any(strcmp(x,{'numPhases','current','lineRating'}))|| (isnumeric(x) && x>0 && x<=10));

p.addParamValue('SubstationMarker', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('SubEquipmentMarker', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('PVMarker', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('GeneratorMarker', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('LoadMarker', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('RegulatorMarker', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('MVTransformerMarker', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('BoosterMarker', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('ServiceTransformerMarker', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('CapacitorMarker', 'on', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('CapacitorLabel', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('EndOfFeederMarker', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('EndOfFeederLabel', 'off', @(x)any(strcmp(x,{'on','off'})));
p.addParamValue('CustomMarker', 'off', @(x)isstr(x)||iscellstr(x));
p.addParamValue('CustomLegend', '', @(x)isstr(x)||iscellstr(x));

p.addParamValue('EnergyMeter', 'all', @(x)ischar(x)||iscell(x));
p.addParamValue('NumPhases', [1,2,3], @(x)isnumeric(x) && all(x>=1) && all(x<=3));
p.addParamValue('PhasesToPlot', [1,1,1], @(x)isnumeric(x) && length(x)==3 && all(x<=1));
p.addParamValue('BusName', 'all', @ischar);
p.addParamValue('Downstream', 'off', @(x)any(strcmp(x,{'on','off'})));

p.addParamValue('MappingBackground', 'none', @(x)any(strcmp(x,{'none','hybrid','satellite','roadmap','terrain'})));

p.addParamValue('Lines', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('Generators', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('PV', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('Transformers', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('Capacitors', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('Loads', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));


p.parse(DSSCircObj, varargin{:}); %parse inputs
nargoutchk(0,1);

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;

%% If there are no input parameters, show plotCircuitLinesOptions GUI to allow the user to select plotting styles
if isempty(varargin)
    lineHandles = plotCircuitLinesOptions(DSSCircObj,{'GiveOptions'});
    uiwait(lineHandles);
    drawnow;
end


%% If Lines structure wasn't send in, get it from the COM
if ischar(Lines) && strcmp(Lines,'noInput')
    Lines = getLineInfo(DSSCircObj);
end
LinesToPlot = Lines;

% Check Lines Structure
if isempty(LinesToPlot)
    error('The Lines structure is empty. There is nothing to plot.');
end


%% If specific energy meter was selected, filter to that
if ischar(EnergyMeter) %if only one name sent in as a string, convert to cell
    EnergyMeter = {EnergyMeter};
end
if ~strcmp('all',EnergyMeter)
    %Check that all energy meter names exist in the circuit
    allMeterNames = lower(DSSCircuit.Meters.AllNames);
    if any(~ismember(lower(EnergyMeter),allMeterNames))
        error(sprintf('The EnergyMeter %s not found in the circuit',EnergyMeter{~ismember(lower(EnergyMeter),allMeterNames)}));
    end
    
    %Loop through each energy meter name sent in to get all lines to plot
    allZoneBranches = [];
    for ii=1:length(EnergyMeter)
        DSSCircuit.Meters.Name = EnergyMeter{ii};
        allZoneBranches = [allZoneBranches; DSSCircuit.Meters.AllBranchesInZone];
    end
    
    %Only keep the lines in the selected energy meter zones
    allZoneLines = allZoneBranches(strncmpi('line',allZoneBranches,4));
    allZoneLines = cellfun(@(x)x(6:end),allZoneLines,'UniformOutput',false);
    LinesToPlot = LinesToPlot(ismember(lower({LinesToPlot.name}),lower(allZoneLines)));
    if isempty(LinesToPlot)
        error('There are no lines to plot that meet the selected criteria.');
    end
else
    EnergyMeter = DSSCircuit.Meters.AllNames;
end

%% Filter lines that have the right number of phases
condition = ismember([LinesToPlot.numPhases],NumPhases);
LinesToPlot = LinesToPlot(condition);
if isempty(LinesToPlot)
    error('There are no lines to plot that meet the selected criteria.');
end


%% Filter lines that do not have the correct phase (A,B,C)
phaseA = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.1)','Tokens'));
phaseB = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.2)','Tokens'));
phaseC = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.3)','Tokens'));
isThreePhase = sum([phaseA;phaseB;phaseC])==0;
phaseA(isThreePhase) = 1;
phaseB(isThreePhase) = 1;
phaseC(isThreePhase) = 1;
LinesToPlot = LinesToPlot(phaseA*PhasesToPlot(1) | phaseB*PhasesToPlot(2) | phaseC*PhasesToPlot(3));
if isempty(LinesToPlot)
    error('There are no lines to plot that meet the selected criteria.');
end


%% Remove all lines except between BusName and substation (if BusName is set)
if ~strcmp(BusName,'all')
    Transformers = getTransformerInfo(DSSCircObj);
    selectedBuses = findUpstreamBuses(DSSCircObj,BusName,'Lines',Lines,'Transformers',Transformers);
    if strcmp(Downstream,'on') %If Downstream property is on, also plot all buses downstream of BusName
        downstreamBuses = findDownstreamBuses(DSSCircObj,BusName,'Lines',Lines,'Transformers',Transformers);
        selectedBuses = [selectedBuses,downstreamBuses'];
    end
    linesBus1 = regexprep({LinesToPlot.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    linesBus2 = regexprep({LinesToPlot.bus2},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    LinesToPlot = LinesToPlot(ismember(linesBus1,selectedBuses) & ismember(linesBus2,selectedBuses));
    if isempty(LinesToPlot)
        error('There are no lines to plot that meet the selected criteria.');
    end
end


%% Get List of Bus Names that will be in the circuit plot (everything will be filtered on this)
allBusesToPlot = [{LinesToPlot.bus1},{LinesToPlot.bus2}];
allBusesToPlot = lower(regexprep(allBusesToPlot,'(\.[0-9]+)','')); %take out the phase numbers on buses if they have them


%% Remove all lines that do not have coordinates
condition = [LinesToPlot.bus1CoordDefined] & [LinesToPlot.bus2CoordDefined];%condition = sum(reshape([LinesToPlot.bus1Coordinates],2,[]))~=0 & sum(reshape([LinesToPlot.bus2Coordinates],2,[]))~=0;
LinesToPlot = LinesToPlot(condition);

if isempty(LinesToPlot)
   error('Unable to find any buses with coordinates. Make sure Buscoords is defined in OpenDSS.'); 
end
bus1Coord = reshape([LinesToPlot.bus1Coordinates],2,[])';
bus2Coord = reshape([LinesToPlot.bus2Coordinates],2,[])';


%% plot lines
lineHandles = plot([bus1Coord(:,2),bus2Coord(:,2)]',[bus1Coord(:,1),bus2Coord(:,1)]');
set(lineHandles,{'DisplayName'},strcat('line.',{LinesToPlot.name})');
cbar_axes = colorbar;%('Peer',gca,'Location','East');
set(gca,'FontWeight','bold','FontSize',10);
hold on;
FigHandles.lines = lineHandles;

legendHandles = []; %empty variabes to store legend information
legendText = {};


%% Plot Line Coloring (Based on input)
colormap('jet');
cmap = colormap;

if ischar(Coloring)
    switch Coloring
        
        case {'numPhases'}
            % number of phases
            colorMatrix = zeros(length(lineHandles),3);
            oneOrTwoPhase = [LinesToPlot.numPhases]==1 | [LinesToPlot.numPhases]==2;
            colorMatrix(oneOrTwoPhase,:) = 0.5*ones(sum(oneOrTwoPhase),3);
            for ii=1:length(lineHandles)
                set(lineHandles(ii),'Color',colorMatrix(ii,:));
            end
            colorbar('off');
            
        case {'perPhase'}
            % phase number
            phaseA = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.1)','Tokens'));
            phaseB = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.2)','Tokens'));
            phaseC = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.3)','Tokens'));
            colorMatrix = [phaseA' phaseB' phaseC'];
            threePhase = [LinesToPlot.numPhases]==3;
            colorMatrix(threePhase,:) = zeros(sum(threePhase),3); %return 3 phase back to black
            differentColors = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
            colorHandles = zeros(6,1);
            colorLegends = {'Phase A','Phase B','Phase C','Phase AB','Phase AC','Phase BC'};
            for ii=1:length(lineHandles)
                set(lineHandles(ii),'Color',colorMatrix(ii,:));
                [tf, loc] = ismember(colorMatrix(ii,:),differentColors,'rows');
                if tf
                    colorHandles(loc) = lineHandles(ii);
                end
            end
            colorLegends = colorLegends(colorHandles~=0);
            colorHandles = colorHandles(colorHandles~=0);
            legendHandles = [legendHandles,colorHandles'];
            legendText = [legendText,colorLegends];
%             [legend_h,object_h,plot_h,text_strings] = legend;
%             legend([plot_h;colorHandles],[text_strings,colorLegends]);
            colorbar('off');
            %text = annotation('textbox','String',sprintf('Red - Phase A\nGreen - Phase B\nBlue - Phase C\nYellow - Phase AB\nPurple - Phase AC\nCyan - Phase BC'),'Position',[0.2 0.82 0.1 0.1],'HorizontalAlignment','center','FontWeight','bold','FitBoxToText','on','Margin',1,'VerticalAlignment','middle');
            title('Circuit Plot by Phase','FontWeight','bold','FontSize',12);
            
        case {'energyMeter'}
            % energy meter zones
            allMeterNames = EnergyMeter; %DSSCircuit.Meters.AllNames;
            zoneNumber = zeros(length(LinesToPlot),1);
            for ii=1:length(allMeterNames) %Loop through each energy meter name to get zone number for each line
                DSSCircuit.Meters.Name = allMeterNames{ii};
                allZoneBranches = DSSCircuit.Meters.AllBranchesInZone;
                allZoneLines = allZoneBranches(strncmpi('line',allZoneBranches,4));
                allZoneLines = cellfun(@(x)x(6:end),allZoneLines,'UniformOutput',false);
                TFinZone = ismember(lower({LinesToPlot.name}),lower(allZoneLines));
                if any(TFinZone)
                    zoneNumber(TFinZone) = TFinZone(TFinZone)*ii;
                    legendHandles = [legendHandles,lineHandles(find(TFinZone,1,'first'))];
                    legendText = [legendText,allMeterNames(ii)];
                end
            end
            colormap('HSV');
            cmap = colormap;
            for ii=1:length(lineHandles)
                lineColors = round(63*(zoneNumber(ii)-1)/length(allMeterNames)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end

            colorbar('off');
            title('Circuit Plot by EnergyMeter','FontWeight','bold','FontSize',12);
            
        case {'voltage120'}
            % voltage
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minVoltage = min([[LinesToPlot.bus1VoltagePU],[LinesToPlot.bus2VoltagePU]]);
                maxVoltage = max([[LinesToPlot.bus1VoltagePU],[LinesToPlot.bus2VoltagePU]]);
            else
                minVoltage = ContourScale(1)/120; %they are told to put it in on 120V base scale
                maxVoltage = ContourScale(2)/120;
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1VoltagePU-minVoltage)/(maxVoltage-minVoltage)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minVoltage*120 maxVoltage*120]);
            title('Voltage (120 V Base)','FontWeight','bold','FontSize',12);
            
        case {'voltage'}
            % voltage
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minVoltage = min([[LinesToPlot.bus1Voltage],[LinesToPlot.bus2Voltage]]);
                maxVoltage = max([[LinesToPlot.bus1Voltage],[LinesToPlot.bus2Voltage]]);
            else
                minVoltage = ContourScale(1)*1000; %they are told to put it in kV
                maxVoltage = ContourScale(2)*1000;
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1Voltage-minVoltage)/(maxVoltage-minVoltage)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minVoltage/1000 maxVoltage/1000]);
            title('Voltage (kV)','FontWeight','bold','FontSize',12);
            
        case {'voltagePU'}
            % voltage
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minVoltage = min([[LinesToPlot.bus1VoltagePU],[LinesToPlot.bus2VoltagePU]]);
                maxVoltage = max([[LinesToPlot.bus1VoltagePU],[LinesToPlot.bus2VoltagePU]]);
            else
                minVoltage = ContourScale(1); %they are told to put it in pu
                maxVoltage = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1VoltagePU-minVoltage)/(maxVoltage-minVoltage)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minVoltage maxVoltage]);
            title('Voltage (pu)','FontWeight','bold','FontSize',12);
            
        case {'voltage120LL'}
            % voltage
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minVoltage = min([[LinesToPlot.bus1VoltageLLPU],[LinesToPlot.bus2VoltageLLPU]]);
                maxVoltage = max([[LinesToPlot.bus1VoltageLLPU],[LinesToPlot.bus2VoltageLLPU]]);
            else
                minVoltage = ContourScale(1)/120; %they are told to put it in on 120V base scale
                maxVoltage = ContourScale(2)/120;
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1VoltagePU-minVoltage)/(maxVoltage-minVoltage)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minVoltage*120 maxVoltage*120]);
            title('Voltage L-L (120 V Base)','FontWeight','bold','FontSize',12);
            
        case {'voltageLL'}
            % voltage
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minVoltage = min([[LinesToPlot.bus1VoltageLL],[LinesToPlot.bus2VoltageLL]]);
                maxVoltage = max([[LinesToPlot.bus1VoltageLL],[LinesToPlot.bus2VoltageLL]]);
            else
                minVoltage = ContourScale(1)*1000; %they are told to put it in kV
                maxVoltage = ContourScale(2)*1000;
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1Voltage-minVoltage)/(maxVoltage-minVoltage)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minVoltage/1000 maxVoltage/1000]);
            title('Voltage L-L (kV)','FontWeight','bold','FontSize',12);
            
        case {'voltagePULL'}
            % voltage
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minVoltage = min([[LinesToPlot.bus1VoltageLLPU],[LinesToPlot.bus2VoltageLLPU]]);
                maxVoltage = max([[LinesToPlot.bus1VoltageLLPU],[LinesToPlot.bus2VoltageLLPU]]);
            else
                minVoltage = ContourScale(1); %they are told to put it in pu
                maxVoltage = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1VoltagePU-minVoltage)/(maxVoltage-minVoltage)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minVoltage maxVoltage]);
            title('Voltage L-L (pu)','FontWeight','bold','FontSize',12);
            
        case {'lineLoading'}
            % line loading
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minLineLoading = min([LinesToPlot.bus1Current]./[LinesToPlot.lineRating]);
                maxLineLoading = max([LinesToPlot.bus1Current]./[LinesToPlot.lineRating]);
            else
                minLineLoading = ContourScale(1)/100; %they are told to put it in percent
                maxLineLoading = ContourScale(2)/100;
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1Current/LinesToPlot(ii).lineRating-minLineLoading)/(maxLineLoading-minLineLoading)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minLineLoading*100 maxLineLoading*100]);
            title('Line Loading Percent (100*Current/Rating)','FontWeight','bold','FontSize',12);
            
        case {'realLosses'}
            % line losses
            lineLength = abs([LinesToPlot.bus2Distance]-[LinesToPlot.bus1Distance]);
            realLosses = real([LinesToPlot.losses])./lineLength;
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minLineLosses = min(realLosses(~isinf(realLosses)));
                maxLineLosses = max(realLosses(~isinf(realLosses)));
            else
                minLineLosses = ContourScale(1);
                maxLineLosses = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(realLosses(ii)-minLineLosses)/(maxLineLosses-minLineLosses)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minLineLosses maxLineLosses]);
            title('Real Line Losses (kW/km)','FontWeight','bold','FontSize',12);
            
        case {'reactiveLosses'}
            % line losses
            lineLength = abs([LinesToPlot.bus2Distance]-[LinesToPlot.bus1Distance]);
            reactiveLosses = imag([LinesToPlot.losses])./lineLength;
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minLineLosses = min(reactiveLosses(~isinf(reactiveLosses)));
                maxLineLosses = max(reactiveLosses(~isinf(reactiveLosses)));
            else
                minLineLosses = ContourScale(1);
                maxLineLosses = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(reactiveLosses(ii)-minLineLosses)/(maxLineLosses-minLineLosses)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minLineLosses maxLineLosses]);
            title('Reactive Line Losses (kVAR/km)','FontWeight','bold','FontSize',12);
            
        case {'distance'}
            % distance from substation
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minDistance = min([LinesToPlot.bus1Distance]);
                maxDistance = max([LinesToPlot.bus1Distance]);
            else
                minDistance = ContourScale(1);
                maxDistance = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1Distance-minDistance)/(maxDistance-minDistance)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minDistance maxDistance]);
            title('Distance from Substation (km)','FontWeight','bold','FontSize',12);
            
        case {'unbalance'}
            % unbalance of phases
            for ii=1:length(LinesToPlot)
                if LinesToPlot(ii).numPhases ~=1
                    phases = [~isempty(regexp(LinesToPlot(ii).bus1,'.^?(\.1)','Tokens')),~isempty(regexp(LinesToPlot(ii).bus1,'.^?(\.2)','Tokens')),~isempty(regexp(LinesToPlot(ii).bus1,'.^?(\.3)','Tokens'))];
                    if ~any(phases ~= 0) %If all phases are blank (ie they were not put into OpenDSS), then default to all 3 phases
                        phases = [true true true];
                    end
                    power = sqrt(LinesToPlot(ii).bus1PhasePowerReal.^2 + LinesToPlot(ii).bus1PhasePowerReactive.^2);
                    power = power(phases);
                    powerUnbalance(ii) = abs( max(power)-min(power) ); %maximum difference in magnitude of power flow between phases
                else
                    powerUnbalance(ii) = 0;
                end
            end
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minUnbalance = min(powerUnbalance);
                maxUnbalance = max(powerUnbalance);
            else
                minUnbalance = ContourScale(1);
                maxUnbalance = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(powerUnbalance(ii)-minUnbalance)/(maxUnbalance-minUnbalance)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minUnbalance maxUnbalance]);
            title('Power (kVA) Unbalance','FontWeight','bold','FontSize',12);
            
        case {'voltageAngle'}
            % average voltage angle
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minAngle = min([LinesToPlot.bus1VoltageAngle]);
                maxAngle = max([LinesToPlot.bus1VoltageAngle]);
            else
                minAngle = ContourScale(1)/180*pi;
                maxAngle = ContourScale(2)/180*pi;
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(LinesToPlot(ii).bus1VoltageAngle-minAngle)/(maxAngle-minAngle)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minAngle*180/pi maxAngle*180/pi]);
            title('Voltage Angle (degrees)','FontWeight','bold','FontSize',12);
            
        case {'powerFactor'}
            % power factor
            totalPower = sqrt([LinesToPlot.bus1PowerReal].^2+[LinesToPlot.bus1PowerReactive].^2);
            powerFactor = ([LinesToPlot.bus1PowerReal])./totalPower.*sign([LinesToPlot.bus1PowerReactive]);
            powerFactor(powerFactor<0) = powerFactor(powerFactor<0)+2; %move negative power factor to above 1
            powerFactor = 2-powerFactor; %flip axis
            %initial extreme power factors
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minPF = 0.5;
                maxPF = 1.5;
            else
                minPF = -ContourScale(1);
                maxPF = 2-ContourScale(2);
            end
            
            powerFactor(powerFactor<minPF | powerFactor>maxPF) = 0; %set extreme power factors to zero
            if isstr(ContourScale) && strcmp(ContourScale,'auto') %recalculate max/min for auto scale
                minPF = min(powerFactor(powerFactor~=0));
                maxPF = max(powerFactor(powerFactor~=0));
            end
            for ii=1:length(lineHandles)
                if powerFactor(ii)~=0
                    lineColors = round(63*(powerFactor(ii)-minPF)/(maxPF-minPF)+1);
                    set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
                else
                    set(lineHandles(ii),'Color',[0.1 0.1 0.1]);
                end
            end
            caxis(gca,[minPF maxPF]);
            title('Power Factor','FontWeight','bold','FontSize',12);
            colorLabels = str2num(get(cbar_axes,'YTickLabel'));
            colorLabels(colorLabels<1) = -colorLabels(colorLabels<1);
            colorLabels(colorLabels>1) = 2-colorLabels(colorLabels>1);
            set(cbar_axes,'YTickLabel',num2str(colorLabels));
            
        case {'powerFlowDirection'}
            bus1Power = [LinesToPlot.bus1PowerReal];
            bus2Power = [LinesToPlot.bus2PowerReal];
            powerDownstream = bus1Power;
            powerDownstream([LinesToPlot.bus1Distance]>[LinesToPlot.bus2Distance]) = bus2Power([LinesToPlot.bus1Distance]>[LinesToPlot.bus2Distance]);
            maxPowerDownstream = max(abs(powerDownstream));
            
            cmap = [[((96.5:0.5:128)/128)';((83:-1:20)/128)'],[((20:83)/128)';((128:-0.5:96.5)/128)'],[((20:83)/128)';((83:-1:20)/128)']];
            colormap(cmap);
            
            for ii=1:length(lineHandles)
                if abs(powerDownstream(ii))<0.01 %small flow less than 0.01 kW
                    set(lineHandles(ii),'Color',[0.7,0.7,0.7]);
                elseif powerDownstream(ii)<0
                    set(lineHandles(ii),'Color',cmap(round(64+powerDownstream(ii)/maxPowerDownstream*63),:));
                else
                    set(lineHandles(ii),'Color',cmap(round(65+powerDownstream(ii)/maxPowerDownstream*63),:));
                end
            end
            caxis(gca,[-maxPowerDownstream maxPowerDownstream]);
            title('Power Flow Direction (kW)','FontWeight','bold','FontSize',12);
            
            [B,IX] = sort(powerDownstream);
            downstreamPowerGroup = hggroup; set(lineHandles(IX(B>=0)),'Parent',downstreamPowerGroup); %group together for one item in the legend
            legendHandles = [legendHandles,downstreamPowerGroup];
            legendText = [legendText,'Downstream Flow'];
            reversePowerGroup = hggroup; set(flipud(lineHandles(IX(B<0))),'Parent',reversePowerGroup); %group together for one item in the legend
            legendHandles = [legendHandles,reversePowerGroup];
            legendText = [legendText,'Reverse Flow'];
            legend(legendHandles,legendText,'Interpreter','none');
            figureChildren = get(gca,'children');
            set(gca,'children',[figureChildren(3:end);figureChildren(1:2)]); %grouping moved the lines to the top, so move them back to the bottom
    
        case {'impedance'}
            if all([LinesToPlot.bus1Zsc1]==0) %if there is no impedance data, a faultstudy will need to be run to get the results.
                button = questdlg('Plotting is about to modify the most recent OpenDSS solution and change the mode to faultstudy.','OpenDSS Solution Mode','Ok','Cancel','Ok');
                if strcmp(button,'Ok')
                    DSSText.Command = 'Disable Load.*'; % Disable all loads
                    DSSText.Command = 'Disable Capacitor.*'; % Disable all caps
                    DSSText.command = 'solve mode=faultstudy';
                    DSSText.Command = 'Enable Load.*'; % Enable all loads
                    DSSText.Command = 'Enable Capacitor.*'; % Enable all caps
                    LinesToPlot = getLineInfo(DSSCircObj,{LinesToPlot.name});
                else
                    return;
                end
            end
            
            bus1Zsc1 = vertcat(LinesToPlot.bus1Zsc1);
            bus1Zsc1 = hypot(bus1Zsc1(:,1),bus1Zsc1(:,2));
            bus2Zsc1 = vertcat(LinesToPlot.bus2Zsc1);
            bus2Zsc1 = hypot(bus2Zsc1(:,1),bus2Zsc1(:,2));
            
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minImpedance = min([bus1Zsc1;bus2Zsc1]);
                maxImpedance = max([bus1Zsc1;bus2Zsc1]);
            else
                minImpedance = ContourScale(1);
                maxImpedance = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(mean([bus1Zsc1(ii),bus2Zsc1(ii)])-minImpedance)/(maxImpedance-minImpedance)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minImpedance maxImpedance]);
            title('Positive-Sequence Short-Circuit Impedance Magnitude','FontWeight','bold','FontSize',12);
            
        case {'resistance'}
            if all([LinesToPlot.bus1Zsc1]==0) %if there is no impedance data, a faultstudy will need to be run to get the results.
                button = questdlg('Plotting is about to modify the most recent OpenDSS solution and change the mode to faultstudy.','OpenDSS Solution Mode','Ok','Cancel','Ok');
                if strcmp(button,'Ok')
                    DSSText.Command = 'Disable Load.*'; % Disable all loads
                    DSSText.Command = 'Disable Capacitor.*'; % Disable all caps
                    DSSText.command = 'solve mode=faultstudy';
                    DSSText.Command = 'Enable Load.*'; % Enable all loads
                    DSSText.Command = 'Enable Capacitor.*'; % Enable all caps
                    LinesToPlot = getLineInfo(DSSCircObj,{LinesToPlot.name});
                else
                    return;
                end
            end
            
            bus1Zsc1 = vertcat(LinesToPlot.bus1Zsc1);
            bus1Zsc1 = bus1Zsc1(:,1);
            bus2Zsc1 = vertcat(LinesToPlot.bus2Zsc1);
            bus2Zsc1 = bus2Zsc1(:,1);
            
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minResistance = min([bus1Zsc1;bus2Zsc1]);
                maxResistance = max([bus1Zsc1;bus2Zsc1]);
            else
                minResistance = ContourScale(1);
                maxResistance = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(mean([bus1Zsc1(ii),bus2Zsc1(ii)])-minResistance)/(maxResistance-minResistance)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minResistance maxResistance]);
            title('Positive-Sequence Short-Circuit Resistance','FontWeight','bold','FontSize',12);
            
        case {'reactance'}
            if all([LinesToPlot.bus1Zsc1]==0) %if there is no impedance data, a faultstudy will need to be run to get the results.
                button = questdlg('Plotting is about to modify the most recent OpenDSS solution and change the mode to faultstudy.','OpenDSS Solution Mode','Ok','Cancel','Ok');
                if strcmp(button,'Ok')
                    DSSText.Command = 'Disable Load.*'; % Disable all loads
                    DSSText.Command = 'Disable Capacitor.*'; % Disable all caps
                    DSSText.command = 'solve mode=faultstudy';
                    DSSText.Command = 'Enable Load.*'; % Enable all loads
                    DSSText.Command = 'Enable Capacitor.*'; % Enable all caps
                    LinesToPlot = getLineInfo(DSSCircObj,{LinesToPlot.name});
                else
                    return;
                end
            end
            
            bus1Zsc1 = vertcat(LinesToPlot.bus1Zsc1);
            bus1Zsc1 = bus1Zsc1(:,2);
            bus2Zsc1 = vertcat(LinesToPlot.bus2Zsc1);
            bus2Zsc1 = bus2Zsc1(:,2);
            
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minReactance = min([bus1Zsc1;bus2Zsc1]);
                maxReactance = max([bus1Zsc1;bus2Zsc1]);
            else
                minReactance = ContourScale(1);
                maxReactance = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(mean([bus1Zsc1(ii),bus2Zsc1(ii)])-minReactance)/(maxReactance-minReactance)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minReactance maxReactance]);
            title('Positive-Sequence Short-Circuit Reactance','FontWeight','bold','FontSize',12);
            
        case {'faultCurrent3P'}
            if all([LinesToPlot.bus1Zsc1]==0) %if there is no impedance data, a faultstudy will need to be run to get the results.
                button = questdlg('Plotting is about to modify the most recent OpenDSS solution and change the mode to faultstudy.','OpenDSS Solution Mode','Ok','Cancel','Ok');
                if strcmp(button,'Ok')
                    DSSText.Command = 'Disable Load.*'; % Disable all loads
                    DSSText.Command = 'Disable Capacitor.*'; % Disable all caps
                    DSSText.command = 'solve mode=faultstudy';
                    DSSText.Command = 'Enable Load.*'; % Enable all loads
                    DSSText.Command = 'Enable Capacitor.*'; % Enable all caps
                    LinesToPlot = getLineInfo(DSSCircObj,{LinesToPlot.name});
                else
                    return;
                end
            end
            
            DSSText.command = 'export faultstudy';
            raw = importdata(DSSText.result);
            faultBusNames = raw.textdata(2:end,1);
%             faultBusNames = cellfun(@num2str,faultBusNames,'UniformOutput',false);
            faultBusNames = regexprep(faultBusNames,' ','');
            faultCurrents = raw.data(:,:);
            
            [tf loc] = ismember(lower(regexprep({LinesToPlot.bus1},'(\.[0-9]+)','')),lower(faultBusNames));
            bus1FaultCurrents = faultCurrents(loc,1);
            bus1FaultCurrents(faultCurrents(loc,3)==0) = 0; %if the line to line fault current is zero, then obviously there is no 3-phase fault
            [tf loc] = ismember(lower(regexprep({LinesToPlot.bus2},'(\.[0-9]+)','')),lower(faultBusNames));
            bus2FaultCurrents = faultCurrents(loc,1);
            bus2FaultCurrents(faultCurrents(loc,3)==0) = 0; %if the line to line fault current is zero, then obviously there is no 3-phase fault

            
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minFaultCurrent = min([bus1FaultCurrents(bus1FaultCurrents~=0);bus2FaultCurrents(bus2FaultCurrents~=0)]);
                maxFaultCurrent = max([bus1FaultCurrents;bus2FaultCurrents]);
            else
                minFaultCurrent = ContourScale(1);
                maxFaultCurrent = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                faultCurrent = mean([bus1FaultCurrents(ii),bus2FaultCurrents(ii)]);
                if faultCurrent==0
                    set(lineHandles(ii),'Color',[0,0,0]);
                else
                    lineColors = round(63*(faultCurrent-minFaultCurrent)/(maxFaultCurrent-minFaultCurrent)+1);
                    set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
                end
            end
            caxis(gca,[minFaultCurrent maxFaultCurrent]);
            title('3-Phase Fault Current (A)','FontWeight','bold','FontSize',12);
            
        case {'faultCurrent1P'}
            if all([LinesToPlot.bus1Zsc1]==0) %if there is no impedance data, a faultstudy will need to be run to get the results.
                button = questdlg('Plotting is about to modify the most recent OpenDSS solution and change the mode to faultstudy.','OpenDSS Solution Mode','Ok','Cancel','Ok');
                if strcmp(button,'Ok')
                    DSSText.Command = 'Disable Load.*'; % Disable all loads
                    DSSText.Command = 'Disable Capacitor.*'; % Disable all caps
                    DSSText.command = 'solve mode=faultstudy';
                    DSSText.Command = 'Enable Load.*'; % Enable all loads
                    DSSText.Command = 'Enable Capacitor.*'; % Enable all caps
                    LinesToPlot = getLineInfo(DSSCircObj,{LinesToPlot.name});
                else
                    return;
                end
            end
            
            DSSText.command = 'export faultstudy';
            raw = importdata(DSSText.result);
            faultBusNames = raw.textdata(2:end,1);
%             faultBusNames = cellfun(@num2str,faultBusNames,'UniformOutput',false);
            faultBusNames = regexprep(faultBusNames,' ','');
            faultCurrents = raw.data(:,:);
            
            [tf loc] = ismember(lower(regexprep({LinesToPlot.bus1},'(\.[0-9]+)','')),lower(faultBusNames));
            bus1FaultCurrents = faultCurrents(loc,2);
            [tf loc] = ismember(lower(regexprep({LinesToPlot.bus2},'(\.[0-9]+)','')),lower(faultBusNames));
            bus2FaultCurrents = faultCurrents(loc,2);
            
            
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minFaultCurrent = min([bus1FaultCurrents;bus2FaultCurrents]);
                maxFaultCurrent = max([bus1FaultCurrents;bus2FaultCurrents]);
            else
                minFaultCurrent = ContourScale(1);
                maxFaultCurrent = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                lineColors = round(63*(mean([bus1FaultCurrents(ii),bus2FaultCurrents(ii)])-minFaultCurrent)/(maxFaultCurrent-minFaultCurrent)+1);
                set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
            end
            caxis(gca,[minFaultCurrent maxFaultCurrent]);
            title('1-Phase Fault Current (A)','FontWeight','bold','FontSize',12);
            
        case {'faultCurrentLL'}
            if all([LinesToPlot.bus1Zsc1]==0) %if there is no impedance data, a faultstudy will need to be run to get the results.
                button = questdlg('Plotting is about to modify the most recent OpenDSS solution and change the mode to faultstudy.','OpenDSS Solution Mode','Ok','Cancel','Ok');
                if strcmp(button,'Ok')
                    DSSText.Command = 'Disable Load.*'; % Disable all loads
                    DSSText.Command = 'Disable Capacitor.*'; % Disable all caps
                    DSSText.command = 'solve mode=faultstudy';
                    DSSText.Command = 'Enable Load.*'; % Enable all loads
                    DSSText.Command = 'Enable Capacitor.*'; % Enable all caps
                    LinesToPlot = getLineInfo(DSSCircObj,{LinesToPlot.name});
                else
                    return;
                end
            end
            
            DSSText.command = 'export faultstudy';
            raw = importdata(DSSText.result);
            faultBusNames = raw.textdata(2:end,1);
%             faultBusNames = cellfun(@num2str,faultBusNames,'UniformOutput',false);
            faultBusNames = regexprep(faultBusNames,' ','');
            faultCurrents = raw.data(:,:);
            
            [tf loc] = ismember(lower(regexprep({LinesToPlot.bus1},'(\.[0-9]+)','')),lower(faultBusNames));
            bus1FaultCurrents = faultCurrents(loc,3);
            [tf loc] = ismember(lower(regexprep({LinesToPlot.bus2},'(\.[0-9]+)','')),lower(faultBusNames));
            bus2FaultCurrents = faultCurrents(loc,3);
            
            
            if isstr(ContourScale) && strcmp(ContourScale,'auto')
                minFaultCurrent = min([bus1FaultCurrents(bus1FaultCurrents~=0);bus2FaultCurrents(bus2FaultCurrents~=0)]);
                maxFaultCurrent = max([bus1FaultCurrents;bus2FaultCurrents]);
            else
                minFaultCurrent = ContourScale(1);
                maxFaultCurrent = ContourScale(2);
            end
            for ii=1:length(lineHandles)
                faultCurrent = mean([bus1FaultCurrents(ii),bus2FaultCurrents(ii)]);
                if faultCurrent==0
                    set(lineHandles(ii),'Color',[0,0,0]);
                else
                    lineColors = round(63*(faultCurrent-minFaultCurrent)/(maxFaultCurrent-minFaultCurrent)+1);
                    set(lineHandles(ii),'Color',cmap(min(64,max(1,lineColors)),:));
                end
                
            end
            caxis(gca,[minFaultCurrent maxFaultCurrent]);
            title('Line-To-Line Fault Current (A)','FontWeight','bold','FontSize',12);
    end
else
    set(lineHandles,'Color',Coloring);
    colorbar('off');
end


%% Plot Line Thickness (based on input)
switch Thickness

    case {'numPhases'}
        % number of phases
        thicknessMatrix = 4*ones(length(lineHandles),1);
        oneOrTwoPhase = [LinesToPlot.numPhases]==1 | [LinesToPlot.numPhases]==2;
        thicknessMatrix(oneOrTwoPhase,:) = 2.5;
        for ii=1:length(lineHandles)
            set(lineHandles(ii),'LineWidth',thicknessMatrix(ii));
        end

    case {'lineRating'}
        % line rating
        maxRating = max([LinesToPlot(abs([LinesToPlot.bus1Distance]-[LinesToPlot.bus2Distance])>0.001 & [LinesToPlot.bus1Distance]>0 & [LinesToPlot.bus2Distance]>0 & [LinesToPlot.bus1Voltage]>1000).lineRating]); %only not substation, not secondary, and cables with length
        for ii=1:length(lineHandles)
            set(lineHandles(ii),'LineWidth',8*LinesToPlot(ii).lineRating/maxRating);
        end
        titleHandle = get(gca,'title');
        titleString = get(titleHandle,'String');
        title([titleString,' with line thickness by line rating'],'FontWeight','bold','FontSize',12);

    case {'current'}
        % line current
        maxCurrent = max([LinesToPlot([LinesToPlot.bus1Distance]>0).bus1Current]); %not in substation
        for ii=1:length(lineHandles)
            set(lineHandles(ii),'LineWidth',10*LinesToPlot(ii).bus1Current/maxCurrent+0.1);
        end
        titleHandle = get(gca,'title');
        titleString = get(titleHandle,'String');
        title([titleString,' with line thickness by line current'],'FontWeight','bold','FontSize',12);

    otherwise
        % Default
        set(lineHandles,'LineWidth',Thickness);
end
        

%% Plot Substation if property SubstationMarker is on
if strcmp(SubstationMarker,'on')
    coordinates = findSubstationLocation(DSSCircObj);
    sh = plot(coordinates(2),coordinates(1),'-kh','MarkerSize',15,'MarkerFaceColor',[0,0.5,1],'LineStyle','none');
    set(sh,'DisplayName','Substation');
    legendHandles = [legendHandles,sh];
    legendText = [legendText,'Substation'];
    FigHandles.substation = sh;
end


%% Plot Generators if property GeneratorMarker is on
if strcmp(GeneratorMarker,'on')
    if ischar(Generators) && strcmp(Generators,'noInput')
        Generators = getGeneratorInfo(DSSCircObj);
    end
    if ~isempty(Generators) && ~strcmp(Generators(1).name,'NONE')
        Generators = Generators(ismember(lower(regexprep({Generators.busName},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(Generators) && ~strcmp(Generators(1).name,'NONE')
        Generators = Generators([Generators.distance]>0.01);
    end
    if ~isempty(Generators) && ~strcmp(Generators(1).name,'NONE')
        
        GeneratorCoords = reshape([Generators.coordinates],2,[])';
        gh = plot(repmat(GeneratorCoords(:,2)',2,1),repmat(GeneratorCoords(:,1)',2,1),'-y*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',[1 0.55 0],'LineStyle','none');
        set(gh,{'DisplayName'},strcat('generator.',{Generators.name})');
        generatorGroup = hggroup; set(gh,'Parent',generatorGroup); %group together for one item in the legend
        legendHandles = [legendHandles,generatorGroup];
        legendText = [legendText,'Generator'];
        FigHandles.generators = gh;
        FigHandles.generatorGroup = generatorGroup;
    else
        FigHandles.generators = '';
    end
end


%% Plot PV plants if property PVMarker is on
if strcmp(PVMarker,'on')
    if ischar(PV) && strcmp(PV,'noInput')
        PV = getPVInfo(DSSCircObj);
    end
    if ~isempty(PV) && ~strcmp(PV(1).name,'NONE')
        PV = PV(ismember(lower(regexprep({PV.busName},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(PV) && ~strcmp(PV(1).name,'NONE')
        PV = PV([PV.distance]>0.01);
    end
    if ~isempty(PV) && ~strcmp(PV(1).name,'NONE')
        PVcoords = reshape([PV.coordinates],2,[])';
        ph = plot(repmat(PVcoords(:,2)',2,1),repmat(PVcoords(:,1)',2,1),'-kp','MarkerSize',11,'MarkerFaceColor','y','LineStyle','none');
        set(ph,{'DisplayName'},strcat('pvsystem.',{PV.name})');
        pvGroup = hggroup; set(ph,'Parent',pvGroup); %group together for one item in the legend
        legendHandles = [legendHandles,pvGroup];
        legendText = [legendText,'PV PCC'];
        FigHandles.PV = ph;
        FigHandles.PVGroup = pvGroup;
    else
        FigHandles.PV = '';
    end
end


%% Plot loads if property LoadMarker is on
if strcmp(LoadMarker,'on')
    if ischar(Loads) && strcmp(Loads,'noInput')
        Loads = getLoadInfo(DSSCircObj);
    end
    if ~isempty(Loads) && ~strcmp(Loads(1).name,'NONE')
        Loads = Loads(ismember(lower(regexprep({Loads.busName},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(Loads) && ~strcmp(Loads(1).name,'NONE')
        Loads = Loads([Loads.distance]>0.01);
    end
    if ~isempty(Loads) && ~strcmp(Loads(1).name,'NONE')
        loadBuses = regexprep({Loads.busName},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        [loadBuses, ia, ic] = unique(loadBuses);
        Loads = Loads(ia);
        Buses = getBusInfo(DSSCircObj,loadBuses,1);
        BusCoords = reshape([Buses.coordinates],2,[])';
        
        lh = plot(repmat(BusCoords(:,2)',2,1),repmat(BusCoords(:,1)',2,1),'-ko','MarkerSize',4,'MarkerFaceColor','w','LineStyle','none');
        set(lh,{'DisplayName'},strcat('load.',{Loads.name})');
        loadGroup = hggroup; set(lh,'Parent',loadGroup); %group together for one item in the legend
        legendHandles = [legendHandles,loadGroup];
        legendText = [legendText,'Loads'];
        FigHandles.loads = lh;
        FigHandles.loadGroup = loadGroup;
    else
        FigHandles.loads = '';
    end
end


%% Plot service transformers if property ServiceTransformerMarker is on
if strcmp(ServiceTransformerMarker,'on')
    if ischar(Transformers) && strcmp(Transformers,'noInput')
        Transformers = getTransformerInfo(DSSCircObj);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers(ismember(lower(regexprep({Transformers.bus1},'(\.[0-9]+)','')),allBusesToPlot) | ismember(lower(regexprep({Transformers.bus2},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers([Transformers.bus1Distance]>0.01);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        if any([Transformers.bus2kV]<1 & ~[Transformers.controlled]) %mark service transformers
            ServiceTransformers = Transformers([Transformers.bus2kV]<1 & ~[Transformers.controlled]);
            transformerBuses = regexprep({ServiceTransformers.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
            Buses = getBusInfo(DSSCircObj,transformerBuses,1);
            BusCoords = reshape([Buses.coordinates],2,[])';
            th = plot(repmat(BusCoords(:,2)',2,1),repmat(BusCoords(:,1)',2,1),'-k>','MarkerSize',7,'MarkerFaceColor',[1,0.5,0.5],'Color',[0.4 0.4 0.4],'LineStyle','none');
            set(th,{'DisplayName'},strcat('transformer.',{ServiceTransformers.name})');
            serviceTransformerGroup = hggroup; set(th,'Parent',serviceTransformerGroup); %group together for one item in the legend
            legendHandles = [legendHandles,serviceTransformerGroup];
            legendText = [legendText,'Service Transformer'];
            FigHandles.serviceTransformer = th;
            FigHandles.serviceTransformerGroup = serviceTransformerGroup;
        else
            FigHandles.serviceTransformer = '';
        end
    end
end


%% Plot transformers if property BoosterMarker is on
if strcmp(BoosterMarker,'on')
    if ischar(Transformers) && strcmp(Transformers,'noInput')
        Transformers = getTransformerInfo(DSSCircObj);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers(ismember(lower(regexprep({Transformers.bus1},'(\.[0-9]+)','')),allBusesToPlot) | ismember(lower(regexprep({Transformers.bus2},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers([Transformers.bus1Distance]>0.01);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        if any([Transformers.bus2kV]>1 & [Transformers.bus2kV]==[Transformers.bus1kV] & ~[Transformers.controlled]) %mark boosters
            Booster = Transformers([Transformers.bus2kV]>1 & [Transformers.bus2kV]==[Transformers.bus1kV] & ~[Transformers.controlled]);
            transformerBuses = regexprep({Booster.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
            Buses = getBusInfo(DSSCircObj,transformerBuses,1);
            BusCoords = reshape([Buses.coordinates],2,[])';
            th = plot(repmat(BusCoords(:,2)',2,1),repmat(BusCoords(:,1)',2,1),'-k<','MarkerSize',10,'MarkerFaceColor',[0.8,0,0],'LineStyle','none');
            set(th,{'DisplayName'},strcat('transformer.',{Booster.name})');
            BoosterGroup = hggroup; set(th,'Parent',BoosterGroup); %group together for one item in the legend
            legendHandles = [legendHandles,BoosterGroup];
            legendText = [legendText,'Booster'];
            FigHandles.booster = th;
            FigHandles.boosterGroup = BoosterGroup;
        else
            FigHandles.booster = '';
        end
    end
end


%% Plot medium voltage transformers if property MVTransformerMarker is on
if strcmp(MVTransformerMarker,'on')
    if ischar(Transformers) && strcmp(Transformers,'noInput')
        Transformers = getTransformerInfo(DSSCircObj);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers(ismember(lower(regexprep({Transformers.bus1},'(\.[0-9]+)','')),allBusesToPlot) | ismember(lower(regexprep({Transformers.bus2},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers([Transformers.bus1Distance]>0.01);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        if any([Transformers.bus2kV]>1 & [Transformers.bus2kV]~=[Transformers.bus1kV] & ~[Transformers.controlled]) %mark medium voltage transformers
            MVTransformer = Transformers([Transformers.bus2kV]>1 & [Transformers.bus2kV]~=[Transformers.bus1kV] & ~[Transformers.controlled]);
            transformerBuses = regexprep({MVTransformer.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
            Buses = getBusInfo(DSSCircObj,transformerBuses,1);
            BusCoords = reshape([Buses.coordinates],2,[])';
            th = plot(repmat(BusCoords(:,2)',2,1),repmat(BusCoords(:,1)',2,1),'-kv','MarkerSize',10,'MarkerFaceColor',[0.8,0,0],'LineStyle','none');
            set(th,{'DisplayName'},strcat('transformer.',{MVTransformer.name})');
            MVTransformerGroup = hggroup; set(th,'Parent',MVTransformerGroup); %group together for one item in the legend
            legendHandles = [legendHandles,MVTransformerGroup];
            legendText = [legendText,'MV Transformer'];
            FigHandles.MVTransformer = th;
            FigHandles.MVTransformerGroup = MVTransformerGroup;
        else
            FigHandles.MVTransformer = '';
        end
    end
end


%% Plot regulators if property RegulatorMarker is on
if strcmp(RegulatorMarker,'on')
    if ischar(Transformers) && strcmp(Transformers,'noInput')
        Transformers = getTransformerInfo(DSSCircObj);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers(ismember(lower(regexprep({Transformers.bus1},'(\.[0-9]+)','')),allBusesToPlot) | ismember(lower(regexprep({Transformers.bus2},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        Transformers = Transformers([Transformers.bus1Distance]>0.01);
    end
    if ~isempty(Transformers) && ~strcmp(Transformers(1).name,'NONE')
        if any([Transformers.controlled]) %mark controlled transformers
            ControlledTransformers = Transformers([Transformers.controlled]);
            transformerBuses = regexprep({ControlledTransformers.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
            Buses = getBusInfo(DSSCircObj,transformerBuses,1);
            BusCoords = reshape([Buses.coordinates],2,[])';
            th = plot(repmat(BusCoords(:,2)',2,1),repmat(BusCoords(:,1)',2,1),'-kd','MarkerSize',9,'MarkerFaceColor','r','LineStyle','none');
            set(th,{'DisplayName'},strcat('transformer.',{ControlledTransformers.name})');
            controlledTransformerGroup = hggroup; set(th,'Parent',controlledTransformerGroup); %group together for one item in the legend
            legendHandles = [legendHandles,controlledTransformerGroup];
            legendText = [legendText,'LTC/VREG'];
            FigHandles.regulators = th;
            FigHandles.regulatorGroup = controlledTransformerGroup;
        else
            FigHandles.regulators = '';
        end
    end
end


%% Plot capacitors if property CapacitorMarker is on
PlotEditToolbarOn = false;
if strcmp(CapacitorMarker,'on')
    if ischar(Capacitors) && strcmp(Capacitors,'noInput')
        Capacitors = getCapacitorInfo(DSSCircObj);
    end
    if ~isempty(Capacitors) && ~strcmp(Capacitors(1).name,'NONE')
        Capacitors = Capacitors(ismember(lower(regexprep({Capacitors.busName},'(\.[0-9]+)','')),allBusesToPlot)); %remove if not on the path of lines being plotted
    end
    if strcmp(SubEquipmentMarker,'off') && ~isempty(Capacitors) && ~strcmp(Capacitors(1).name,'NONE')
        Capacitors = Capacitors([Capacitors.distance]>0.01);
    end
    if ~isempty(Capacitors) && ~strcmp(Capacitors(1).name,'NONE')
        CapacitorCoords = reshape([Capacitors.coordinates],2,[])';
        if any(~[Capacitors.switching]) %plot fixed capacitors
            ch = plot(repmat(CapacitorCoords(~[Capacitors.switching],2)',2,1),repmat(CapacitorCoords(~[Capacitors.switching],1)',2,1),'-ks','MarkerSize',10,'MarkerFaceColor','g','LineStyle','none');
            set(ch,{'DisplayName'},strcat('capacitor.',{Capacitors(~[Capacitors.switching]).name})');
            fixedCapacitorGroup = hggroup; set(ch,'Parent',fixedCapacitorGroup); %group together for one item in the legend
            legendHandles = [legendHandles,fixedCapacitorGroup];
            legendText = [legendText,'Fixed Capacitor'];
            FigHandles.fixedCapacitor = ch;
            FigHandles.fixedCapacitorGroup = fixedCapacitorGroup;
        else
            FigHandles.fixedCapacitor = '';
        end
        if any([Capacitors.switching]) %plot switching capacitors
            ch = plot(repmat(CapacitorCoords([Capacitors.switching],2)',2,1),repmat(CapacitorCoords([Capacitors.switching],1)',2,1),'-k^','MarkerSize',10,'MarkerFaceColor','g','LineStyle','none');
            set(ch,{'DisplayName'},strcat('capacitor.',{Capacitors([Capacitors.switching]).name})');
            switchingCapacitorGroup = hggroup; set(ch,'Parent',switchingCapacitorGroup); %group together for one item in the legend
            legendHandles = [legendHandles,switchingCapacitorGroup];
            legendText = [legendText,'Switching Capacitor'];
            FigHandles.switchingCapacitor = ch;
            FigHandles.switchingCapacitorGroup = switchingCapacitorGroup;
        else
            FigHandles.switchingCapacitor = '';
        end
        if strcmp(CapacitorLabel,'on')
            gcaYlim = ylim;
            gcaHeight = abs(diff(gcaYlim));
            for ii=1:length(Capacitors)
                [arrowx,arrowy] = dsxy2figxy(gca, [CapacitorCoords(ii,2);CapacitorCoords(ii,2)], [CapacitorCoords(ii,1)+gcaHeight/10;CapacitorCoords(ii,1)+gcaHeight/100]);
                ah(ii) = annotation('textarrow',arrowx, arrowy,'HeadStyle','plain','HeadLength',6,'HeadWidth',9,'Color','r','LineWidth',4,...
                    'String',sprintf('%i kVAr',Capacitors(ii).kvar),'TextColor','r','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
                %set(ah,'parent',gca);
                %set(ah,'position',[CapacitorCoords(ii,2),CapacitorCoords(ii,1)+gcaHeight/10 0 -gcaHeight/10]);
            end
            FigHandles.capacitorLabels = ah';
            viewmenufcn('PloteditToolbar');
            PlotEditToolbarOn = true;
            warningState = warning('query','all');
            if ~all(strcmp({warningState.state},'off')) %display popup with more information if warnings aren't turned off
                warndlg('Arrows and text will not automatically move with the figure during panning or zooming.  To link the arrows to the axis, select the ''Pin to axes'' from the toolbar, click an arrow, and repeat for each arrow.');
            end
        end
    end
end


%% Plot end of feeder if property EndOfFeederMarker is on
if strcmp(EndOfFeederMarker,'on')
    Buses = getBusInfo(DSSCircObj,allBusesToPlot);
    Buses = Buses([Buses.numPhases]==3);
    Buses = Buses([Buses.voltage]>50);
    [longestDistance index] = max([Buses.distance]);
    toBus = Buses(index);
    Buses = getBusInfo(DSSCircObj,{toBus.name},1);
    BusCoords = reshape([Buses.coordinates],2,[])';
    eofh = plot(repmat(BusCoords(:,2)',2,1),repmat(BusCoords(:,1)',2,1),'-ko','MarkerSize',9,'MarkerFaceColor','c','LineStyle','none');
    set(eofh,'DisplayName','End Of Feeder');
    legendHandles = [legendHandles,eofh];
    legendText = [legendText,'End of Feeder (3-phase)'];
    FigHandles.endOfFeeder = eofh;
    if strcmp(EndOfFeederLabel,'on')
        gcaYlim = ylim;
        gcaHeight = abs(diff(gcaYlim));
        [arrowx,arrowy] = dsxy2figxy(gca, [BusCoords(1,2);BusCoords(1,2)], [BusCoords(1,1)+gcaHeight/10;BusCoords(1,1)+gcaHeight/100]);
        ah = annotation('textarrow',arrowx, arrowy,'HeadStyle','plain','HeadLength',6,'HeadWidth',9,'Color','r','LineWidth',4,...
            'String',sprintf('%2.2f km',longestDistance),'TextColor','r','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
        FigHandles.endOfFeederLabel = ah;
        if ~PlotEditToolbarOn
            viewmenufcn('PloteditToolbar','on');
            warningState = warning('query','all');
            if ~all(strcmp({warningState.state},'off')) %display popup with more information if warnings aren't turned off
                warndlg('Arrows and text will not automatically move with the figure during panning or zooming.  To link the arrows to the axis, select the ''Pin to axes'' from the toolbar, click an arrow, and repeat for each arrow.');
            end
        end
    end
end


%% Plot custom bus marker CustomMarker CustomLegend 
if ~strcmp(CustomMarker,'off')
    if ischar(CustomMarker)
        CustomMarker = {CustomMarker};
    end
    Buses = getBusInfo(DSSCircObj,CustomMarker,1);
    BusesCoords = reshape([Buses.coordinates],2,[])';
    csh = plot(repmat(BusesCoords(:,2)',2,1),repmat(BusesCoords(:,1)',2,1),'-ko','MarkerSize',10,'MarkerFaceColor','c','LineStyle','none');
    set(csh,'DisplayName',CustomLegend);
    legendHandles = [legendHandles,csh'];
    if isempty(CustomLegend)
        CustomLegend = CustomMarker;
    end
    CustomLegend = reshape(CustomLegend,1,length(CustomLegend));
    legendText = [legendText,CustomLegend];
    FigHandles.customMarker = csh;
end


%% Create legend
if ~isempty(legendHandles)
    lgh = legend(legendHandles,legendText,'Interpreter','none');
    FigHandles.legend = lgh;
    FigHandles.legendHandles = legendHandles;
    FigHandles.legendText = legendText;
end


%% Plot satellite image if property MappingBackground is on
if ~strcmp(MappingBackground,'none')
    if all(bus1Coord(:,1)>-90 & bus1Coord(:,1)<90 & bus1Coord(:,2)>-180 & bus1Coord(:,2)<180)
        axis equal;
        plotGoogleMap('MapType',MappingBackground,'Alpha',0.85);
    else
        warning('Distribution system bus coordinates are not Latitude/Longitude values. Mapping background cannot be displayed');
    end
end


%% Add button to toolbar
if isempty(findall(gcf,'tag','NodeViewToggle'))
    tbh = findall(gcf,'Type','uitoolbar','Tag','FigureToolBar');
    buttonIcon = ones(15,15); buttonIcon(6:9,2:5)=0; buttonIcon(6:9,11:14)=0; buttonIcon(7:8,5:11)=0; buttonIcon = repmat(buttonIcon,[1,1,3]);
    if ~isempty(tbh)
        tth = uitoggletool(tbh,'CData',buttonIcon,'Separator','on','Tag','NodeViewToggle','TooltipString','Turn On Node View');
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
groupChildren =get(findobj(gca,'Type','hggroup')','Children');
if iscell(groupChildren) groupChildren=vertcat(groupChildren{:}); end
set(groupChildren,'buttondownfcn',[removeString,'set(gco,''Selected'',''on''); xCoords = get(gco,''XData''); yCoords = get(gco,''YData'');' addText]);
%if clicking not on an object, remove any text objects and selections
set(gca,'buttondownfcn',removeString)


%% Right Click button
% Set AllowForms to true to allow visualizations (for future versions of OpenDSS)
DSSCircObj.AllowForms = 1;

set(gca,'UserData',{DSSCircObj,Lines});
hcmenu = uicontextmenu;
hcb1 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; if DSSCircObj.ActiveCircuit.SetActiveElement(get(gco,''DisplayName''))>0 DSSCircObj.Text.Command = ''formedit''; else errordlg(sprintf(''Did not find element %s in the circuit.'',get(gco,''DisplayName''))); end'];
hcb2 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.Text.Command = [''visualize voltages element='', get(gco,''DisplayName'')];'];
hcb3 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.Text.Command = [''visualize currents element='', get(gco,''DisplayName'')];'];
hcb4 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; DSSCircObj.Text.Command = [''visualize powers element='', get(gco,''DisplayName'')];'];
hcb5 = ['temp = get(gca,''UserData''); DSSCircObj=temp{1}; Lines=temp{2}; disp(''Plotting circuit lines with marker...''); objectName=get(gco,''DisplayName''); DSSCircObj.ActiveCircuit.SetActiveElement(objectName); figure; plotCircuitLines(DSSCircObj,''CustomMarker'',DSSCircObj.ActiveCircuit.ActiveElement.BusNames{1},''CustomLegend'',objectName,''Lines'',Lines);'];
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






%% Function from MATLAB documentation [docroot '/techdoc/creating_plots/examples']
function varargout = dsxy2figxy(varargin)
% dsxy2figxy -- Transform point or position from data space 
% coordinates into normalized figure coordinates 
% Transforms [x y] or [x y width height] vectors from data space
% coordinates to normalized figure coordinates in order to locate
% annotation objects within a figure. These objects are: arrow, 
% doublearrow, textarrow, ellipse line, rectangle, textbox 
%
% Syntax:
%    [figx figy] = dsxy2figxy([x1 y1],[x2 y2])  % GCA is used
%    figpos      = dsxy2figxy([x1 y1 width height])
%    [figx figy] = dsxy2figxy(axes_handle, [x1 y1],[x2 y2])
%    figpos      = dsxy2figxy(axes_handle, [x1 y1 width height])
%
% Usage: Obtain a position on a plot in data space and  
%        apply this function to locate an annotation there, e.g., 
%   [axx axy] = ginput(2); (input is in data space)
%   [figx figy] = dsxy2figxy(gca, axx, axy);  (now in figure space)
%   har = annotation('textarrow',figx,figy); 
%   set(har,'String',['(' num2str(axx(2)) ',' num2str(axy(2)) ')']) 
%
%   Copyright 2006-2009 The MathWorks, Inc. 

% Obtain arguments (limited argument checking is done)
% Determine if axes handle is specified
if length(varargin{1}) == 1 && ishandle(varargin{1}) ...
                            && strcmp(get(varargin{1},'type'),'axes')	
	hAx = varargin{1};
	varargin = varargin(2:end); % Remove arg 1 (axes handle)
else
	hAx = gca;
end;

% Remaining args are either two point locations or a position vector
if length(varargin) == 1        % Assume a 4-element position vector
	pos = varargin{1};
else
	[x,y] = deal(varargin{:});  % Assume two pairs (start, end points)
end

% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Make axes units normalized 
axpos = get(hAx,'Position');    % Get axes position
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));

% Transform from data space coordinates to normalized figure coordinates 
if exist('x','var')     % Transform a and return pair of points
	varargout{1} = (x - axlim(1)) * axpos(3) / axwidth + axpos(1);
	varargout{2} = (y - axlim(3)) * axpos(4) / axheight + axpos(2);
else                    % Transform and return a position rectangle
	pos(1) = (pos(1) - axlim(1)) / axwidth * axpos(3) + axpos(1);
	pos(2) = (pos(2) - axlim(3)) / axheight * axpos(4) + axpos(2);
	pos(3) = pos(3) * axpos(3) / axwidth;
	pos(4) = pos(4) * axpos(4 )/ axheight;
	varargout{1} = pos;
end

% Restore axes units
set(hAx,'Units',axun)
end