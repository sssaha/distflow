%% plotAmpProfile
% Plots the line currents profile and line rating vs. distance
%
%% Syntax
%  plotAmpProfile(DSSCircObj,BusName);
%  plotAmpProfile(DSSCircObj,BusName, _'PropertyName'_ ,PropertyValue);
%
%% Description
% Function to plot line currents in in each between the selected bus and the substation.  
% The line current and line rating is plotted vs. distance from the substation.  Clicking on objects in
% the figure will display the name of the object, and right clicking will give a menu for viewing
% properties of the object.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|BusName|* - Property for the name of the bus (string) that the current (amp)
% profile should be plotted to.  Only the direct line between the bus and the substation
% will be plotted
% * *|Properties|* - optional properties as one or more name-value pairs in any order
% * -- *|'AveragePhase'|* - Property for if the average power should be
% plotted alone or in addition to the phase plots |'on' | {'off'} | 'addition'|
% * -- *|'Lines'|* - Structure of the circuit lines from getLineInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
%
%% Outputs
% * *none* - a figure is displayed with the plot
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
% Example of an Amp profile plot to a bus
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% figure; plotAmpProfile(DSSCircObj, 'G2101JK1400_N300995_sec_1','AveragePhase','addition')
% ylim([0 1000])
%

function plotAmpProfile(DSSCircObj, BusName, varargin)

%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addRequired('BusName', @ischar);
p.addParamValue('AveragePhase', 'off', @(x)any(strcmp(x,{'on','off','addition'})));
p.addParamValue('Lines', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));

p.parse(DSSCircObj, BusName, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

%% get circuit information
if ischar(Lines) && strcmp(Lines,'noInput')
    Lines = getLineInfo(DSSCircObj);
end
LinesToPlot = Lines;

%% Remove all lines except between BusName and substation
selectedBuses = findUpstreamBuses(DSSCircObj,BusName,'Lines',Lines);
linesBus1 = regexprep({LinesToPlot.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
linesBus2 = regexprep({LinesToPlot.bus2},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
LinesToPlot = LinesToPlot(ismember(linesBus1,selectedBuses) & ismember(linesBus2,selectedBuses));


%% Plot line rating lines
legendHandles = []; %empty variabes to store legend information
legendText = {};

LongLines = LinesToPlot(abs([LinesToPlot.bus2Distance]-[LinesToPlot.bus1Distance])>0.01);

% find and plot all vertical lines at bus connections
distances = [[LongLines.bus1Distance],[LongLines.bus2Distance]];
lineRatings = [[LongLines.lineRating],[LongLines.lineRating]];
[B IA] = sortrows([distances;lineRatings]');
distances = distances(IA);
lineRatings = lineRatings(IA);
[B IA] = sortrows([distances;[0,0,diff(abs(diff(lineRatings)))]]');
plot(distances(IA),lineRatings(IA),'--m','LineWidth',2);
hold all;

h4 = plot([[LongLines.bus1Distance]',[LongLines.bus2Distance]']',[[LongLines.lineRating]',[LongLines.lineRating]']','--m','LineWidth',2);
set(h4,{'DisplayName'},strcat('line.',{LongLines.name})');
ratingsGroup = hggroup; set(h4,'Parent',ratingsGroup); %group together for one item in the legend
legendHandles = [legendHandles,ratingsGroup];
legendText = [legendText,'Line Rating'];


%% Plot average phase lines
if strcmp(AveragePhase,'on') || strcmp(AveragePhase,'addition') %plot average power if requested in properties
    
    % find and plot all vertical lines at bus connections
    currents = [[LinesToPlot.bus1Current],[LinesToPlot.bus1Current]];
    busNames = [{LinesToPlot.bus1},{LinesToPlot.bus2}];
    distances = [[LinesToPlot.bus1Distance],[LinesToPlot.bus2Distance]];
    uniqueBusNames = unique(busNames);
    newLineX = [];
    newLineY = [];
    for ii=1:length(uniqueBusNames)
        [tf loc] = ismember(busNames,uniqueBusNames(ii));
        if sum(tf)>1
            busDistance = mean(distances(tf));
            busCurrents = currents(tf);
            for jj=1:sum(tf)-1
                newLineX = [newLineX; [busDistance,busDistance]];
                newLineY = [newLineY; [busCurrents(1), busCurrents(2)]];
                busCurrents = busCurrents(2:end);
            end
        end
    end
    plot(newLineX',newLineY','-g','LineWidth',2);
    
    h0 = plot([[LinesToPlot.bus1Distance]',[LinesToPlot.bus2Distance]']',[[LinesToPlot.bus1Current]',[LinesToPlot.bus2Current]']','-g','LineWidth',2);
    set(h0,{'DisplayName'},strcat('line.',{LinesToPlot.name})');
    averageGroup = hggroup; set(h0,'Parent',averageGroup); %group together for one item in the legend
    legendHandles = [legendHandles,averageGroup];
    legendText = [legendText,'Average'];
end


%% Plot each phase line
if strcmp(AveragePhase,'off') || strcmp(AveragePhase,'addition') %plot phases
    
    currentPhase = reshape([LinesToPlot.bus1PhaseCurrent],3,[])';
    
    % find the phases on each line
    phaseA = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.1)','Tokens'));
    phaseB = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.2)','Tokens'));
    phaseC = ~cellfun(@isempty,regexp({LinesToPlot.bus1},'.^?(\.3)','Tokens'));
    threePhase = ~phaseA & ~phaseB & ~phaseC; %if no .1.2.3, then three phase
    phaseA(threePhase) = true;
    phaseB(threePhase) = true;
    phaseC(threePhase) = true;
    
    % find and plot all vertical lines at bus connections
    currentPhases = [currentPhase;currentPhase];
    busNames = [{LinesToPlot.bus1},{LinesToPlot.bus2}];
    distances = [[LinesToPlot.bus1Distance],[LinesToPlot.bus2Distance]];
    uniqueBusNames = unique(busNames);
    newLineX = [];
    newLineY = [];
    for ii=1:length(uniqueBusNames)
        [tf loc] = ismember(busNames,uniqueBusNames(ii));
        if sum(tf)>1
            busDistance = mean(distances(tf));
            busPowers = currentPhases(tf,:);
            for jj=1:sum(tf)-1
                newLineX = [newLineX; [busDistance,busDistance]];
                newLineY = [newLineY; [busPowers(1,:), busPowers(2,:)]];
                busPowers = busPowers(2:end,:);
            end
        end
    end
    newLineX = [newLineX,newLineX,newLineX];
    plot(newLineX(:,[1,4])',newLineY(:,[1,4])','-k','LineWidth',2);
    plot(newLineX(:,[2,5])',newLineY(:,[2,5])','-r','LineWidth',2);
    plot(newLineX(:,[3,6])',newLineY(:,[3,6])','-b','LineWidth',2);
    
    % plot
    h1 = plot([[LinesToPlot(phaseA).bus1Distance]',[LinesToPlot(phaseA).bus2Distance]']',[currentPhase(phaseA,1),currentPhase(phaseA,1)]','-k','LineWidth',2);
    set(h1,{'DisplayName'},strcat('line.',{LinesToPlot(phaseA).name})');
    h2 = plot([[LinesToPlot(phaseB).bus1Distance]',[LinesToPlot(phaseB).bus2Distance]']',[currentPhase(phaseB,2),currentPhase(phaseB,2)]','-r','LineWidth',2);
    set(h2,{'DisplayName'},strcat('line.',{LinesToPlot(phaseB).name})');
    h3 = plot([[LinesToPlot(phaseC).bus1Distance]',[LinesToPlot(phaseC).bus2Distance]']',[currentPhase(phaseC,3),currentPhase(phaseC,3)]','-b','LineWidth',2);
    set(h3,{'DisplayName'},strcat('line.',{LinesToPlot(phaseC).name})');
    phaseAGroup = hggroup; set(h1,'Parent',phaseAGroup); %group together for one item in the legend
    phaseBGroup = hggroup; set(h2,'Parent',phaseBGroup); %group together for one item in the legend
    phaseCGroup = hggroup; set(h3,'Parent',phaseCGroup); %group together for one item in the legend
    legendHandles = [legendHandles,phaseAGroup];
    legendText = [legendText,'PhaseA'];
    legendHandles = [legendHandles,phaseBGroup];
    legendText = [legendText,'PhaseB'];
    legendHandles = [legendHandles,phaseCGroup];
    legendText = [legendText,'PhaseC'];
end


%% Plot Edits
grid on;
set(gca,'FontSize',10,'FontWeight','bold')
xlabel(gca,'Distance from Substation (km)','FontSize',12,'FontWeight','bold')
ylabel(gca,'Line Current (Amps)','FontSize',12,'FontWeight','bold')
title('Feeder Amp Profile','FontWeight','bold','FontSize',12);
legend(legendHandles,legendText);


%% Add button to toolbar
if isempty(findall(gcf,'tag','NodeViewToggle'))
    tbh = findall(gcf,'Type','uitoolbar');
    buttonIcon = ones(15,15); buttonIcon(6:9,2:5)=0; buttonIcon(6:9,11:14)=0; buttonIcon(7:8,5:11)=0; buttonIcon = repmat(buttonIcon,[1,1,3]);
    tth = uitoggletool(tbh,'CData',buttonIcon,'Separator','on','Tag','NodeViewToggle','TooltipString','Turn On Node View');
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

end