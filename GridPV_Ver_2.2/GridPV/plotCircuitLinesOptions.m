%% plotCircuitLinesOptions
% GUI for providing options for how to plot the feeder circuit diagram
%
%% Syntax
%  plotCircuitLinesOptions(DSSCircObj);
%
%% Description
% plotCircuitLines plots the feeder circuit diagram and has many different input argument parameters for changing coloring, line
% thickness, background, etc.  This function provides a GUI for selecting
% the plotting styles for plotCircuitLines instead of through text arguments.  This function can be called directly with the OpenDSS circuit object,
% or plotCircuitLines.m will call this function if no input arguments were selected.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
%
%% Outputs
% * *none* - a figure of the circuit is displayed based on the option inputs
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
% Examples of calling the GUI for plotCircuitLinesOptions
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% figure; plotCircuitLinesOptions(DSSCircObj)
%

function varargout = plotCircuitLinesOptions(varargin)
global calledFromPlotCircuitLines DSSCircObj
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotCircuitLinesOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @plotCircuitLinesOptions_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargin==2 && iscell(varargin{2})
    if strcmp('GiveOptions',varargin{2})
        calledFromPlotCircuitLines=true;
    end
end

if nargin>=1 && strcmp(class(varargin{1}),'COM.OpenDSSEngine_DSS')
    DSSCircObj=varargin{1};
    if nargin==1
        calledFromPlotCircuitLines = false;
    end
end


if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

%% --- Executes just before plotCircuitLinesOptions is made visible.
function plotCircuitLinesOptions_OpeningFcn(hObject, eventdata, handles, varargin)
global DSSCircObj

% Choose default command line output for plotCircuitLinesOptions
handles.output = hObject;

set(handles.uipanelColoring,'SelectionChangeFcn',@coloring_Callback);

energyMeterNames = DSSCircObj.ActiveCircuit.Meters.AllNames;
set(handles.listboxEnergyMeters,'String',energyMeterNames);
set(handles.listboxEnergyMeters,'Value',1:length(energyMeterNames));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plotCircuitLinesOptions wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

%% --- Outputs from this function are returned to the command line.
function varargout = plotCircuitLinesOptions_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%% --- Executes on button press in pushbuttonPlot.
function pushbuttonPlot_Callback(hObject, eventdata, handles)
global Coloring ContourScale Thickness SubstationMarker SubEquipmentMarker PVMarker GeneratorMarker LoadMarker RegulatorMarker MVTransformerMarker BoosterMarker ServiceTransformerMarker CapacitorMarker CapacitorLabel EndOfFeederMarker EndOfFeederLabel CustomMarker CustomLegend EnergyMeter NumPhases PhasesToPlot MappingBackground BusName Downstream calledFromPlotCircuitLines DSSCircObj

%% Coloring
if get(handles.radiobuttonColorSpec,'Value')
    Coloring = str2num(get(handles.editRGB,'String'));
elseif get(handles.radiobuttonNumPhases,'Value')
    Coloring = 'numPhases';
elseif get(handles.radiobuttonPerPhase,'Value')
    Coloring = 'perPhase';
elseif get(handles.radiobuttonEnergyMeter,'Value')
    Coloring = 'energyMeter';
elseif get(handles.radiobuttonContour,'Value')
    switch get(handles.popupmenuColorContour,'Value')
        case 1
            Coloring = 'voltage120';
        case 2
            Coloring = 'voltagePU';
        case 3
            Coloring = 'voltage';
        case 4
            Coloring = 'voltage120LL';
        case 5
            Coloring = 'voltagePULL';
        case 6
            Coloring = 'voltageLL';
        case 7
            Coloring = 'lineLoading';
        case 8
            Coloring = 'realLosses';
        case 9
            Coloring = 'reactiveLosses';
        case 10
            Coloring = 'distance';
        case 11
            Coloring = 'unbalance';
        case 12
            Coloring = 'voltageAngle';
        case 13
            Coloring = 'powerFactor';
        case 14
            Coloring = 'powerFlowDirection';
        case 15
            Coloring = 'impedance';
        case 16
            Coloring = 'resistance';
        case 17
            Coloring = 'reactance';
        case 18
            Coloring = 'faultCurrent3P';
        case 19
            Coloring = 'faultCurrent1P';
        case 20
            Coloring = 'faultCurrentLL';
    end
end

%% Contour options
if get(handles.checkboxAuto,'Value')
    ContourScale = 'auto';
else
    ContourScale = [str2double(get(handles.editMin,'String')),str2double(get(handles.editMax,'String'))];
end

%% Thickness
if get(handles.radiobuttonFixed,'Value')
    Thickness = str2num(get(handles.editThickness,'String'));
elseif get(handles.radiobuttonNumPhasesThick,'Value')
    Thickness = 'numPhases';
elseif get(handles.radiobuttonCurrent,'Value')
    Thickness = 'current';
elseif get(handles.radiobuttonLineRating,'Value')
    Thickness = 'lineRating';
end

%% Markers
if get(handles.checkboxSubstation,'Value') SubstationMarker = 'on'; else SubstationMarker = 'off'; end
if get(handles.checkboxSubEquipment,'Value') SubEquipmentMarker = 'off'; else SubEquipmentMarker = 'on'; end

if get(handles.checkboxCapacitor,'Value')
    CapacitorMarker = 'on';
    if get(handles.checkboxCapacitorLabel,'Value') CapacitorLabel = 'on'; else CapacitorLabel = 'off'; end
else
    CapacitorMarker = 'off';
    CapacitorLabel = 'off';
end
if ~get(handles.checkboxTransformers,'Value')
    RegulatorMarker = 'off';
    MVTransformerMarker = 'off';
    BoosterMarker = 'off';
    ServiceTransformerMarker = 'off';
else
    if get(handles.checkboxRegulators,'Value') RegulatorMarker = 'on'; else RegulatorMarker = 'off'; end
    if get(handles.checkboxMVTransformers,'Value') MVTransformerMarker = 'on'; else MVTransformerMarker = 'off'; end
    if get(handles.checkboxBoosters,'Value') BoosterMarker = 'on'; else BoosterMarker = 'off'; end
    if get(handles.checkboxServiceTransformers,'Value') ServiceTransformerMarker = 'on'; else ServiceTransformerMarker = 'off'; end
end
if get(handles.checkboxLoads,'Value') LoadMarker = 'on'; else LoadMarker = 'off'; end
if get(handles.checkboxGenerators,'Value') GeneratorMarker = 'on'; else GeneratorMarker = 'off'; end
if get(handles.checkboxPV,'Value') PVMarker = 'on'; else PVMarker = 'off'; end
if get(handles.checkboxEndOfFeeder,'Value')
    EndOfFeederMarker = 'on';
    if get(handles.checkboxEndOfFeederLabel,'Value') EndOfFeederLabel = 'on'; else EndOfFeederLabel = 'off'; end
else
    EndOfFeederMarker = 'off';
    EndOfFeederLabel = 'off';
end

if get(handles.checkboxCustom,'Value') 
    CustomMarker = get(handles.editCustomBus,'String');
    CustomLegend = get(handles.editCustomLabel,'String');
else
    CustomMarker = 'off';
    CustomLegend = '';
end

%% Background
if get(handles.popupmenuBackground,'Value')==1
    MappingBackground = 'none';
elseif get(handles.popupmenuBackground,'Value')==2
    MappingBackground = 'hybrid';
elseif get(handles.popupmenuBackground,'Value')==3
    MappingBackground = 'satellite';
elseif get(handles.popupmenuBackground,'Value')==4
    MappingBackground = 'roadmap';
elseif get(handles.popupmenuBackground,'Value')==5
    MappingBackground = 'terrain';
end

%% EnergyMeters to plot
EnergyMeterStrings = get(handles.listboxEnergyMeters,'String');
energyMeterIndexes = get(handles.listboxEnergyMeters,'Value');
if length(EnergyMeterStrings) == length(energyMeterIndexes)
    EnergyMeter = 'all';
else
    EnergyMeter = EnergyMeterStrings(energyMeterIndexes);
end

%% Phases to plot
NumPhases = [];
if get(handles.checkboxNumPhase1,'Value') NumPhases=[NumPhases,1]; end
if get(handles.checkboxNumPhase2,'Value') NumPhases=[NumPhases,2]; end
if get(handles.checkboxNumPhase3,'Value') NumPhases=[NumPhases,3]; end
PhasesToPlot = [get(handles.checkboxPhaseA,'Value'),get(handles.checkboxPhaseB,'Value'),get(handles.checkboxPhaseC,'Value')];

%% Upstream/downstream
if get(handles.checkboxBusName,'Value') BusName=get(handles.editBusName,'String'); end
if isempty(BusName) BusName='all'; end
if get(handles.checkboxDownstream,'Value') Downstream='on'; else Downstream='off'; end

if calledFromPlotCircuitLines %if called from plotCircuitLines, no need to call it again
    close;
else
    close;
    drawnow;
    figure;
    plotCircuitLines(DSSCircObj,'Coloring',Coloring,'ContourScale',ContourScale,'Thickness',Thickness,'SubstationMarker',SubstationMarker,'SubEquipmentMarker',SubEquipmentMarker,'PVMarker',PVMarker,'GeneratorMarker',GeneratorMarker,'LoadMarker',LoadMarker,...
        'RegulatorMarker',RegulatorMarker,'MVTransformerMarker',MVTransformerMarker,'BoosterMarker', BoosterMarker,'ServiceTransformerMarker',ServiceTransformerMarker, 'CapacitorMarker',CapacitorMarker,'CapacitorLabel',CapacitorLabel,...
        'EndOfFeederMarker',EndOfFeederMarker,'EndOfFeederLabel',EndOfFeederLabel,'CustomMarker',CustomMarker,'CustomLegend',CustomLegend,'EnergyMeter',EnergyMeter,'NumPhases',NumPhases,'PhasesToPlot',PhasesToPlot,'MappingBackground',MappingBackground,'BusName',BusName,'Downstream',Downstream);
end

end


%%
function editMax_Callback(hObject, eventdata, handles)
set(handles.checkboxAuto,'Value',0)
end


function editMin_Callback(hObject, eventdata, handles)
set(handles.checkboxAuto,'Value',0)
end


%% --- Executes on button press in checkboxBusName.
function checkboxBusName_Callback(hObject, eventdata, handles)
if get(handles.checkboxBusName,'Value')
    set(handles.checkboxDownstream,'Enable','on');
    set(handles.editBusName,'Enable','on');
else
    set(handles.checkboxDownstream,'Enable','off');
    set(handles.editBusName,'Enable','off');
end
end


%% --- Executes on button press in checkboxBusName.
function coloring_Callback(hObject, eventdata, handles)
selectedObject = get(get(hObject,'SelectedObject'),'tag');
if strcmp(selectedObject,'radiobuttonColorSpec') || strcmp(selectedObject,'radiobuttonNumPhases') || strcmp(selectedObject,'radiobuttonPerPhase') || strcmp(selectedObject,'radiobuttonEnergyMeter')
    set(get(findobj(gcf,'Tag','uipanelContour'),'Children'),'Enable','off');
    set(findobj(gcf,'Tag','popupmenuColorContour'),'Enable','off');

else
    set(get(findobj(gcf,'Tag','uipanelContour'),'Children'),'Enable','on');
    set(findobj(gcf,'Tag','popupmenuColorContour'),'Enable','on');
end
if strcmp(selectedObject,'radiobuttonColorSpec')
    set(findobj(gcf,'Tag','text1'),'Enable','on'); set(findobj(gcf,'Tag','editRGB'),'Enable','on');
else
    set(findobj(gcf,'Tag','text1'),'Enable','off'); set(findobj(gcf,'Tag','editRGB'),'Enable','off');
end
end


%% --- Executes on button press in checkboxCustom.
function checkboxCustom_Callback(hObject, eventdata, handles)
if ~get(handles.checkboxCustom,'Value')
    set(handles.text10,'Enable','off'); set(handles.editCustomBus,'Enable','off');
    set(handles.text11,'Enable','off'); set(handles.editCustomLabel,'Enable','off');
else
    set(handles.text10,'Enable','on'); set(handles.editCustomBus,'Enable','on');
    set(handles.text11,'Enable','on'); set(handles.editCustomLabel,'Enable','on');
end
end


%% --- Executes on button press in checkboxTransformers.
function checkboxTransformers_Callback(hObject, eventdata, handles)
if ~get(handles.checkboxTransformers,'Value')
    set(handles.checkboxRegulators,'Enable','off');
    set(handles.checkboxMVTransformers,'Enable','off');
    set(handles.checkboxBoosters,'Enable','off');
    set(handles.checkboxServiceTransformers,'Enable','off');
else
    set(handles.checkboxRegulators,'Enable','on');
    set(handles.checkboxMVTransformers,'Enable','on');
    set(handles.checkboxBoosters,'Enable','on');
    set(handles.checkboxServiceTransformers,'Enable','on');
end
end


%% --- Executes on button press in checkboxCapacitor.
function checkboxCapacitor_Callback(hObject, eventdata, handles)
if ~get(handles.checkboxCapacitor,'Value')
    set(handles.checkboxCapacitorLabel,'Enable','off');
else
    set(handles.checkboxCapacitorLabel,'Enable','on');
end
end


%% --- Executes on button press in checkboxEndOfFeeders.
function checkboxEndOfFeeder_Callback(hObject, eventdata, handles)
if ~get(handles.checkboxEndOfFeeder,'Value')
    set(handles.checkboxEndOfFeederLabel,'Enable','off');
else
    set(handles.checkboxEndOfFeederLabel,'Enable','on');
end
end

