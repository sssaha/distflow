%% placePVplant
% Draw PV on the circuit diagram and save plant info for WVM input
%
%% Syntax
%  placePVplant()
%
%% Description
% This function is a user interface where the PV plant can be drawn on the
% circuit diagram.  The user will setup all the PV plant info
% and save it to a file for running WVM.
%
%% Inputs
% * *none*
%
%% Outputs
% * *none* saves a *.mat file with the structure plantinfo for input to WVM
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
% Showing the user interface:
%%
% placePVplant()
%

function varargout = placePVplant(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @placePVplant_OpeningFcn, ...
                   'gui_OutputFcn',  @placePVplant_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

%% --- Executes just before placePVplant is made visible.
function placePVplant_OpeningFcn(hObject, eventdata, handles, varargin)
global DSSCircObj DSSCircuit DSSText figHandles

% Choose default command line output for placePVplant
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

% get circuit file
[FileName,PathName] = uigetfile('*.dss','Basecase: Select the OpenDSS file with the circuit');
location = cd;
DSSText.command = sprintf('Compile (%s)',[PathName,FileName]);
cd(location);
DSSText.command = 'solve';

axes(handles.axes1);
plotCircuitLines(DSSCircObj,'Coloring','perPhase','Thickness',3);

latRange = ylim;
lonRange = xlim;
if latRange(1)<-90 || latRange(2)>90 || lonRange(1)<-180 || lonRange(2)>180
    error('Distribution system bus coordinates are not Latitude/Longitude values. Use initCoordConversion to convert the feeder to Lat/Lon');
end
plotGoogleMap('MapType','hybrid','Alpha',0.85);

figHandles = handles;
% UIWAIT makes placePVplant wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end


%% --- Outputs from this function are returned to the command line.
function varargout = placePVplant_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%%
function editCentralMW_Callback(hObject, eventdata, handles)
drawrectangle();
end


%%
function editCentralDensity_Callback(hObject, eventdata, handles)
drawrectangle();
end



%%
function editDistributedMW_Callback(hObject, eventdata, handles)
global plantinfo

plantinfo.MW = str2num(get(handles.editDistributedMW,'String'));
calcDistributedSize(0);
end



%% --- Executes on button press in radiobuttonCentral.
function radiobuttonCentral_Callback(hObject, eventdata, handles)
global plantinfo figHandles

if isfield(figHandles,'rectangle')
    delete(figHandles.rectangle);
    figHandles = rmfield(figHandles,'rectangle');
end
if isfield(figHandles,'polyHandle')
    delete(figHandles.polyHandle);
    figHandles = rmfield(figHandles,'polyHandle');
end

h = helpdlg('Click on the bus that the central plant will be connected.');
uiwait(h);

plantinfo.type='square';

set(handles.radiobuttonCentral,'Value',1);

% get user click
[Lon,Lat] = ginput(1);
plantinfo.Lat = Lat;
plantinfo.Lon = Lon;

drawrectangle();
end


%% --- Executes on button press in radiobuttonDistributed.
function radiobuttonDistributed_Callback(hObject, eventdata, handles)
global plantinfo figHandles

if isfield(figHandles,'rectangle')
    delete(figHandles.rectangle);
    figHandles = rmfield(figHandles,'rectangle');
end

plantinfo.type='polygon';

set(handles.radiobuttonDistributed,'Value',1);

if exist('impoly') %if has image processing toolbox
    h = helpdlg('Click on each vertex of the polygon surrounding the area.');
    uiwait(h);
    if isfield(figHandles,'polyHandle')
        delete(figHandles.polyHandle);
        figHandles = rmfield(figHandles,'polyHandle');
    end
    figHandles.polyHandle = impoly;
    addNewPositionCallback(figHandles.polyHandle,@calcDistributedSize);
    polygons = getPosition(figHandles.polyHandle);
    plantinfo.Lon = polygons(:,1);
    plantinfo.Lat = polygons(:,2);

else
    h = helpdlg('Click on each vertex of the polygon surrounding the area.  Hit enter when done.');
    uiwait(h);
    
    plantinfo.Lon = [];
    plantinfo.Lat = [];
    [x,y] = ginput(1);
    while ~isempty(x)
        plantinfo.Lon = [plantinfo.Lon; x];
        plantinfo.Lat = [plantinfo.Lat; y];
        if isfield(figHandles,'polyHandle')
            delete(figHandles.polyHandle);
            figHandles = rmfield(figHandles,'polyHandle');
        end
        figHandles.polyHandle = plot(plantinfo.Lon,plantinfo.Lat,'LineWidth',2);
        [x,y] = ginput(1);
    end
end


plantinfo.MW = str2num(get(handles.editDistributedMW,'String'));
plantinfo.PVdensity = str2num(get(handles.editDistributedDensity,'String'));
calcDistributedSize(0);
end


%% Calculates the size of the distributed system
function calcDistributedSize(position)
global plantinfo figHandles

if isobject(figHandles.polyHandle)
    polygons = getPosition(figHandles.polyHandle);
    plantinfo.Lon = polygons(:,1);
    plantinfo.Lat = polygons(:,2);
end
[x,y,utmzone] = deg2utm(plantinfo.Lat,plantinfo.Lon);
area = polyarea(x,y);

plantinfo.PVdensity = plantinfo.MW*1e6/(area*100);
set(figHandles.editDistributedDensity,'String',num2str(plantinfo.PVdensity));

end


%% place the central pv plant square
function drawrectangle
global plantinfo figHandles

[x,y,utmzone] = deg2utm(plantinfo.Lat,plantinfo.Lon);

plantinfo.MW = str2num(get(figHandles.editCentralMW,'String'));
plantinfo.PVdensity = str2num(get(figHandles.editCentralDensity,'String'));

area = plantinfo.MW*1e6/(plantinfo.PVdensity*100); %m^2

squareLength = sqrt(area);

[LatStart,LonStart] = utm2deg(x-squareLength/2,y-squareLength/2,utmzone);
[LatEnd,LonEnd] = utm2deg(x+squareLength/2,y+squareLength/2,utmzone);

if isfield(figHandles,'rectangle')
    delete(figHandles.rectangle);
    figHandles = rmfield(figHandles,'rectangle');
end
if isfield(figHandles,'polyHandle')
    delete(figHandles.polyHandle);
    figHandles = rmfield(figHandles,'polyHandle');
end
figHandles.rectangle = rectangle('Position',[LonStart,LatStart,LonEnd-LonStart,LatEnd-LatStart],'LineWidth',2,'EdgeColor','b');

end


%% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global plantinfo figHandles

plantinfo.tilt = str2num(get(figHandles.editTilt,'String'));
plantinfo.azimuth = str2num(get(figHandles.editAzimuth,'String'));

if get(figHandles.radiobuttonPFfixed,'Value')
    plantinfo.powerFactor.type = 'fixed';
    plantinfo.powerFactor.value = str2num(get(figHandles.editPowerFactor,'String'));
elseif get(figHandles.radiobuttonPFschedule,'Value')
    plantinfo.powerFactor.type = 'schedule';
    plantinfo.powerFactor.filepath = get(handles.editPFscheduleFile,'string');
elseif get(figHandles.radiobuttonPFfunction,'Value')
    plantinfo.powerFactor.type = 'function';
    plantinfo.powerFactor.filepath = get(handles.editPFfunctionFile,'string');
elseif get(figHandles.radiobuttonVVControl,'Value')
    plantinfo.powerFactor.type = 'VVControl';
    plantinfo.powerFactor.filepath = get(handles.editVVControlFile,'string');
end

[file,path] = uiputfile('*.mat','Save Plant Info As');
save([path,file],'plantinfo');

msgbox('PV plant info file created. Save was successful.', 'Done');

end


%% --- Executes on button press in radiobuttonPFfixed.
function radiobuttonPFfixed_Callback(hObject, eventdata, handles)

end


%% --- Executes on button press in radiobuttonPFschedule.
function radiobuttonPFschedule_Callback(hObject, eventdata, handles)

end


%% --- Executes on button press in radiobuttonPFfunction.
function radiobuttonPFfunction_Callback(hObject, eventdata, handles)

end


%% --- Executes on button press in pushbuttonMakeSchedule.
function pushbuttonMakeSchedule_Callback(hObject, eventdata, handles)
makePFschedule();
end


%% --- Executes on button press in pushbuttonMakeFunction.
function pushbuttonMakeFunction_Callback(hObject, eventdata, handles)
makePFoutputFunction();
end


%%
function editPFscheduleFile_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile('*.mat','Select the Power Factor Schedule file');
if FileName~=0
    set(handles.editPFscheduleFile,'string',[PathName,FileName]);
end
end

%%
function editPFfunctionFile_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile('*.mat','Select the Power Factor Function file');
if FileName~=0
    set(handles.editPFfunctionFile,'string',[PathName,FileName]);
end
end


%% --- Executes on button press in pushbuttonMakeVVControl.
function pushbuttonMakeVVControl_Callback(hObject, eventdata, handles)
makeVVCcurve();
end


%%
function editVVControlFile_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile('*.mat','Select the VV Control Info file');
if FileName~=0
    set(handles.editVVControlFile,'string',[PathName,FileName]);
end
end


% --- Executes on button press in checkboxTracking.
function checkboxTracking_Callback(hObject, eventdata, handles)

end