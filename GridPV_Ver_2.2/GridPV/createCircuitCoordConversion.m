%% createCircuitCoordConversion
% Function to create conversion of circuit coordinates to GPS coordinates
%
%% Syntax
%  createCircuitCoordConversion();
%
%% Description
% Function is a user interface to map the Google map and the circuit 
% drawing on top of each other.  The user aligns the two images and the 
% conversion is created for getting GPS Lat/Lon for the OpenDSS bus 
% coordinates. This is used when the OpenDSS coordinate system is unknown 
% and not any standard coordinate systems like UTM.
%
%% Inputs
% * *none* - user will select the OpenDSS circuit file along with the coordinates file through the GUI
%
%% Outputs
% * *none*  - a new OpenDSS bus coordinates files is saved out for the circuit
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
% Starts the user interface.  Directions are in the interface.
%%
% createCircuitCoordConversion();
%

function varargout = createCircuitCoordConversion(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @createCircuitCoordConversion_OpeningFcn, ...
                   'gui_OutputFcn',  @createCircuitCoordConversion_OutputFcn, ...
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
end

%% --- Executes just before createCircuitCoordConversion is made visible.
function createCircuitCoordConversion_OpeningFcn(hObject, eventdata, handles, varargin)
global DSSCircObj DSSCircuit DSSText

% Choose default command line output for createCircuitCoordConversion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

axes(handles.axes1);
plotGoogleMap('MapType','hybrid');

% UIWAIT makes createCircuitCoordConversion wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end


%% --- Outputs from this function are returned to the command line.
function varargout = createCircuitCoordConversion_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
end


%% --- Executes on button press in pushPlotCircuit.
function pushPlotCircuit_Callback(hObject, eventdata, handles)
global DSSCircObj DSSCircuit DSSText xlimits ylimits FileName PathName

% get circuit file
[FileName,PathName] = uigetfile('*.dss','Basecase: Select the OpenDSS file with the circuit');
location = cd;
DSSText.command = sprintf('Compile (%s)',[PathName,FileName]);
DSSText.command = 'solve';
cd(location);

% plot over current map
axes(handles.axes2);
position = get(gca,'Position');
plotCircuitLines(DSSCircObj,'Coloring','perPhase');
set(gca,'Color','none');
set(gca,'Position',position);
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
colorbar('off');

end






%% --- Executes on slider movement.
function sliderXZoom_Callback(hObject, eventdata, handles)
global xlimits ylimits

set(gca,'xlim',[xlimits(1) 1/get(handles.sliderXZoom,'Value')*(xlimits(2)-xlimits(1))+xlimits(1)]-get(handles.sliderXposition,'Value')*(xlimits(2)-xlimits(1)));
set(gca,'ylim',[ylimits(1) 1/get(handles.sliderYzoom,'Value')*(ylimits(2)-ylimits(1))+ylimits(1)]-get(handles.sliderYposition,'Value')*(ylimits(2)-ylimits(1)));

end

%% --- Executes on slider movement.
function sliderXposition_Callback(hObject, eventdata, handles)
global xlimits ylimits

set(gca,'xlim',[xlimits(1) 1/get(handles.sliderXZoom,'Value')*(xlimits(2)-xlimits(1))+xlimits(1)]-get(handles.sliderXposition,'Value')*(xlimits(2)-xlimits(1)));

end


%% --- Executes on slider movement.
function sliderYposition_Callback(hObject, eventdata, handles)
global xlimits ylimits

set(gca,'ylim',[ylimits(1) 1/get(handles.sliderYzoom,'Value')*(ylimits(2)-ylimits(1))+ylimits(1)]-get(handles.sliderYposition,'Value')*(ylimits(2)-ylimits(1)));

end


%% --- Executes on slider movement.
function sliderYzoom_Callback(hObject, eventdata, handles)
global xlimits ylimits

set(gca,'xlim',[xlimits(1) 1/get(handles.sliderXZoom,'Value')*(xlimits(2)-xlimits(1))+xlimits(1)]-get(handles.sliderXposition,'Value')*(xlimits(2)-xlimits(1)));
set(gca,'ylim',[ylimits(1) 1/get(handles.sliderYzoom,'Value')*(ylimits(2)-ylimits(1))+ylimits(1)]-get(handles.sliderYposition,'Value')*(ylimits(2)-ylimits(1)));

end



%% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global DSSCircObj DSSCircuit DSSText FileName PathName

axes(handles.axes1);
xlimitsGPS = get(gca,'xlim');
ylimitsGPS = get(gca,'ylim');

axes(handles.axes2);
xlimitsCircuit = get(gca,'xlim');
ylimitsCircuit = get(gca,'ylim');

CircuitConversion.xScale = (xlimitsGPS(2)-xlimitsGPS(1))/(xlimitsCircuit(2)-xlimitsCircuit(1));
CircuitConversion.xIntercept = xlimitsGPS(1)-CircuitConversion.xScale*xlimitsCircuit(1);

CircuitConversion.yScale = (ylimitsGPS(2)-ylimitsGPS(1))/(ylimitsCircuit(2)-ylimitsCircuit(1));
CircuitConversion.yIntercept = ylimitsGPS(1)-CircuitConversion.yScale*ylimitsCircuit(1);

%Get the current coordinates
[busCoordNames coordinates] = getBusCoordinatesArray(DSSCircObj);

latitude = CircuitConversion.yScale*coordinates(:,2)+CircuitConversion.yIntercept;
longitude = CircuitConversion.xScale*coordinates(:,1)+CircuitConversion.xIntercept;


[fnm pnm] = uigetfile({'*.dss; *.txt; *.xls; *.csv'}, 'Select your bus coordinates file');
if ~(isequal(fnm, 0) || isequal(pnm, 0))
    [SUCCESS,MESSAGE,MESSAGEID] = movefile([pnm fnm], [pnm 'BACKUP_' fnm], 'f');
    if ~SUCCESS
        disp(MESSAGE);
        waitfor(warndlg(sprintf('A backup of your coordinates file was not made (see command window for the message).\n\nPlease manually make a backup copy of your coordinates file before clicking ''OK'' to this message.'), 'No Backup Made!'));
    end

    %Write to file
    %Could not find a way to avoid iteration
    fid = fopen([pnm fnm], 'w');
    for ii = 1:length(busCoordNames)
        fprintf(fid, '%s %3.15f %3.15f\r\n', busCoordNames{ii}, longitude(ii), latitude(ii));
    end
    fclose(fid);

    location = cd;
    DSSText.command = sprintf('Compile (%s)',[PathName,FileName]);
    DSSText.command = 'solve';
    cd(location);
    
    msgbox('Conversion complete. Save was successful.', 'Done');
    uiwait;
    close(handles.figure1);
    
end

end




%% --------------------------------------------------------------------
function uitoggletool1_OffCallback(hObject, eventdata, handles)
axes(handles.axes2)
end

%% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles)
axes(handles.axes1);
end

%% --------------------------------------------------------------------
function uitoggletool2_OffCallback(hObject, eventdata, handles)
axes(handles.axes2)
end

%% --------------------------------------------------------------------
function uitoggletool2_OnCallback(hObject, eventdata, handles)
axes(handles.axes1);
end

%% --------------------------------------------------------------------
function uitoggletool3_OffCallback(hObject, eventdata, handles)
axes(handles.axes2)
end

%% --------------------------------------------------------------------
function uitoggletool3_OnCallback(hObject, eventdata, handles)
axes(handles.axes1);
end
