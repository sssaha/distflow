%% makeVVCcurve
% GUI for setting up the OpenDSS VVControl function parameters
%
%% Syntax
%  makeVVCcurve()
%
%% Description
% This function is a user interface to create the Volt/Var control function
% in OpenDSS.  The required parameters are entered into the interface and a
% mat file is saved with the parameters.  This function is often called
% from placePVplant.m when the PV plant power factor control is selected.
% The saved mat file is used in createPVscenarioFiles.m when the solar scenario OpenDSS
% generators are created.
%
%% Inputs
% * *none*
%
%% Outputs
% * *none* saves a *.mat file with the VVControl parameters
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
% makeVVCcurve()
%

function varargout = makeVVCcurve(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @makeVVCcurve_OpeningFcn, ...
                   'gui_OutputFcn',  @makeVVCcurve_OutputFcn, ...
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

%% --- Executes just before makeVVCcurve is made visible.
function makeVVCcurve_OpeningFcn(hObject, eventdata, handles, varargin)
global voltage VarOutput

voltage = [0.9,0.95,0.98,1.02,1.05,1.10];
VarOutput = [1,1,0,0,-1,-1];

set(handles.uitable1,'data',[voltage',VarOutput']);

plotLine();
% Choose default command line output for makeVVCcurve
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

end


%% --- Outputs from this function are returned to the command line.
function varargout = makeVVCcurve_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
end



%%
function plotLine
global voltage VarOutput


plot(voltage,VarOutput,'LineWidth',2,'Marker','square','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','b');
grid on;
xlim([0.85 1.15]);
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('Voltage');
ylabel('VARs Generated');
end


%% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global voltage VarOutput

inverterKVApu = str2num(get(handles.editInverterKVA,'string'));

[file,path] = uiputfile('*.mat','Save VVControl As','VVControlInfo.mat');
save([path,file],'voltage','VarOutput','inverterKVApu');

end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
global voltage VarOutput

table = get(handles.uitable1,'data')';

voltage = table(1,:);
VarOutput = table(2,:);

plotLine();
end
