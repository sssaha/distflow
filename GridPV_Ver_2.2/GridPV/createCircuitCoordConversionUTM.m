%% createCircuitCoordConversionUTM
% Function to create conversion of circuit coordinates in UTM to GPS coordinates
%
%% Syntax
%  createCircuitCoordConversionUTM();
%
%% Description
% Function is a user interface that allows the user to select the UTM zone 
% the circuit coordinates are currently in. The 
% conversion is created for getting GPS Lat/Lon for the OpenDSS bus 
% coordinates and the new Lat/Lon OpenDSS buscoords are saved.
%
%% Inputs
% * *none* 
%
%% Outputs
% * *none*  - a Circuit Conversion file is saved for any future plotting
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
% createCircuitCoordConversionUTM();
%

function varargout = createCircuitCoordConversionUTM(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @createCircuitCoordConversionUTM_OpeningFcn, ...
                   'gui_OutputFcn',  @createCircuitCoordConversionUTM_OutputFcn, ...
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


%% --- Executes just before createCircuitCoordConversionUTM is made visible.
function createCircuitCoordConversionUTM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to createCircuitCoordConversionUTM (see VARARGIN)

% Choose default command line output for createCircuitCoordConversionUTM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Populate listboxes
set(handles.listbox1, 'String', [{'C'}; {'D'}; {'E'}; {'F'}; ...
     {'G'}; {'H'}; {'J'}; {'K'}; {'L'}; {'M'}; {'N'}; {'P'}; {'Q'}; ...
     {'R'}; {'S'}; {'T'}; {'U'}; {'V'}; {'W'}; {'X'};]);
 set(handles.listbox2, 'String', [{'01'}; {'02'}; {'03'}; {'04'}; {'05'}; {'06'}; ...
     {'07'}; {'08'}; {'09'}; {'10'}; {'11'}; {'12'}; {'13'}; {'14'}; {'15'}; {'16'}; ...
     {'17'}; {'18'}; {'19'}; {'20'}; {'21'}; {'22'}; {'23'}; {'24'}; {'25'}; {'26'}; ...
     {'27'}; {'28'}; {'29'}; {'30'}; {'31'}; {'32'}; {'33'}; {'34'}; {'35'}; {'36'}; ...
     {'37'}; {'38'}; {'39'}; {'40'}; {'41'}; {'42'}; {'43'}; {'44'}; {'45'}; {'46'}; ...
     {'47'}; {'48'}; {'49'}; {'50'}; {'51'}; {'52'}; {'53'}; {'54'}; {'55'}; {'56'}; ...
     {'57'}; {'58'}; {'59'}; {'60'};]);
 set(handles.listbox2, 'value', 1);
 set(handles.listbox1, 'value', 1);

end


%% --- Outputs from this function are returned to the command line.
function varargout = createCircuitCoordConversionUTM_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

end


%% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)

end


%% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


%% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
end


%% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text9u.
function text9u_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text9u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 9);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text10u.
function text10u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 10);
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text11u.
function text11u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 11);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text12u.
function text12u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 12);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text13u.
function text13u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 13);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text14u.
function text14u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 14);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text15u.
function text15u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 15);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text16u.
function text16u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 16);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text17u.
function text17u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 17);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text18u.
function text18u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 18);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text19u.
function text19u_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 17);
set(handles.listbox2, 'Value', 19);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text10t.
function text10t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 10);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text11t.
function text11t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 11);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text12t.
function text12t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 12);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text13t.
function text13t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 13);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text14t.
function text14t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 14);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text15t.
function text15t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 15);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text16t.
function text16t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 16);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text17t.
function text17t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 17);
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text18t.
function text18t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 18);
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text19t.
function text19t_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 16);
set(handles.listbox2, 'Value', 19);
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text10s.
function text10s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 10);
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text11s.
function text11s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 11);
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text12s.
function text12s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 12);
end

%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text13s.
function text13s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 13);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text14s.
function text14s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 14);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text15s.
function text15s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 15);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text16s.
function text16s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 16);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text17s.
function text17s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 17);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text18s.
function text18s_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 15);
set(handles.listbox2, 'Value', 18);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text11r.
function text11r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 11);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text12r.
function text12r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 12);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text13r.
function text13r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 13);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text14r.
function text14r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 14);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text15r.
function text15r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 15);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text16r.
function text16r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 16);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text17r.
function text17r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 17);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text18r.
function text18r_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 14);
set(handles.listbox2, 'Value', 18);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text11q.
function text11q_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 13);
set(handles.listbox2, 'Value', 11);
end


%% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text18q.
function text18q_ButtonDownFcn(hObject, eventdata, handles)
set(handles.listbox1, 'Value', 13);
set(handles.listbox2, 'Value', 18);
end


%% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

[fnmMain pnmMain] = uigetfile({'*.dss'}, 'Select the OpenDSS file with the circuit');
if ~(isequal(fnmMain, 0) || isequal(pnmMain, 0))
    location = cd;
    DSSText.command = sprintf('Compile (%s)',[pnmMain, fnmMain]);
    DSSText.command = sprintf('solve');
    cd(location);

    [fnm pnm] = uigetfile({'*.dss; *.txt; *.xls; *.csv'}, 'Select your bus coordinates file');
    if ~(isequal(fnm, 0) || isequal(pnm, 0))
        [SUCCESS,MESSAGE,MESSAGEID] = movefile([pnm fnm], [pnm 'BACKUP_' fnm], 'f');
        if ~SUCCESS
            disp(MESSAGE);
            waitfor(warndlg(sprintf('A backup of your coordinates file was not made (see command window for the message).\n\nPlease manually make a backup copy of your coordinates file before clicking ''OK'' to this message.'), 'No Backup Made!'));
        end
        
        %Get the current coordinates
        [busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);
        
        %Create cell array of list box contents
        numbers = get(handles.listbox2, 'string');
        letters = get(handles.listbox1, 'string');
        
        %Prepare input to utm2deg
        UTMzone = repmat([numbers{get(handles.listbox2, 'value')} ' ' letters{get(handles.listbox1, 'value')}], length(busCoordArray), 1);
        
        %Convert the coordinates
        [lat lon] = utm2deg(busCoordArray(:, 1), busCoordArray(:, 2), UTMzone);
        
        %Write to file
        %Could not find a way to avoid iteration
        fid = fopen([pnm fnm], 'w');
        for ii = 1:length(busCoordNames)
            fprintf(fid, '%s %3.15f %3.15f\r\n', busCoordNames{ii}, lon(ii), lat(ii));
        end
        fclose(fid);
        
        location = cd;
        DSSText.command = sprintf('Compile (%s)',[pnmMain, fnmMain]);
        DSSText.command = sprintf('solve');
        cd(location);
        
        msgbox('Conversion complete. Save was successful.', 'Done');
        uiwait;
        
        close(handles.figure1);
    end
    
end
end

    
