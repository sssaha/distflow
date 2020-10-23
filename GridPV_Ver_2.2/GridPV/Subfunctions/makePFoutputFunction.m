%% makePFoutputFunction
% GUI for creating power factor as a function of PV power output
%
%% Syntax
%  makePFoutputFunction()
%
%% Description
% This function is a user interface to create the Power Factor
% as a function of PV power output.  The user draws the function
% and then saves it to a .mat file.  This function is often called
% from placePVplant.m when the PV plant power factor control is selected.
% The saved mat file is used in createPVscenarioFiles.m when the solar scenario OpenDSS
% generators are created.
%
%% Inputs
% * *none*
%
%% Outputs
% * *none* saves a *.mat file with the power factor function of PV power output
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
% makePFoutputFunction()
%

function varargout = makePFoutputFunction(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @makePFoutputFunction_OpeningFcn, ...
                   'gui_OutputFcn',  @makePFoutputFunction_OutputFcn, ...
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

%% --- Executes just before makePFoutputFunction is made visible.
function makePFoutputFunction_OpeningFcn(hObject, eventdata, handles, varargin)
global powerFactor MWoutput

powerFactor = [1,1,1,1,1,1,0.975,0.975,0.975,0.95,0.95];
MWoutput = (0:0.1:1)*100;

plotLine();
% Choose default command line output for makePFoutputFunction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

end


%% --- Outputs from this function are returned to the command line.
function varargout = makePFoutputFunction_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
end


%% --- Executes on button press in plus1.
function plus1_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(1) = powerFactor(1)+0.005;
plotLine();

end

%% --- Executes on button press in plus2.
function plus2_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(2) = powerFactor(2)+0.005;
plotLine();

end


%% --- Executes on button press in plus5.
function plus5_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(5) = powerFactor(5)+0.005;
plotLine();

end


%% --- Executes on button press in plus8.
function plus8_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(8) = powerFactor(8)+0.005;
plotLine();

end


%% --- Executes on button press in plus9.
function plus9_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(9) = powerFactor(9)+0.005;
plotLine();

end


%% --- Executes on button press in plus10.
function plus10_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(10) = powerFactor(10)+0.005;
plotLine();

end


%% --- Executes on button press in plus11.
function plus11_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(11) = powerFactor(11)+0.005;
plotLine();

end


%% --- Executes on button press in plus12.
function plus12_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(12) = powerFactor(12)+0.005;
plotLine();

end


%% --- Executes on button press in plus13.
function plus13_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(13) = powerFactor(13)+0.005;
plotLine();

end


%% --- Executes on button press in plus14.
function plus14_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(14) = powerFactor(14)+0.005;
plotLine();

end


%% --- Executes on button press in plus15.
function plus15_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(15) = powerFactor(15)+0.005;
plotLine();

end


%% --- Executes on button press in minus1.
function minus1_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(1) = powerFactor(1)-0.005;
plotLine();

end


%% --- Executes on button press in minus2.
function minus2_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(2) = powerFactor(2)-0.005;
plotLine();

end


%% --- Executes on button press in minus5.
function minus5_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(5) = powerFactor(5)-0.005;
plotLine();

end


%% --- Executes on button press in minus8.
function minus8_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(8) = powerFactor(8)-0.005;
plotLine();

end


%% --- Executes on button press in minus9.
function minus9_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(9) = powerFactor(9)-0.005;
plotLine();

end


%% --- Executes on button press in minus10.
function minus10_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(10) = powerFactor(10)-0.005;
plotLine();

end


%% --- Executes on button press in minus11.
function minus11_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(11) = powerFactor(11)-0.005;
plotLine();

end


%% --- Executes on button press in minus12.
function minus12_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(12) = powerFactor(12)-0.005;
plotLine();

end


%% --- Executes on button press in minus13.
function minus13_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(13) = powerFactor(13)-0.005;
plotLine();

end


%% --- Executes on button press in minus14.
function minus14_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(14) = powerFactor(14)-0.005;
plotLine();

end


%% --- Executes on button press in minus15.
function minus15_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(15) = powerFactor(15)-0.005;
plotLine();

end


%% --- Executes on button press in plus3.
function plus3_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(3) = powerFactor(3)+0.005;
plotLine();

end


%% --- Executes on button press in minus3.
function minus3_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(3) = powerFactor(3)-0.005;
plotLine();

end


%% --- Executes on button press in plus4.
function plus4_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(4) = powerFactor(4)+0.005;
plotLine();

end


%% --- Executes on button press in minus4.
function minus4_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(4) = powerFactor(4)-0.005;
plotLine();

end


%% --- Executes on button press in plus6.
function plus6_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(6) = powerFactor(6)+0.005;
plotLine();

end


%% --- Executes on button press in minus6.
function minus6_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(6) = powerFactor(6)-0.005;
plotLine();

end


%% --- Executes on button press in plus7.
function plus7_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(7) = powerFactor(7)+0.005;
plotLine();

end


%% --- Executes on button press in minus7.
function minus7_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

powerFactor(7) = powerFactor(7)-0.005;
plotLine();

end

%%
function plotLine
global powerFactor MWoutput

powerFactor(powerFactor>1) = 1;
interpMWoutput = 0:0.01:100;

powerFactorExpand = interp1(MWoutput,powerFactor,interpMWoutput,'nearest');

hold off;
plot(interpMWoutput,powerFactorExpand,'LineWidth',2);
hold all;
plot(MWoutput,powerFactor,'LineStyle','none','Marker','square','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','b');
xlim([min(MWoutput) max(MWoutput)]);
set(gca,'XTick',MWoutput);
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('PV Output (% of rated)');
ylabel('Power Factor');
end


%% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global powerFactor MWoutput

percentPVoutput = MWoutput/100;

[file,path] = uiputfile('*.mat','Save Power Factor Function As','PowerFactorFunction.mat');
save([path,file],'percentPVoutput','powerFactor');

end
