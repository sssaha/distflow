%% getCoordinates
% Gets the coordinates for the buses in busNames
%
%% Syntax
%  coordinates = getCoordinates(DSSCircObj);
%  coordinates = getCoordinates(DSSCircObj,busNames);
%
%% Description
% Function to get coordinates for the buses in busNames. If optional input busNames
% contains a cell array, the function will return a structure for each busName, otherwise coordinates will
% contain all buses in the circuit.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|busNames|* - optional cell array of bus names to find locations for
%
%% Outputs
% * *|coordinates|* is the array of bus coordinates corresponding to
% busNames. The first column is the y values, and second column is x values
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
% Returns the coordinates for buses
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% coordinates = getCoordinates(DSSCircObj); %Get all bus coordinates
% coordinates = getCoordinates(DSSCircObj,{'N1311915'}) %Get coordinates for bus N1311915
% coordinates = getCoordinates(DSSCircObj,[{'N1311915'}; {'n284022'}]) %Get coordinates for two buses
%

function coordinates = getCoordinates(DSSCircObj,varargin)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addOptional('busNames', 'noInput', @iscellstr);

p.parse(DSSCircObj, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

if strcmp(busNames,'noInput')
    busNames = DSSCircObj.ActiveCircuit.AllBusNames;
end

%%
% Define the text interface
DSSText = DSSCircObj.Text;

[busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);

busNames = (regexprep(busNames,'(\.[0-9]+)','')); %take out the phase numbers on buses if they have them

[tf, index] = ismember(upper(busNames),upper(busCoordNames));

index = index(index~=0);

if any(tf)
    xCoord = busCoordArray(index,1);
    yCoord = busCoordArray(index,2);
    coordinates = [yCoord xCoord];
else
    error('None of the bus names entered match a bus in the circuit. Please double check your input.');
end