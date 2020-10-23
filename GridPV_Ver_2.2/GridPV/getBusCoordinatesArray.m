%% getBusCoordinatesArray
% Gets the coordinates for all buses that have a location in OpenDSS
%
%% Syntax
%  [busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);
%
%% Description
% Function to get the buses and their coordinates for all buses that have a
% location in OpenDSS.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
%
%% Outputs
% * *|busCoordNames|* is the array of the bus names
% * *|busCoordArray|* is the matrix of bus coordinates (X,Y) corresponding
% to the bus name in busCoordNames.
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
% Returns the bus names and coordinates for the active circuit in OpenDSS
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% [busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);
% size(busCoordArray)
%

function [busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);

p.parse(DSSCircObj); %parse inputs

%%
% Define the text interface
DSSText = DSSCircObj.Text;

DSSText.Command = 'Export Buscoords';
fid = fopen(DSSText.Result);
busCoordStruct = textscan(fid, '%s %f %f', 'delimiter', ',');
fclose(fid);
busCoordNames = busCoordStruct{1};
busCoordArray = [busCoordStruct{2} busCoordStruct{3}];

delete(DSSText.Result);

if all(busCoordArray==0)
    error('Unable to find any buses with coordinates. Make sure Buscoords is defined in OpenDSS.');
end

end