%% findSubstationLocation
% Locates the substation coordinates
%
%% Syntax
%  coordinates = findSubstationLocation(DSSCircObj);
%
%% Description
% Function to find the coordinates of the substation.  This is used for
% plotting the substation on circuit diagrams.  The substation is located
% at the bus coordinate with the shortest "distance".
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
%
%% Outputs
% * *|coordinates|* is the [Y X] coordinates for the substation bus
% location
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
% Returns the substation location
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% coordinates = findSubstationLocation(DSSCircObj)
%

function coordinates = findSubstationLocation(DSSCircObj)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);

p.parse(DSSCircObj); %parse inputs

%%
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
% Define the text interface
DSSText = DSSCircObj.Text;

busNames = DSSCircuit.AllBusNames;
busDistances = DSSCircuit.AllBusDistances;

if all(busDistances==0)
    error('All bus distances are equal to zero.  Make sure the circuit contains an energy meter at the substation.');
end

% only keep buses that have coordinates
[busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);
condition = ismember(upper(busNames),upper(busCoordNames));
busNames = busNames(condition);
busDistances = busDistances(condition);

if isempty(busDistances)
    error('Unable to find any buses with coordinates. Make sure Buscoords is defined in OpenDSS.');
end

% find minimum distance bus
[distance index] = min(busDistances);
minBusName = busNames{index};

% find coordinates of minimum distance bus
condition = ismember(upper(busCoordNames),upper(minBusName));
coordinates = fliplr(busCoordArray(condition,:));

end



