%% findLongestDistanceBus
% Finds the bus for each phase that is farthest distance from the source bus
%
%% Syntax
%  [longestDistance toBus] = findLongestDistanceBus(DSSCircObj, phaseOption);
%
%% Description
% Function to find the bus for each phase that is farthest distance from the source bus.
% This can be run to find the farthest bus for each phase (generally single
% phase) or farthest 3 phase bus.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|phaseOption|* - 'perPhase' for the farthest bus on each phase or
% '3phase' for the farthest 3 phase bus
%
%% Outputs
% * *|longestDistance|* - distance between fromBus to toBus
% * *|toBus|* - name of bus with highest impedance to the energy monitor
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
% Returns the bus names and distance for the farthest bus
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% [longestDistance toBus] = findLongestDistanceBus(DSSCircObj, 'perPhase')
%

function [longestDistance toBus] = findLongestDistanceBus(DSSCircObj, phaseOption)
%% Parse Inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addRequired('phaseOption',@(x) ischar(x) && (strcmpi(x, 'perPhase') || strcmpi(x, '3phase')));

p.parse(DSSCircObj, phaseOption); %parse inputs

% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
% Define the text interface
DSSText = DSSCircObj.Text;

%% If 3 phase is required, filter all non-3-phase buses
if strcmpi(phaseOption,'3phase')
    busNames = DSSCircuit.AllBusNames;
    Buses = getBusInfo(DSSCircObj,busNames);
    Buses = Buses([Buses.numPhases]==3);
    Buses = Buses([Buses.voltage]>50);
    distances = [Buses.distance];
    busNames = {Buses.name};
    [longestDistance index] = max(distances);
    toBus = busNames{index};
else
    distances = DSSCircuit.AllNodeDistancesByPhase(1);
    busNames = DSSCircuit.AllNodeNamesByPhase(1);
    busVoltages = DSSCircuit.AllNodeVmagByPhase(1);
    busNames = busNames(busVoltages>50); distances = distances(busVoltages>50);
    [longestDistance(1) index] = max(distances);
    toBus(1) = busNames(index);
    
    distances = DSSCircuit.AllNodeDistancesByPhase(2);
    busNames = DSSCircuit.AllNodeNamesByPhase(2);
    busVoltages = DSSCircuit.AllNodeVmagByPhase(2);
    busNames = busNames(busVoltages>50); distances = distances(busVoltages>50);
    [longestDistance(2) index] = max(distances);
    toBus(2) = busNames(index);
    
    distances = DSSCircuit.AllNodeDistancesByPhase(3);
    busNames = DSSCircuit.AllNodeNamesByPhase(3);
    busVoltages = DSSCircuit.AllNodeVmagByPhase(3);
    busNames = busNames(busVoltages>50); distances = distances(busVoltages>50);
    [longestDistance(3) index] = max(distances);
    toBus(3) = busNames(index);
end

end
