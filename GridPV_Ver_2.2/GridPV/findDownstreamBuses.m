%% findDownstreamBuses
% Finds all buses downstream of the busName
%
%% Syntax
%  downstreamBuses = findDownstreamBuses(DSSCircObj,busName);
%  downstreamBuses = findDownstreamBuses(DSSCircObj,busName, _'PropertyName'_ ,PropertyValue);
%
%% Description
% Function to get all the bus names for buses that are downstream of the
% busName.  The downstream buses are defined as buses that are farther from the
% substation on the electrical path of busName.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|busName|* - string of the bus name to start search downstream
% * *|Properties|* - optional properties as one or more name-value pairs in any order
% * -- *|'Lines'|* - Structure of the circuit lines from getLineInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Transformers'|* - Structure of the circuit transformers from getTransformerInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
%
%% Outputs
% * *|downstreamBuses|* is a cell array of the bus names downstream from busName
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
% Returns downstream buses
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% downstreamBuses = findDownstreamBuses(DSSCircObj,'N292792')
%

function downstreamBuses = findDownstreamBuses(DSSCircObj,busName,varargin)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addRequired('busName', @ischar);
p.addParamValue('Lines', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));
p.addParamValue('Transformers', 'noInput', @(x)(ischar(x)&strcmp(x,'noInput'))|isstruct(x));

p.parse(DSSCircObj, busName, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

%% Grab circuit information
if ischar(Lines) && strcmp(Lines,'noInput')
    Lines = getLineInfo(DSSCircObj);
end
if ischar(Transformers) && strcmp(Transformers,'noInput')
    Transformers = getTransformerInfo(DSSCircObj);
end
linesBus1 = regexprep({Lines.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
linesBus2 = regexprep({Lines.bus2},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
transformerBus1 = regexprep({Transformers.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
transformerBus2 = regexprep({Transformers.bus2},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them

%% Initialize variables
busName = regexprep(busName,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
Buses = getBusInfo(DSSCircObj,{busName});

%% Initialize variables
downstreamBuses = {};
searchBusList = {lower(busName)}; %start search at selected bus
searchBusDistances = Buses.distance;

%% Continue looping one level further away (distance) until the end of each feeder branch is reached
while ~isempty(searchBusList)
    
    nextLevelBuses = {};
    nextLevelDistances = [];
    
    % For each bus in the search list, find the downstream buses connected
    % to it, and add those downstream buses to the next level bus 
    % list that will will become the searchBusList
    for ii=1:length(searchBusList)
        connectedBuses = {}; %store list of buses connected to searchBus
        connectedDistances = [];

        % Find all lines connected to the search bus               
        index = find(strcmpi(linesBus1,searchBusList{ii})==1);
        connectedBuses = [connectedBuses {Lines(index).bus2}];
        connectedDistances = [connectedDistances [Lines(index).bus2Distance]];

        index = find(strcmpi(linesBus2,searchBusList{ii})==1);
        connectedBuses = [connectedBuses {Lines(index).bus1}];
        connectedDistances = [connectedDistances [Lines(index).bus1Distance]];
    
        % Find all transformers connected to the bus
        index = find(strcmpi(transformerBus1,searchBusList{ii})==1);
        if ~isempty(index) % Make sure we are not including high side of xfmr and its subsequent downtream buses
            index = index([Transformers(index).bus2kV] <= [Transformers(index).bus1kV]);
        end
        connectedBuses = [connectedBuses {Transformers(index).bus2}];
        connectedDistances = [connectedDistances [Transformers(index).bus2Distance]];

        index = find(strcmpi(transformerBus2,searchBusList{ii})==1);
        if ~isempty(index) % Make sure we are not including high side of xfmr and its subsequent downtream buses
            index = index([Transformers(index).bus1kV] <= [Transformers(index).bus2kV]);
        end
        connectedBuses = [connectedBuses {Transformers(index).bus1}];
        connectedDistances = [connectedDistances [Transformers(index).bus1Distance]];
        
        % find all connected bus further away (distance)
        connectedBuses = connectedBuses(connectedDistances>=searchBusDistances(ii));
        connectedDistances = connectedDistances(connectedDistances>=searchBusDistances(ii));
        connectedBuses = regexprep(connectedBuses,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        
        % don't keep the connected bus if it was already found previously
        if ~isempty(connectedBuses)
            condition = ~(ismember(connectedBuses,downstreamBuses) | ismember(connectedBuses,searchBusList));
            connectedDistances = connectedDistances(condition);
            connectedBuses = connectedBuses(condition);
        end
        
        if ~isempty(connectedBuses)
            % remove any duplicates
            [connectedBuses,ia,ic] = unique(connectedBuses);
            connectedDistances = connectedDistances(ia);
            
            % save all further connected buses to the next level buses
            nextLevelBuses = [nextLevelBuses connectedBuses];
            nextLevelDistances = [nextLevelDistances connectedDistances];
        end
    end
    
    % save bus to downstreamBuses
    downstreamBuses = [downstreamBuses searchBusList];
    
    % Set the next level buses to the current search buses
    searchBusList = nextLevelBuses;
    searchBusDistances = nextLevelDistances;
    
end
downstreamBuses = downstreamBuses';

end



