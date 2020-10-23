%% findUpstreamBuses
% Finds all buses upstream of the busName
%
%% Syntax
%  upstreamBuses = findUpstreamBuses(DSSCircObj,busName);
%  upstreamBuses = findUpstreamBuses(DSSCircObj,busName, _'PropertyName'_ ,PropertyValue);
%
%% Description
% Function to get all the bus names for buses that are upstream of the
% busName.  The upstream buses are defined as buses that are closer to the
% substation on the electrical path to busName.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|busName|* - string of the bus name to start search upstream
% * *|Properties|* - optional properties as one or more name-value pairs in any order
% * -- *|'Lines'|* - Structure of the circuit lines from getLineInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
% * -- *|'Transformers'|* - Structure of the circuit transformers from getTransformerInfo.  If no input is given, the structure is filled from the most current power flow solution in DSSCircObj COM.
%
%% Outputs
% * *|upstreamBuses|* is a cell array of the bus names upstream from busName
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
% Returns upstream buses
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% upstreamBuses = findUpstreamBuses(DSSCircObj,'n292286')
%

function upstreamBuses = findUpstreamBuses(DSSCircObj,busName,varargin)
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

inputBus = getBusInfo(DSSCircObj, {busName});
if inputBus.distance == 0
    error('Please choose a bus that is a non-zero distance from the substation.')
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
upstreamBuses = {busName};
searchBus = busName; %start search at selected bus
searchBus = regexprep(searchBus,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
minDistance=1e9;

%% Continue looping closer (distance) until the head of the feeder is reached
while minDistance>0
    connectedBuses = {}; %store list of buses connected to searchBus
    connectedDistances = [];
    
    % Find all lines connected to the bus
    index = find(strcmpi(upper(linesBus1),upper(searchBus))==1);
    connectedBuses = [connectedBuses {Lines(index).bus2}];
    connectedDistances = [connectedDistances [Lines(index).bus2Distance]];
    
    index = find(strcmpi(upper(linesBus2),upper(searchBus))==1);
    connectedBuses = [connectedBuses {Lines(index).bus1}];
    connectedDistances = [connectedDistances [Lines(index).bus1Distance]];
    
    % Find all transformers connected to the bus
    index = find(strcmpi(upper(transformerBus1),upper(searchBus))==1);
    connectedBuses = [connectedBuses {Transformers(index).bus2}];
    connectedDistances = [connectedDistances [Transformers(index).bus2Distance]];
    
    index = find(strcmp(upper(transformerBus2),upper(searchBus))==1);
    connectedBuses = [connectedBuses {Transformers(index).bus1}];
    connectedDistances = [connectedDistances [Transformers(index).bus1Distance]];
    
    %remove downstream buses that are already in the list
    connectedBuses = regexprep(connectedBuses,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    notFound = ~ismember(connectedBuses,upstreamBuses);
    connectedBuses = connectedBuses(notFound);
    connectedDistances = connectedDistances(notFound);
    
    % check if any connections found
    if isempty(connectedBuses)
        error(sprintf('Did not find anything connected to bus %s.  Please check that bus %s is in the compiled circuit.',searchBus, searchBus))
    end
    
    % find closest connected bus (distance)
    [minDistance busIndex] = min(connectedDistances);
    if sum(minDistance==connectedDistances)>1 %if there is more than one bus with the same distance away, check to see if one of those has a higher voltage
        Buses = getBusInfo(DSSCircObj,connectedBuses(minDistance==connectedDistances));
        [B,IX] = sort([Buses.kVBase],'descend');
        busIndex = find(strcmp(connectedBuses,Buses(IX(1)).name));
    end
    
    % set closest bus as new search bus
    searchBus = connectedBuses{busIndex};
    searchBus = regexprep(searchBus,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    
    % save bus to upstreamBuses
    upstreamBuses = [upstreamBuses searchBus];
end


end



