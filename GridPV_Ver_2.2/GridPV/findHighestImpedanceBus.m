%% findHighestImpedanceBus
% Finds the highest impedance bus from the substation
%
%% Syntax
%  [highestImpedance highestImpedanceBus] = findHighestImpedanceBus(DSSCircObj, requiredLineRating);
%  [highestImpedance highestImpedanceBus] = findHighestImpedanceBus(DSSCircObj, requiredLineRating, threePhase);
%
%% Description
% Function to find highest impedance bus from the substation.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|requiredLineRating|* - the minimum allowed conductor size (amps) line 
% rating for PV placement. A larger plant requires a higher required line 
% rating.  To not restrict the search algorithm, set this to zero.
% * *|threePhase|* - optional input, logical value for if the bus must be 3
% phase.  If the input is a logical true, only 3 phase buses will be
% returned.
%
%% Outputs
% * *|highestImpedance|* - impedance rating between fromBus to toBus
% * *|highestImpedanceBus|* - name of bus with highest impedance to the source bus
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
% Returns the bus names for the highest impedance bus in the circuit
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% [highestImpedance highestImpedanceBus] = findHighestImpedanceBus(DSSCircObj, 220)
%

function [highestImpedance highestImpedanceBus] = findHighestImpedanceBus(DSSCircObj, requiredLineRating, varargin)
%% Parse Inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addRequired('requiredLineRating',@(x)isnumeric(x) && all(x>=0));
p.addOptional('threePhase', 0,  @(x)isnumeric(x) && ((x==1) || (x==0)));

p.parse(DSSCircObj, requiredLineRating, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
% Define the text interface
DSSText = DSSCircObj.Text;

%% Get sequence impedances
oldSolutionMode = DSSCircObj.ActiveCircuit.Solution.ModeID; %save old solution mode

DSSText.command = 'solve mode=faultstudy';
DSSText.command = 'Export SeqZ';
    while isempty(DSSText.Result) end %wait for export to finish
file = DSSText.Result;
MyCSV = importdata(file);
delete(file);

DSSText.command = ['set mode=',oldSolutionMode];
DSSText.command = 'solve';

Z1 = MyCSV.data(:,6);
Z0 = MyCSV.data(:,7);
busNames = MyCSV.textdata(2:end,1);
Za = Z1+(Z0-Z1)/3;

%% find which buses can support the required line rating

Lines = getLineInfo(DSSCircObj);

condition = [Lines.lineRating]>requiredLineRating;
Lines = Lines(condition);

bus1 = regexprep({Lines.bus1},'(\.[0-9]+)',''); %take out the phase numbers on buses
bus2 = regexprep({Lines.bus2},'(\.[0-9]+)',''); %take out the phase numbers on buses
busesToKeep = unique([bus1 bus2]);

%% filter buses that are not active
Buses = getBusInfo(DSSCircObj,busesToKeep);
Buses = Buses([Buses.voltagePU]>0.6);
busesToKeep = {Buses.name};
    
    
%% if the threePhase argument is true, filter all buses that are not 3-phase
if nargin==3 && threePhase==1
    Buses = getBusInfo(DSSCircObj,busesToKeep);
    Buses = Buses([Buses.numPhases]==3);
    busesToKeep = {Buses.name};
end


%% figure out which buses from the sequence impedances meet the criteria
condition = ismember(regexprep(upper(busNames),'(\.[0-9]+)',''),upper(busesToKeep));
busNames = busNames(condition);
Za = Za(condition);

[highestImpedance index] = max(Za);
highestImpedanceBus = busNames(index);

end
