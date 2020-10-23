%% circuitCheck
% Used to error-check the circuit for any obvious abnormalities.
%
%% Syntax
%  warnSt = circuitCheck(DSSCircObj);
%  warnSt = circuitCheck(DSSCircObj,'Warnings','off');
%
%% Description
% Used for checking OpenDSS circuits for errors or abnormalities that do
% not prevent OpenDSS from running but will cause errors during analysis
% (e.g. Phase-a line downstream of a bus with only phases b and c). It is
% capable of performing a complete circuit check with a warning describing
% each error found. Warnings can be turned off. A more comprehensive list of
% elements that cause the errors can be found inside the structure, warnSt,
% that is outputted at the end of the check.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|'Warnings'|* - indicates if the user wants command-prompt
% warnings on or not |{'on'} | 'off'|
%
%% Outputs
% * *|warnSt|* is a structure with parameters relating to the results of
% various validity check. If the circuit failed a check, an entry for that
% check appears in this structure with fields for the check name, a string
% with the description, and a list of offenders that caused the fail.
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
% Example of a circuit test:
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% warnSt = circuitCheck(DSSCircObj)
% warnSt = circuitCheck(DSSCircObj, 'Warnings', 'off');
%

function warnSt = circuitCheck(DSSCircObj, varargin)
%% Parse Inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS_internal);
p.addParamValue('Warnings', 'on', @(x)any(strcmpi(x,{'on','off'})));

p.parse(DSSCircObj, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

%% Define Circuit

% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
% Define the text interface
DSSText = DSSCircObj.Text;

%% Begin Check
warnSt = struct();

%% ************** Define all check thresholds **********************
lineLengthChk = 5; % km max line length
busDistChk = 35; % km max distance for buses
elmntVoltWndw = 0.05; % 0.05 equates to (+/- 5 percent)
busVoltWndw = 0.05; % 0.05 equates to (+/- 5 percent)
maxLineLoad = 1.0; % 100 percent line loading threshold
linePercDiff = 1.3; % A downstream line rating greater than 130% of its upstream line gets flagged
xfmrPowerMax = 1.0; %if total of power through a transformer is greater than 100% of the transformer rated kva, flag
xfmrPowerMin = 0.01; %if total of power through a transformer is less than 1% of the transformer rated kva, flag

%% Get Isolated Loads and Branches
isoBranches = DSSCircuit.Topology.AllIsolatedBranches;
isoLoads = DSSCircuit.Topology.AllIsolatedLoads; %There are two options for isolated loads - they are disabled, or they are not really isolated (they have voltage). If a load is isolated and enabled, the solution will not converge
if strcmp(isoBranches{1}, 'NONE')
    isoBranches = [];
end
if strcmp(isoLoads{1}, 'NONE')
    isoLoads = [];
end
% Combine all isolated
allIsolated = [isoBranches; isoLoads];
% Remove the blanks inserted by the COM
blanks = strcmp(allIsolated, '');
allIsolated(blanks) = [];
if ~isempty(allIsolated)
    % See if any of the branches are enabled
    notIsolated = false(length(allIsolated), 1); % Flag for cataching disabled elements that are present in the COM, but not actually isolated 
    enabled = false(length(allIsolated), 1);
    hasVolt = false(length(allIsolated), 1);
    for ii = 1:length(allIsolated)
       DSSCircuit.SetActiveElement(allIsolated{ii});
       enabled(ii) = DSSCircuit.ActiveElement.Enabled;
       if enabled(ii) % check for voltage
           voltages = DSSCircuit.ActiveElement.Voltages;
           if ~any(isnan(voltages)) && any(voltages~=0)
               hasVolt(ii) = 1;
           end
       else
           % Check if the element is a disabled element attached to an
           % enabled part of the circuit
           busNames = DSSCircuit.ActiveElement.BusNames;
           busNames = regexprep(busNames,'(\.[0-9]+)','');
           notIsolated(ii) = any(ismember(busNames, DSSCircuit.AllBusNames));
       end
    end
    allIsolated = allIsolated(~notIsolated);
    enabled = enabled(~notIsolated);
    hasVolt = hasVolt(~notIsolated);
    disabIsolated = allIsolated(~enabled);
    enabVoltIsolated = allIsolated(hasVolt);
    enabNoVoltIsolated = allIsolated(enabled & ~hasVolt);
    if ~isempty(allIsolated)
        warnSt(1).IsolatedElem.nm = 'IsolatedElements';
        warnSt(1).IsolatedElem.str = sprintf('There are %d isolated elements, %d of which are enabled and without voltage. ', ...
             length(allIsolated), length(enabNoVoltIsolated));
        % prepare cell array for output
        warnSt(1).IsolatedElem.offenders = cell(length(allIsolated), 4);
        % output offenders
        warnSt(1).IsolatedElem.offenders(1:length(enabNoVoltIsolated)+1, 1) = ...
            [{'Isolated&Enabled w/out Voltage'}; enabNoVoltIsolated];
        warnSt(1).IsolatedElem.offenders(1:length(disabIsolated)+1, 2) = ...
            [{'Isolated&Disabled'}; disabIsolated];
        warnSt(1).IsolatedElem.offenders(1:length(enabVoltIsolated)+1, 3) = ...
            [{'Isolated&Enabled w/ Voltage'}; enabVoltIsolated];
        warnSt(1).IsolatedElem.offenders(1:length(allIsolated)+1, 4) = ...
            [{'AllIsolated'}; allIsolated];
    end
end

%% Isolated Nodes
% Only check if a valid solution exists
if DSSCircObj.ActiveCircuit.Solution.Converged == 1
    
    IsolatedNodes = DSSCircuit.AllNodeNames(isnan(DSSCircuit.AllBusVmag));
    if ~isempty(IsolatedNodes)
        warnSt(1).IsolatedNodes.nm = 'IsolatedNodes';
        warnSt(1).IsolatedNodes.str = sprintf('There are %d isolated nodes in the circuit. Investigate the node by using DSSText.command=''show busflow BUSNAME kVA elem''; where BUSNAME is the name of the isolated node without the decimal phasing.', length(IsolatedNodes));
        warnSt(1).IsolatedNodes.offenders = [{'NodeName'}; IsolatedNodes];
    end

end

%% Loops
numLoops = DSSCircuit.Topology.NumLoops;
if numLoops
    warnSt(1).Loops.nm = 'Loops';
    warnSt(1).Loops.str = sprintf('There are %d loops in the circuit. The loops can also be viewed using DSSText.command=''show loops'';', numLoops);
    warnSt(1).Loops.offenders = [{'AllLoopedPairs'}; DSSCircuit.Topology.AllLoopedPairs];
end

%% Get Element Info
% Using internal functions to only pull fields we need, saving time

%Lines
    % The Lines internal function includes the first circuit check: Line names and matching bus phases
    [Lines warnStLines warnFlag] = getLineInfo_Internal(DSSCircObj);
    if warnFlag
        warnSt(1).InvalidLineBusName = warnStLines;
    end
%Transformers
    Transformers = getTransformerInfo_Internal(DSSCircObj);
%Buses
    Buses = getBusInfo_Internal(DSSCircObj, DSSCircObj.ActiveCircuit.AllBusNames);
%Loads
    Loads  = getLoadInfo_Internal(DSSCircObj);
%Capacitors
    Capacitors = getCapacitorInfo_Internal(DSSCircObj);
%Generators   
    Generators = getGeneratorInfo_Internal(DSSCircObj);
%PV
    PV = getPVInfo_Internal(DSSCircObj);
    
%% Check for Bus Coordinates

%BUSCOORDS
busCoordArray = getBusCoordinatesArray_Internal(DSSCircObj);
if all(busCoordArray==0)
    warnSt(1).NoBusCoords.nm = 'NoBusCoords';
    warnSt(1).NoBusCoords.str = 'There are no bus coordinates with this compiled circuit. Toolbox functionality will be severely limited.';
else
    NoCoordBuses = Buses(sum(reshape([Buses.coordinates],2,[]))==0);
    NoCoordPrimaryBuses = NoCoordBuses([NoCoordBuses.voltage]>600);
    if ~isempty(NoCoordPrimaryBuses)
        warnSt(1).MissingBusCoords.nm = 'MissingBusCoords';
        warnSt(1).MissingBusCoords.str = sprintf('There are %d buses above 600V that are missing coordinates.', length(NoCoordPrimaryBuses));
        warnSt(1).MissingBusCoords.offenders = cell(length(NoCoordPrimaryBuses)+1, 1);
        warnSt(1).MissingBusCoords.offenders(1, 1) = [{'BusName'}];
        warnSt(1).MissingBusCoords.offenders(2:end, 1) = {NoCoordPrimaryBuses.name}';
    end
    
    if ~all(busCoordArray(:,2)>-90 & busCoordArray(:,2)<90 & busCoordArray(:,1)>-180 & busCoordArray(:,1)<180)
        warnSt(1).BusCoordLatLon.nm = 'BusCoordLatLon';
        warnSt(1).BusCoordLatLon.str = sprintf('The bus coordinates are not in latitude/longitude values. The mapping background will not be able to display for the circuit plot. Use initCoordConversion(); to convert coordinates to Lat/Lon.');
    end
end
    
%% Check parameters

% LINE LENGTHS
% Filter by bus distance difference instead of length because of
% consistency of units (km) - line lengths can be anything
badLines = Lines(abs([Lines.bus2Distance] - [Lines.bus1Distance])> lineLengthChk);
if ~isempty(badLines)
    warnSt(1).LineLength.nm = 'LineLength';
    warnSt(1).LineLength.str = sprintf('%d of the line lengths exceed %d kilometers.', length({badLines.name}), lineLengthChk);
    warnSt(1).LineLength.offenders = cell(length(badLines)+1, 2);
    warnSt(1).LineLength.offenders(1, :) = [{'LineName'}, {'Line Length (km)'}];
    warnSt(1).LineLength.offenders(2:end, 1) = {badLines.name}';
    warnSt(1).LineLength.offenders(2:end, 2) = num2cell(abs([badLines.bus2Distance] - [badLines.bus1Distance])');
end

% LINE RATINGS
badLinei = 1; %initialize offenders structure counter
% iterate through Lines
for ii = 1:length(Lines)
    % If it is overloaded greater than thershold, store it and increment our counter
    if Lines(ii).bus1Current/Lines(ii).lineRating > maxLineLoad
        badLines(badLinei) = Lines(ii);
        badLinei = badLinei + 1;
    end
end
%if some were found, create entry in the warning structure
if badLinei > 1
    [junk sortIndex] = sort([100*[badLines.bus1Current]./[badLines.lineRating]]','descend');
    warnSt(1).LineOverLoading.nm = 'LineOverLoading';
    warnSt(1).LineOverLoading.str = sprintf('%d of the lines are loaded more than %d percent. Visualize using plotCircuitLines(DSSCircObj,''Coloring'',''lineLoading'')', length({badLines.name}), maxLineLoad*100);
    warnSt(1).LineOverLoading.offenders = cell(length(badLines)+1, 2);
    warnSt(1).LineOverLoading.offenders(1, :) = [{'LineName'}, {'Line Loading (%)'}];
    warnSt(1).LineOverLoading.offenders(2:end, 1) = {badLines(sortIndex).name}';
    warnSt(1).LineOverLoading.offenders(2:end, 2) = num2cell([100*[badLines(sortIndex).bus1Current]./[badLines(sortIndex).lineRating]]');
end

% BUS DISTANCES
% Check to see if any buses are greater than the threshold away from the substation, which may indicate a problem with the circuit set up
% Logical AND with a check for > 0.05 to filter out any buses that may be in a disabled section of the circuit
badBuses = Buses([Buses.distance] > busDistChk & [Buses.voltagePU] > 0.05);
if ~isempty(badBuses)
    warnSt(1).BusDistance.nm = 'BusDistance';
    warnSt(1).BusDistance.str = sprintf('%d of the bus distances exceeds %d km from the substation.', length({badBuses.name}), busDistChk);
    warnSt(1).BusDistance.offenders = cell(length(badBuses)+1, 2);
    warnSt(1).BusDistance.offenders(1, :) = [{'Bus Name'}, {'Bus Distance (km)'}];
    warnSt(1).BusDistance.offenders(2:end, 1) = {badBuses.name}';
    warnSt(1).BusDistance.offenders(2:end, 2) = num2cell([badBuses.distance]');
    warnSt(1).BusDistance.offenders = {badBuses.name}';
end

%% Check Element Voltages and Phases
% Only check if a valid solution exists
if DSSCircObj.ActiveCircuit.Solution.Converged == 1
    % CAPACITORS
    if ~isempty(Capacitors) && ~strcmp(Capacitors(1).name, 'NONE')
        Capi = 1; % initialize our offenders counter
        Busi = 1; % initialize our offenders' buses counter
        % Iterate over all capacitor elements
        for ii = 1:length(Capacitors)
            capBusTmp = Buses(strcmpi(regexprep(Capacitors(ii).busName,'(\.[0-9]+)',''), {Buses.name})); %temporary bus to be save if element is bad
            if Capacitors(ii).numPhases >= 2 % determine if kV needs to be scaled to l-n value
                LL_LN = sqrt(3); % if it is 3 phase, convert l-l value to l-n
            else
                if Capacitors(ii).isDelta
                    LL_LN = sqrt(3); % even if it is not 3 phase but it is delta connection, convert
                else
                    LL_LN = 1; % if it is 1 phase and not delta connection, leave it alone
                end
            end
            % If our current capacitor is outside of the capacitor's bus's rating - save the capacitor, save the bus, and increment both counters
            if abs(capBusTmp.kVBase - Capacitors(ii).kV/LL_LN) > elmntVoltWndw*capBusTmp.kVBase && capBusTmp.voltagePU ~= 0
                badCap(Capi) = Capacitors(ii);
                badCap(Capi).kV = badCap(Capi).kV/LL_LN;
                Capi = Capi + 1;
                capBus(Busi) = capBusTmp;
                Busi = Busi + 1;
            end
        end
        % If some offenders were found, create entry in warning structure
        if Capi > 1
            warnSt(1).CapacitorRatingMismatch.nm = 'CapacitorRatingMismatch';
            warnSt(1).CapacitorRatingMismatch.str = sprintf('%d of the capacitor kV ratings differ from its bus kV rating by more than %d percent.', length({badCap.name}), elmntVoltWndw*100);
            warnSt(1).CapacitorRatingMismatch.offenders = cell(Capi, 4);
            warnSt(1).CapacitorRatingMismatch.offenders(1, :) = [{'Capacitor Name'} {'Capacitor Rating (kV l-l)'} {'Bus Name'} {'Bus Rating (kV l-l)'}];
            warnSt(1).CapacitorRatingMismatch.offenders(2:end, 1) = {badCap.name}';
            warnSt(1).CapacitorRatingMismatch.offenders(2:end, 2) = num2cell([badCap.kV]'*sqrt(3));
            warnSt(1).CapacitorRatingMismatch.offenders(2:end, 3) = {capBus.name}';
            warnSt(1).CapacitorRatingMismatch.offenders(2:end, 4) = num2cell([capBus.kVBase]'*sqrt(3));
        end
    end

    % LOADS
    if ~isempty(Loads) && ~strcmp(Loads(1).name, 'NONE')
        Loadi = 1; % initialize our offenders counter
        Busi = 1; % initialize our offenders' buses counter
        for ii = 1:length(Loads)
            %temporary bus to be save if element is bad
            loadBusTmp = Buses(strcmpi(regexprep(Loads(ii).busName,'(\.[0-9]+)',''), {Buses.name})); %remove the decimal notation from load's bus name
            if Loads(ii).numPhases >= 2 %determine if kV needs to be scaled to l-n value
                LL_LN = sqrt(3); % if it is 3 phase, convert l-l value to l-n
            else
                if Loads(ii).isDelta
                    LL_LN = sqrt(3); % even if it is not 3 phase but it is delta connection, convert
                    if Loads(ii).deltaVoltOutOfPhase %180 out of phase, so 2 instead of 120 degrees and sqrt(3)
                        LL_LN = 2;
                    end
                else
                    LL_LN = 1; % if it is 1 phase, leave it alone
                end
            end
            % If our current load is outside of the load's bus's rating - save the load, save the bus, and increment both counters
            if abs(loadBusTmp.kVBase - Loads(ii).kV/LL_LN) > elmntVoltWndw*loadBusTmp.kVBase && loadBusTmp.voltagePU ~= 0
                badLoad(Loadi) = Loads(ii);
                badLoad(Loadi).kV = badLoad(Loadi).kV/LL_LN;
                Loadi = Loadi + 1;
                loadBus(Busi) = loadBusTmp;
                Busi = Busi + 1;
            end
        end
     %if some offenders were found, create entry in the warning structure
        if Loadi > 1
            warnSt(1).LoadRatingMismatch.nm = 'LoadRatingMismatch';
            warnSt(1).LoadRatingMismatch.str = sprintf('%d of the load kV ratings differ from its bus kV rating by more than %d percent.', length({badLoad.name}), elmntVoltWndw*100);
            warnSt(1).LoadRatingMismatch.offenders = cell(Loadi, 4);
            warnSt(1).LoadRatingMismatch.offenders(1, :) = [{'Load Name'} {'Load Rating (kV l-l)'} {'Bus Name'} {'Bus Rating (kV l-l)'}];
            warnSt(1).LoadRatingMismatch.offenders(2:end, 1) = {badLoad.name}';
            warnSt(1).LoadRatingMismatch.offenders(2:end, 2) = num2cell([badLoad.kV]'*sqrt(3));
            warnSt(1).LoadRatingMismatch.offenders(2:end, 3) = {loadBus.name}';
            warnSt(1).LoadRatingMismatch.offenders(2:end, 4) = num2cell([loadBus.kVBase]'*sqrt(3));
        end
    end

    %GENERATORS
    if ~strcmp(Generators(1).name, 'NONE')
        PVi = 1;  % initialize our offenders counter
        Busi = 1;  % initialize our offenders' buses counter
        for ii = 1:length(Generators)
            pvBusTmp = Buses(strcmp(regexprep(Generators(ii).busName,'(\.[0-9]+)',''), {Buses.name})); %temp bus to save if element is bad
            
            %determine if kV needs to be scaled to l-n value
            if Generators(ii).numPhases >= 2 
                LL_LN = sqrt(3); % if it is 3 phase, convert l-l value to l-n
            else
                if Generators(ii).isDelta
                   LL_LN = sqrt(3); % even if it is not 3 phase but it is delta connection, convert
                   if Generators(ii).deltaVoltOutOfPhase %180 out of phase, so 2 instead of 120 degrees and sqrt(3)
                        LL_LN = 2;
                   end
                else
                    LL_LN = 1; % if it is 1 phase, leave it alone
                end
            end
            
            % If our current capacitor is outside of the Generators's bus's rating - save the Generators, save the bus, and increment both counters
            if abs(pvBusTmp.kVBase - Generators(ii).kV/LL_LN) > elmntVoltWndw*pvBusTmp.kVBase && pvBusTmp.voltagePU ~= 0
                badPV(PVi) = PV(ii);
                badPV(PVi).kV = badPV(PVi).kV/LL_LN;
                PVi = PVi + 1;
                pvBus(Busi) = pvBusTmp;
                Busi = Busi + 1;
            end
        end
     %if some offenders were found, create entry in the warning structure
        if PVi > 1
            warnSt(1).GeneratorRatingMismatch.nm = 'GeneratorRatingMismatch';
            warnSt(1).GeneratorRatingMismatch.str = sprintf('%d of the Generator kV ratings differ from its bus kV rating by more than %d percent.', length({badPV.name}), elmntVoltWndw*100);
            warnSt(1).GeneratorRatingMismatch.offenders = cell(PVi, 4);
            warnSt(1).GeneratorRatingMismatch.offenders(1, :) = [{'Generator Name'} {'Generator Rating (kV l-l)'} {'Bus Name'} {'Bus Rating(kV l-l)'}];
            warnSt(1).GeneratorRatingMismatch.offenders(2:end, 1) = {badPV.name}';
            warnSt(1).GeneratorRatingMismatch.offenders(2:end, 2) = num2cell([badPV.kV]'*sqrt(3));
            warnSt(1).GeneratorRatingMismatch.offenders(2:end, 3) = {pvBus.name}';
            warnSt(1).GeneratorRatingMismatch.offenders(2:end, 4) = num2cell([pvBus.kVBase]'*sqrt(3));
        end
    end
    
    %PV
    if ~strcmp(PV(1).name, 'NONE')
        PVi = 1;  % initialize our offenders counter
        Busi = 1;  % initialize our offenders' buses counter
        for ii = 1:length(PV)
            pvBusTmp = Buses(strcmp(regexprep(PV(ii).busName,'(\.[0-9]+)',''), {Buses.name})); %temp bus to save if element is bad
            
            %determine if kV needs to be scaled to l-n value
            if PV(ii).numPhases >= 2 
                LL_LN = sqrt(3); % if it is 3 phase, convert l-l value to l-n
            else
                if PV(ii).isDelta
                   LL_LN = sqrt(3); % even if it is not 3 phase but it is delta connection, convert
                   if PV(ii).deltaVoltOutOfPhase %180 out of phase, so 2 instead of 120 degrees and sqrt(3)
                        LL_LN = 2;
                   end
                else
                    LL_LN = 1; % if it is 1 phase, leave it alone
                end
            end
            
            % If our current capacitor is outside of the PV's bus's rating - save the PV, save the bus, and increment both counters
            if abs(pvBusTmp.kVBase - PV(ii).kV/LL_LN) > elmntVoltWndw*pvBusTmp.kVBase && pvBusTmp.voltagePU ~= 0
                badPV(PVi) = PV(ii);
                badPV(PVi).kV = badPV(PVi).kV/LL_LN;
                PVi = PVi + 1;
                pvBus(Busi) = pvBusTmp;
                Busi = Busi + 1;
            end
        end
     %if some offenders were found, create entry in the warning structure
        if PVi > 1
            warnSt(1).PVRatingMismatch.nm = 'PVRatingMismatch';
            warnSt(1).PVRatingMismatch.str = sprintf('%d of the PV kV ratings differ from its bus kV rating by more than %d percent.', length({badPV.name}), elmntVoltWndw*100);
            warnSt(1).PVRatingMismatch.offenders = cell(PVi, 4);
            warnSt(1).PVRatingMismatch.offenders(1, :) = [{'PV Name'} {'PV Rating (kV l-l)'} {'Bus Name'} {'Bus Rating(kV l-l)'}];
            warnSt(1).PVRatingMismatch.offenders(2:end, 1) = {badPV.name}';
            warnSt(1).PVRatingMismatch.offenders(2:end, 2) = num2cell([badPV.kV]'*sqrt(3));
            warnSt(1).PVRatingMismatch.offenders(2:end, 3) = {pvBus.name}';
            warnSt(1).PVRatingMismatch.offenders(2:end, 4) = num2cell([pvBus.kVBase]'*sqrt(3));
        end
    end

    %TRANSFORMERS
    % Xfmr Voltage Rating
    if ~strcmp(Transformers(1).name, 'NONE')
        Transi = 1; % initialize our offenders counter
        Busi = 1; % initialize our offenders' buses counter
        for ii = 1:length(Transformers)
                       
            % If our current xfmr is outside of the xfmr's bus's rating - save the xfmr, save the bus, and increment both counters
            % Repeat twice - once for each bus
            
            % Check bus 1 ratings
            
            %determine if kV needs to be scaled to l-n value
            if Transformers(ii).numPhases >= 2
                LL_LN = sqrt(3); % if it is 3 phase, convert l-l value to l-n
            else
                if Transformers(ii).isDelta1
                    LL_LN = sqrt(3); % even if it is not 3 phase but it is delta connection, convert
                    if Transformers(ii).bus1DeltaVoltOutOfPhase %180 out of phase, so 2 instead of 120 degrees and sqrt(3)
                        LL_LN = 2;
                    end
                else
                    LL_LN = 1; % if it is 1 phase, leave it alone
                end
            end
            
            transBusTmp = Buses(strcmp(regexprep(Transformers(ii).bus1,'(\.[0-9]+)',''), {Buses.name}));
            if abs(transBusTmp.kVBase - Transformers(ii).bus1kV/LL_LN) > elmntVoltWndw*transBusTmp.kVBase && transBusTmp.voltagePU ~= 0
                badTrans(Transi).name = Transformers(ii).name;
                badTrans(Transi).buskV = Transformers(ii).bus1kV/LL_LN;
                Transi = Transi + 1;
                transBus(Busi) = transBusTmp;
                Busi = Busi + 1;
            end
            
            % Check bus 2 ratings
            
            %determine if kV needs to be scaled to l-n value
            if Transformers(ii).numPhases >= 2
                LL_LN = sqrt(3); % if it is 3 phase, convert l-l value to l-n
            else
                if Transformers(ii).isDelta2
                    LL_LN = sqrt(3); % even if it is not 3 phase but it is delta connection, convert
                    if Transformers(ii).bus2DeltaVoltOutOfPhase %180 out of phase, so 2 instead of 120 degrees and sqrt(3)
                        LL_LN = 2;
                    end
                else
                    LL_LN = 1; % if it is 1 phase, leave it alone
                end
            end
            
            transBusTmp = Buses(strcmp(regexprep(Transformers(ii).bus2,'(\.[0-9]+)',''), {Buses.name}));
            if abs(transBusTmp.kVBase - Transformers(ii).bus2kV/LL_LN) > elmntVoltWndw*transBusTmp.kVBase && transBusTmp.voltagePU ~= 0
                badTrans(Transi).name = Transformers(ii).name;
                badTrans(Transi).buskV = Transformers(ii).bus2kV/LL_LN;
                Transi = Transi + 1;
                transBus(Busi) = transBusTmp;
                Busi = Busi + 1;
            end

        end
     %if some offenders were found, create entry in warning structure
        if Transi > 1
            warnSt(1).TransformerRatingMismatch.nm = 'TransformerRatingMismatch';
            warnSt(1).TransformerRatingMismatch.str = sprintf('%d of the transformer kV rating differs from its bus kV rating by more than %d percent.', length({badTrans.name}), elmntVoltWndw*100);
            warnSt(1).TransformerRatingMismatch.offenders = cell(Transi, 4);
            warnSt(1).TransformerRatingMismatch.offenders(1, :) = [{'Xfmr Name'} {'Xfmr Rating (kV l-l)'} {'Bus Name'} {'Bus Rating (kV l-l)'}];
            warnSt(1).TransformerRatingMismatch.offenders(2:end, 1) = {badTrans.name}';
            warnSt(1).TransformerRatingMismatch.offenders(2:end, 2) = num2cell([badTrans.buskV]'*sqrt(3));
            warnSt(1).TransformerRatingMismatch.offenders(2:end, 3) = {transBus.name}';
            warnSt(1).TransformerRatingMismatch.offenders(2:end, 4) = num2cell([transBus.kVBase]'*sqrt(3));
        end
    end

    %Transformers Overloaded
    if ~strcmp(Transformers(1).name, 'NONE')
        % Using a different indexing method than all the previous checks
        % Store a logical as long as the trnsformers strucuture of all transformers that are offenders
        badTransi = [Transformers.normhkva]*xfmrPowerMax < abs([Transformers.inputkva]);
        % If any of the logicals are 1, we need to store the offender - add to the warning structure
        if any(badTransi)
            badTrans = Transformers(badTransi);
            [junk,sortIndex] = sort([badTrans.inputkva]'./[badTrans.kva]'*100,'descend');
            warnSt(1).TransformerOverloaded.nm = 'TransformerOverloaded';
            warnSt(1).TransformerOverloaded.str = sprintf('%d of the transformers are overloaded, bus1 power is more than %d percent of transformer kVA rating. Check that the loads on that transformer are entered correctly.', length({badTrans.name}), xfmrPowerMax*100);
            warnSt(1).TransformerOverloaded.offenders = cell(length(badTrans)+1, 4);
            warnSt(1).TransformerOverloaded.offenders(1, :) = [{'Xfmr Name'} {'Xfmr Rating (kVA)'} {'Transformer Load (kVA)'} {'Percent Overload (%)'}];
            warnSt(1).TransformerOverloaded.offenders(2:end, 1) = {badTrans(sortIndex).name}';
            warnSt(1).TransformerOverloaded.offenders(2:end, 2) = num2cell([badTrans(sortIndex).normhkva]');
            warnSt(1).TransformerOverloaded.offenders(2:end, 3) = num2cell([badTrans(sortIndex).inputkva]');
            warnSt(1).TransformerOverloaded.offenders(2:end, 4) = num2cell([badTrans(sortIndex).inputkva]'./[badTrans(sortIndex).normhkva]'*100);
        end
    end
    
    %Transformers No Power Flow - maybe no loads on it
    if ~strcmp(Transformers(1).name, 'NONE')
        % Using a different indexing method than all the previous checks
        % Store a logical as long as the trnsformers strucuture of all transformers that are offenders
        badTransn = abs([Transformers.inputkva]) < [Transformers.kva]*xfmrPowerMin;
        % If any of the logicals are 1, we need to store the offender - add to the warning structure
        if any(badTransn)
            badTrans = Transformers(badTransn);
            loadBuses = regexprep({Loads.busName},'(\.[0-9]+)',''); % remove node refernces
            badTransTF = ones(length(badTrans), 1); % initialize tf array of noloadxfmrs to keep
            for ii = 1:length(badTrans) % remove any with existing, but really small, loads from badTrans
                if badTrans(ii).bus1kV > badTrans(ii).bus2kV
                    dnStrmBuses = findDownstreamBuses(DSSCircObj, badTrans(ii).bus2, 'Transformers', Transformers, 'Lines', Lines);
                else
                    dnStrmBuses = findDownstreamBuses(DSSCircObj, badTrans(ii).bus1, 'Transformers', Transformers, 'Lines', Lines);
                end
                [ldTF ldLoc] = ismember(dnStrmBuses,loadBuses);
                if any(ldTF)
                    badTransTF(ii) = false;
                    badTrans(ii).ldSum = sum([Loads(ldLoc(ldLoc~=0)).kW]);
                    badTrans(ii).ldNames = sprintf('%s, ', Loads(ldLoc(ldLoc~=0)).name);
                end
                
            end
            badTransLowLoad = badTrans(~badTransTF);
            badTransNoLoad = badTrans(logical(badTransTF));
            
            if ~isempty(badTransNoLoad)
                % No Load Xfmrs
                warnSt(1).TransformerNoLoad.nm = 'TransformerNoLoad';
                warnSt(1).TransformerNoLoad.str = sprintf('%d of the transformer have no load on them. Check that the loads on that transformer are entered correctly.', length({badTransNoLoad.name}));
                warnSt(1).TransformerNoLoad.offenders = cell(length(badTransNoLoad)+1, 2);
                warnSt(1).TransformerNoLoad.offenders(1, :) = [{'Xfmr Name'} {'Xfmr Rating (kVA)'}];
                warnSt(1).TransformerNoLoad.offenders(2:end, 1) = {badTransNoLoad.name}';
                warnSt(1).TransformerNoLoad.offenders(2:end, 2) = num2cell([badTransNoLoad.kva]');
            end
            
            if ~isempty(badTransLowLoad)
                % Low Load Xfmrs
                warnSt(1).TransformerLowLoad.nm = 'TransformerLowLoad';
                warnSt(1).TransformerLowLoad.str = sprintf('%d of the transformer have less than %d percent power flow of their kVA rating. Check that the loads on that transformer are entered correctly.', length({badTransLowLoad.name}), xfmrPowerMin*100);
                warnSt(1).TransformerLowLoad.offenders = cell(length(badTransLowLoad)+1, 5);
                warnSt(1).TransformerLowLoad.offenders(1, :) = [{'Xfmr Name'} {'Xfmr Rating (kVA)'} {'Transformer Load (kVA)'} {'Sum of Load kW Ratings'} {'List of Load Names'}];
                warnSt(1).TransformerLowLoad.offenders(2:end, 1) = {badTransLowLoad.name}';
                warnSt(1).TransformerLowLoad.offenders(2:end, 2) = num2cell([badTransLowLoad.kva]');
                warnSt(1).TransformerLowLoad.offenders(2:end, 3) = num2cell([badTransLowLoad.inputkva]');
                warnSt(1).TransformerLowLoad.offenders(2:end, 4) = num2cell([badTransLowLoad.ldSum])';
                warnSt(1).TransformerLowLoad.offenders(2:end, 5) = {badTransLowLoad.ldNames}';
            end
        end
    end
    
%% Bus Voltage
% BUS VOLTAGES
% Check if any bus voltages are outside the defined range
% Logical AND with a check for > 0.05 to filter out any buses that may be in a disabled section of the circuit
badBuses = Buses(abs([Buses.voltagePU] - 1.0) > busVoltWndw & [Buses.voltagePU] > 0.05);
if ~isempty(badBuses)
    warnSt(1).BusVoltage.nm = 'BusVoltage';
    warnSt(1).BusVoltage.str = sprintf('%d of the enabled bus voltages are outside of the range 1 +/- %1.2f pu. Visualize using plotVoltageProfile(DSSCircObj)', length({badBuses.name}), busVoltWndw);
    warnSt(1).BusVoltage.offenders = {badBuses.name}';
    warnSt(1).BusVoltage.offenders = cell(length(badBuses)+1, 4);
    warnSt(1).BusVoltage.offenders(1, :) = [{'Bus Name'}, {'Bus Voltage (pu)'}, {'Bus Voltage (kV)'}, {'Bus Base Voltage (kV)'}];
    warnSt(1).BusVoltage.offenders(2:end, 1) = {badBuses.name}';
    warnSt(1).BusVoltage.offenders(2:end, 2) = num2cell([badBuses.voltagePU]');
    warnSt(1).BusVoltage.offenders(2:end, 3) = num2cell([badBuses.voltage]'/1000);
    warnSt(1).BusVoltage.offenders(2:end, 4) = num2cell([badBuses.kVBase]');
end

end

%% LINES
    [tf loc] = ismember({Lines.parentObject},strcat('Line.',{Lines.name})); %only look at lines where the parent object is a line object (not transformers)
    LinesWithParents = Lines(loc>0); %if parent object is a line, store in LinesWithParents
    LineParents = Lines(loc(loc>0)); %structure of parents matching and ordered with LinesWithParents
    
    % find situations where the downstream line has more phases than the upstream (parent) line
%         issues = [LineParents.numPhases] < [LinesWithParents.numPhases]; %find locations where the parent (upstream) has fewer phases
%         if any(issues)
%             warnSt(warnCnt).nm = 'LinePhasesMismatch';
%             warnSt(warnCnt).str = sprintf('%d of the lines have more phases than possible given the number of upstream phases.', sum(issues));
%             warnSt(warnCnt).offenders = cell(length(LineParents(issues))+1, 4);
%             warnSt(warnCnt).offenders(:, 1:2) = [{'Upstream Line Names'}, {'Downstream Line Names'}; [{LineParents(issues).name}',{LinesWithParents(issues).name}']]; %table of line names with issues (upstream, downstream)
%             warnSt(warnCnt).offenders(:, 3:4) = [{'Upstream Phases'}, {'Downstream Phases'}; num2cell([[LineParents(issues).numPhases]',[LinesWithParents(issues).numPhases]'])]; %table of line phases with issues (upstream, downstream)
%             warnCnt = warnCnt + 1;
%         end
        
    % find situations where the downstream line has a significantly higher rating than the upstream (parent) line
        % Filter for greater than 0.01 length - B/C if the downstream object is a busbar, don't worry about it
            %LineParents = LineParents(abs([LinesWithParents.bus2Distance]-[LinesWithParents.bus1Distance])>=0.01); 
            %LinesWithParents = LinesWithParents(abs([LinesWithParents.bus2Distance]-[LinesWithParents.bus1Distance])>=0.01);
        %find locations where the parent (upstream) is less than the define threshold of the downstream line rating
        issues = [LineParents.lineRating].*linePercDiff < [LinesWithParents.lineRating];
        % Add to the warning structure
        if any(issues)
            warnSt(1).LineRatingMismatch.nm = 'LineRatingMismatch';
            warnSt(1).LineRatingMismatch.str = sprintf('%d of the line ratings are %d percent the size of the immediately upstream line. Visualize using plotCircuitLines(DSSCircObj,''Thickness'',''lineRating'')', size(cell(length(LineParents(issues))+1, 6),1)-1, linePercDiff*100-100);
            warnSt(1).LineRatingMismatch.offenders = cell(length(LineParents(issues))+1, 6);
            warnSt(1).LineRatingMismatch.offenders(:, 1:2) = [{'Upstream Line Names'}, {'Downstream Line Names'}; {LineParents(issues).name}',{LinesWithParents(issues).name}']; %table of line names with issues (upstream, downstream)
            warnSt(1).LineRatingMismatch.offenders(:, 3:4) = [{'Upstream Line Ratings'}, {'Downstream Line Ratings'}; num2cell([[LineParents(issues).lineRating]',[LinesWithParents(issues).lineRating]'])]; %table of line ratings with issues (upstream, downstream)
            warnSt(1).LineRatingMismatch.offenders(:, 5:6) = [{'Upstream Line Code'}, {'Downstream Line Code'}; {LineParents(issues).lineCode}',{LinesWithParents(issues).lineCode}']; %table of line codes with issues (upstream, downstream)
        end

%% Show files
%     DSSText.command = 'show eventlog';
%     DSSText.command = 'show isolated';
%     DSSText.command = 'show loops';

%% Print out warnings if the option is on
if strcmpi(Warnings, 'on')
    fNms = fieldnames(warnSt);
    for ii = 1:length(fNms)
        if strcmp(fNms{ii}, 'NoBusCoords')
            strtmp = '';
        else
            strtmp = 'View the output structure for a list of offenders.';
        end
        fprintf(1, '\n');
        warning(['GridPVToolbox:' eval(['warnSt.' fNms{ii} '.nm'])], [eval(['warnSt.' fNms{ii} '.str']) ' \n' strtmp '\n']);
    end
end

%% No errrors output
if isempty(fieldnames(warnSt)) && DSSCircObj.ActiveCircuit.Solution.Converged == 1
    warnSt(1).NoErrors.str = 'No Errors';
elseif isempty(fieldnames(warnSt)) && DSSCircObj.ActiveCircuit.Solution.Converged == 0
    warnSt = [];
    fprintf(1, '\n\n');
    warning('The circuit did not converge. CircuitCheck did not find any errors that can be checked for without a valid powerflow solution. Additional error checks require a valid powerflow solution. Rerun CircuitCheck once the circuit converges.');
end

%% Check for numerous errors
if ~isempty(warnSt) && length(fields(warnSt)) >= 3
    fprintf(1, '\n\n');
    warning('There are at least 3 warnings associated with this circuit. It is recommended you start at the top of the warning structure list when repairing.');
end

end

function isinterfaceOpenDSS_internal(DSSCircObj)

isinterface = false;

if ispc
    
    % Check if it is a COM interface
    if ~strcmp(class(DSSCircObj),'COM.OpenDSSEngine_DSS')
        error(sprintf('The specified DSSCircObj is not a COM object. \n Please refer to the Manual (DSSStartup -> Outputs) for more information.'))
    else
        % Check for OpenDSS Errors
        if DSSCircObj.Error.Number~=0 && DSSCircObj.AllowForms==0
            warning(sprintf('OpenDSS is reporting an error:\n"%s"',DSSCircObj.Error.Description))
        end
        
        % Check if it is compiled
        if DSSCircObj.NumCircuits == 0
            error(sprintf('The specified OpenDSS circuit object does not contain a compiled circuit. \n You need to run DSSText.command = ''compile YourCircuitName.dss'' \n Please refer to the Manual (DSSStartup -> Outputs) for more information.'))
        else
            % Check if solved
            if DSSCircObj.ActiveCircuit.Solution.Totaliterations == 0
                error(sprintf('The specified OpenDSS circuit object has not been solved.\nTo solve the power flow in OpenDSS, use DSSText.command = ''solve'''))
            else
                % Check if converged solution
                if DSSCircObj.ActiveCircuit.Solution.Converged == 0
                    warning(sprintf(['The specified OpenDSS circuit object contains a circuit, but the solution did not converge and is therefore invalid.\n' ...
                    'Try altering the max number of iterations OpenDSS is allowing for solving the power flow.\n Alternatively, there may be an isolated load that is being manually enabled, '...
                    'causing the power flow to not converge.']))
                else
                    isinterface = true;
                    if length(DSSCircObj.ActiveCircuit.Meters.AllNames)==1 && strcmpi(DSSCircObj.ActiveCircuit.Meters.AllNames,'NONE')
                        [msgstr, msgid] = lastwarn;
                        if ~strcmpi(msgid,'GridPV:EnergyMeter') %don't display the warning immediately again
                            warning('GridPV:EnergyMeter','The OpenDSS circuit does not contain an EnergyMeter.  GridPV Toolbox expects an EnergyMeter at the substation.  Refer to OpenDSS documentation for more details.');
                        end
                    end
                end
            end
        end   
    end

end    

end

function [Lines warnStLines warnFlag] = getLineInfo_Internal(DSSCircObj)
%The purpose of this internal function is to get line info and check for
%things in a different manner than the toolbox function, particularly - to
%check that the number of phases matches the number of decimals in the bus
%name notation (eg 'BUSNAME.1.2')

%Implement soon:
% 'One or more line was transposed and is not connected to the same phases on either end.'

    warnFlag = 0; %Indicates which errors that are set to be caught by this function have shown up
    warnStLines = struct('nm', {}, 'str', {}, 'offenders', {});
    warnStLines(1).nm = 'InvalidLineBusName';
    warnStLines(1).str = 'One or more line has a bus name that does not match the number of phases of the line. (e.g. A 2-phase line should have both bus 1 and 2 with names similar to ''BUSNAME.2.3'' with 2 phases indicated in the decimal notation.)';
    warnStLines(1).offenders = [{'Line Name'} {'NumPhases'} {'Bus1 Name'} {'Bus2 Name'}];
    warnCount = 1; %Count for offenders
    %% Define the circuit
    DSSCircuit = DSSCircObj.ActiveCircuit;

    lineName = DSSCircuit.Lines.AllNames;
    Lines = struct('name',lineName);

    %% Get all info
    for ii=1:length(Lines)
        DSSCircuit.SetActiveElement(['line.' Lines(ii).name]);

        lineBusNames = DSSCircuit.ActiveElement.BusNames;
        Lines(ii).bus1 = lineBusNames{1};
        Lines(ii).bus2 = lineBusNames{2};

        Lines(ii).enabled = DSSCircuit.ActiveElement.Enabled;
        if ~Lines(ii).enabled
            continue;
        end
        
        %Bus Distance
        % bus 1
            DSSCircuit.SetActiveBus(Lines(ii).bus1);
            if isempty(DSSCircuit.ActiveBus.Name)
                generalName = regexprep(Lines(ii).bus1,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
                DSSCircuit.SetActiveBus(generalName); 
            end
            Lines(ii).bus1Distance = DSSCircuit.ActiveBus.Distance;
            
            voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
            voltages = reshape(voltages,2,[]); %two rows for real and reactive
            voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
            Lines(ii).bus1Voltage = mean(voltages(1:DSSCircuit.ActiveBus.NumNodes));
            
        % bus 2
            DSSCircuit.SetActiveBus(Lines(ii).bus2);
            if isempty(DSSCircuit.ActiveBus.Name)
                generalName = regexprep(Lines(ii).bus2,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
                DSSCircuit.SetActiveBus(generalName); 
            end
            Lines(ii).bus2Distance = DSSCircuit.ActiveBus.Distance;
        
        %Phases
        bus1phases = length(regexp(Lines(ii).bus1,'\.([1-9]+)','Tokens'));
        if bus1phases==0 bus1phases=3; end %no bus numbering for 3-phase lines
        bus2phases = length(regexp(Lines(ii).bus2,'\.([1-9]+)','Tokens'));
        if bus2phases==0 bus2phases=3; end %no bus numbering for 3-phase lines
        numPhases = DSSCircuit.ActiveElement.NumPhases;
        if (bus1phases ~= numPhases || bus2phases ~= numPhases)
            warnStLines.offenders = [warnStLines.offenders; {Lines(ii).name} {num2str(numPhases)} {Lines(ii).bus1} {Lines(ii).bus2}];
            warnCount = warnCount+1;
            warnFlag = 1;
        end
        
        currents = DSSCircuit.ActiveCktElement.Currents; %complex currents
        currents = reshape(currents,2,[]); %two rows for real and reactive
        currents = hypot(currents(1,:),currents(2,:)); %magnitude
        Lines(ii).bus1Current = mean(currents(1:DSSCircuit.ActiveCktElement.NumPhases));
        Lines(ii).bus2Current = mean(currents(DSSCircuit.ActiveCktElement.NumConductors+1:DSSCircuit.ActiveCktElement.NumConductors+DSSCircuit.ActiveCktElement.NumPhases));

        Lines(ii).numPhases = DSSCircuit.ActiveElement.NumPhases; 
        Lines(ii).lineRating = DSSCircuit.ActiveElement.NormalAmps;
        
        % find parent (upstream) object from current active line
        DSSCircuit.ParentPDElement; %move cursor to upstream parent
        Lines(ii).parentObject = get(DSSCircuit.ActiveDSSElement,'Name');
    end

    %% Remove lines that are not enabled or have zero volts on either side
    condition = [Lines.enabled]==1;
    Lines = Lines(condition); 
    
    %% Get additional miscellaneous parameters that are present
    DSSCircuit.Lines.first;
    for ii=1:length(Lines)
        Lines(ii).lineCode = get(DSSCircuit.Lines, 'LineCode');
        Lines(ii).length = get(DSSCircuit.Lines, 'Length');
        DSSCircuit.Lines.next;
    end
    
end

function Transformers = getTransformerInfo_Internal(DSSCircObj)
%The purpose of this internal function is to get line info and check for
%things in a different manner than the toolbox function
%% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;


%% Get bus connection of each transformer
transformerName = DSSCircuit.Transformers.AllNames;
Transformers = struct('name',transformerName);

count=0;
for ii=1:length(Transformers)
    DSSCircuit.SetActiveElement(['transformer.' Transformers(ii).name]);
        
    numPhases = DSSCircuit.ActiveCktElement.NumPhases;
    
    buses = DSSCircuit.ActiveElement.BusNames;
    Transformers(ii).bus1 = buses{1};
    Transformers(ii).bus2 = buses{2};
    Transformers(ii).numPhases = DSSCircuit.ActiveElement.NumPhases;
    Transformers(ii).powers = DSSCircuit.ActiveElement.Powers;
    Transformers(ii).enabled = DSSCircuit.ActiveElement.Enabled;
    
    power = reshape(Transformers(ii).powers,2,[]); %two rows for real and reactive
    Transformers(ii).inputkva = sqrt(sum(power(1,1:DSSCircuit.ActiveCktElement.NumTerminals))^2 + sum(power(2,1:DSSCircuit.ActiveCktElement.NumTerminals))^2);
    
    DSSCircuit.SetActiveBus(Transformers(ii).bus1);
    voltages = DSSCircuit.ActiveBus.Voltages;
    if DSSCircuit.ActiveBus.NumNodes==2 && all(abs(voltages(1:2) + voltages(3:4)) < 1) %sometimes for single phase delta connected transformers, the secondary voltage phasors come out 180 out of phase instead of 120 apart
        Transformers(ii).bus1DeltaVoltOutOfPhase = 1;
    else
        Transformers(ii).bus1DeltaVoltOutOfPhase = 0;
    end
    Transformers(ii).bus1Distance = DSSCircuit.ActiveBus.Distance;
    
    DSSCircuit.SetActiveBus(Transformers(ii).bus2);
    voltages = DSSCircuit.ActiveBus.Voltages;
    if DSSCircuit.ActiveBus.NumNodes==2 && all(abs(voltages(1:2) + voltages(3:4)) < 1) %sometimes for single phase delta connected transformers, the secondary voltage phasors come out 180 out of phase instead of 120 apart
        Transformers(ii).bus2DeltaVoltOutOfPhase = 1;
    else
        Transformers(ii).bus2DeltaVoltOutOfPhase = 0;
    end
    Transformers(ii).bus2Distance = DSSCircuit.ActiveBus.Distance;
end

%% Remove xfmr that are not enabled or have zero volts on either side
    condition = [Transformers.enabled]==1;
    Transformers = Transformers(condition); 

% Get Transformer kva size
DSSCircuit.Transformers.first;
for ii=1:length(Transformers)
    DSSCircuit.Transformers.Wdg=1;
    Transformers(ii).kva = get(DSSCircuit.Transformers,'kva');
    Transformers(ii).wdg = get(DSSCircuit.Transformers,'Wdg');
    Transformers(ii).bus1kV = get(DSSCircuit.Transformers,'kV');
    Transformers(ii).isDelta1 = get(DSSCircuit.Transformers,'IsDelta');
    DSSCircuit.Transformers.Wdg=2;
    Transformers(ii).bus2kV = get(DSSCircuit.Transformers,'kV');
    Transformers(ii).isDelta2 = get(DSSCircuit.Transformers,'IsDelta');
    
    DSSText.command = sprintf('? transformer.%s.normhkva',Transformers(ii).name);
    normhkva = DSSText.result;
    if ~isempty(normhkva)
        normhkva = str2double(normhkva);
    else
        normhkva = Transformers(ii).kva*1.1;
    end
    Transformers(ii).normhkva = normhkva;
    DSSCircuit.Transformers.next;
end

end

function Buses = getBusInfo_Internal(DSSCircObj,busNames)
%% Get circuit information
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

Buses = struct('name',busNames);

for ii=1:length(busNames)
    DSSCircuit.SetActiveBus(busNames{ii});
    if isempty(DSSCircuit.ActiveBus.Name)
        generalName = regexprep(busNames{ii},'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        DSSCircuit.SetActiveBus(generalName); 
    end
    
    Buses(ii).numPhases = DSSCircuit.ActiveBus.NumNodes;
    
    voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
    busNodes = DSSCircuit.ActiveBus.Nodes;
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    Buses(ii).voltage = mean(voltages(busNodes>=1 & busNodes<=3));
    
    voltagesPU = DSSCircuit.ActiveBus.puVoltages; %complex voltages
    voltagesPU = reshape(voltagesPU,2,[]); %two rows for real and reactive
    voltagesPU = hypot(voltagesPU(1,:),voltagesPU(2,:)); %voltage magnitude
    Buses(ii).voltagePU = mean(voltagesPU(busNodes>=1 & busNodes<=3));
    
    Buses(ii).distance = DSSCircuit.ActiveBus.Distance;
    Buses(ii).kVBase = DSSCircuit.ActiveBus.kVBase;
    Buses(ii).coordinates = [DSSCircuit.ActiveBus.y, DSSCircuit.ActiveBus.x];

end

end

function Loads = getLoadInfo_Internal(DSSCircObj)

%% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

loadName = DSSCircuit.Loads.AllNames;

Loads = struct('name',loadName);

%% Get Load Info
count=0;
for ii=1:length(Loads)
    DSSCircuit.SetActiveElement(['load.' Loads(ii).name]);
    
    buses = DSSCircuit.ActiveElement.BusNames;
    Loads(ii).busName = buses{1};
    Loads(ii).numPhases = DSSCircuit.ActiveElement.NumPhases;
    Loads(ii).enabled = DSSCircuit.ActiveElement.Enabled;
    
    if ~Loads(ii).enabled
        continue
    end
    
    % Get voltages after continue to avoid isolated loads (which cause
    % errors when grabbing voltage)
    voltages = DSSCircuit.ActiveCktElement.Voltages; %complex voltages
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    
    Loads(ii).voltage = mean(voltages(1:DSSCircuit.ActiveCktElement.NumPhases));
    
    DSSCircuit.SetActiveBus(Loads(ii).busName);
    voltages = DSSCircuit.ActiveBus.Voltages;
    if DSSCircuit.ActiveBus.NumNodes==2 && all(abs(voltages(1:2) + voltages(3:4)) < 1) %sometimes for single phase delta connected transformers, the secondary voltage phasors come out 180 out of phase instead of 120 apart
        Loads(ii).deltaVoltOutOfPhase = 1;
    else
        Loads(ii).deltaVoltOutOfPhase = 0;
    end
end

%% Remove loads that are not enabled or have zero volts on either side
    condition = [Loads.enabled]==1;
    Loads = Loads(condition); 

%% Get other
DSSCircuit.Loads.first;
for ii=1:length(Loads)
    Loads(ii).xfkVA = get(DSSCircuit.Loads,'xfkVA');
    Loads(ii).kW = get(DSSCircuit.Loads,'kW');
    Loads(ii).kvar = get(DSSCircuit.Loads,'kvar');
    Loads(ii).kva = get(DSSCircuit.Loads,'kva');
    Loads(ii).kV = get(DSSCircuit.Loads,'kV');
    Loads(ii).PF = get(DSSCircuit.Loads,'PF');
    Loads(ii).isDelta = get(DSSCircuit.Loads,'IsDelta');
    DSSCircuit.Loads.next;
end
end

function Capacitors = getCapacitorInfo_Internal(DSSCircObj)
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

capacitorNames = DSSCircuit.Capacitors.AllNames;
Capacitors = struct('name',capacitorNames);

if ~strcmp(capacitorNames,'NONE')
    for ii=1:length(capacitorNames)
        DSSCircuit.SetActiveElement(['capacitor.' cell2mat(capacitorNames(ii))]);
        
        capacitorBusNames = DSSCircuit.ActiveElement.BusNames;
        Capacitors(ii).busName = capacitorBusNames{1};
        
        Capacitors(ii).numPhases = DSSCircuit.ActiveElement.NumPhases; 
        Capacitors(ii).enabled = DSSCircuit.ActiveElement.Enabled;
    end

%% Remove capacitors that are not enabled or have zero volts on either side
    condition = [Capacitors.enabled]==1;
    Capacitors = Capacitors(condition); 
    
    DSSCircuit.Capacitors.first;
    for ii=1:length(Capacitors)
        Capacitors(ii).kvar = get(DSSCircuit.Capacitors,'kvar');
        Capacitors(ii).kV = get(DSSCircuit.Capacitors,'kV');
        Capacitors(ii).isDelta = get(DSSCircuit.Capacitors,'IsDelta');
        DSSCircuit.Capacitors.next;
    end
    
end
end

function Generators = getGeneratorInfo_Internal(DSSCircObj)
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;

generatorNames = DSSCircuit.Generators.AllNames;
Generators = struct('name',generatorNames);

if ~strcmp(generatorNames,'NONE')
    
    for ii=1:length(generatorNames)
        DSSCircuit.SetActiveElement(['generator.' cell2mat(generatorNames(ii))]);
        
        pvBusNames = DSSCircuit.ActiveElement.BusNames;
        Generators(ii).busName = pvBusNames{1};
        
        Generators(ii).numPhases = DSSCircuit.ActiveElement.NumPhases;
        Generators(ii).enabled = DSSCircuit.ActiveElement.Enabled;
        
        DSSCircuit.SetActiveBus(Generators(ii).busName);
        voltages = DSSCircuit.ActiveBus.Voltages;
        if DSSCircuit.ActiveBus.NumNodes==2 && all(abs(voltages(1:2) + voltages(3:4)) < 1) %sometimes for single phase delta connected transformers, the secondary voltage phasors come out 180 out of phase instead of 120 apart
            Generators(ii).deltaVoltOutOfPhase = 1;
        else
            Generators(ii).deltaVoltOutOfPhase = 0;
        end
    end

%% Remove generators that are not enabled or have zero volts on either side
    condition = [Generators.enabled]==1;
    Generators = Generators(condition); 
    
    DSSCircuit.Generators.first;
    for ii=1:length(Generators)
        Generators(ii).kV = get(DSSCircuit.Generators,'kV');
        Generators(ii).kW = get(DSSCircuit.Generators,'kW');
        Generators(ii).kvar = get(DSSCircuit.Generators,'kvar');
        Generators(ii).PF = get(DSSCircuit.Generators,'PF');
        Generators(ii).numPhases = get(DSSCircuit.Generators,'Phases');
        DSSText.command = ['? generator.' Generators(ii).name '.conn'];
        if strcmp(DSSText.result,'delta')
            Generators(ii).isDelta = 1;
        else
            Generators(ii).isDelta = 0;
        end
        DSSCircuit.Generators.next;
    end

end
end

function PV = getPVInfo_Internal(DSSCircObj)
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;

DSSCircuit.SetActiveClass('PVSystem');
pvNames = DSSCircuit.ActiveClass.AllNames;
if isempty(pvNames)
    pvNames = 'NONE';
    PV = struct('name',pvNames);
    return;
else
    PV = struct('name',pvNames);
end

for ii=1:length(pvNames)
    DSSCircuit.SetActiveElement(PV(ii).name);
    
    pvBusNames = DSSCircuit.ActiveElement.BusNames;
    PV(ii).busName = pvBusNames{1};
    
    PV(ii).numPhases = DSSCircuit.ActiveElement.NumPhases;
    PV(ii).enabled = DSSCircuit.ActiveElement.Enabled;
    
    DSSCircuit.SetActiveBus(PV(ii).busName);
    voltages = DSSCircuit.ActiveBus.Voltages;
    if DSSCircuit.ActiveBus.NumNodes==2 && all(abs(voltages(1:2) + voltages(3:4)) < 1) %sometimes for single phase delta connected transformers, the secondary voltage phasors come out 180 out of phase instead of 120 apart
        PV(ii).deltaVoltOutOfPhase = 1;
    else
        PV(ii).deltaVoltOutOfPhase = 0;
    end
end

%% Remove PV that are not enabled or have zero volts on either side
condition = [PV.enabled]==1;
PV = PV(condition);

for ii=1:length(PV)
    DSSText.command = ['? ',PV(ii).name,'.kV'];
    PV(ii).kV = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.kva'];
    PV(ii).kVA = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.pf'];
    PV(ii).pf = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.pmpp'];
    PV(ii).pmpp = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.conn'];
    if strcmp(DSSText.result,'delta')
        PV(ii).isDelta = 1;
    else
        PV(ii).isDelta = 0;
    end
end

end

function busCoordArray = getBusCoordinatesArray_Internal(DSSCircObj)
    % Define the text interface
    DSSText = DSSCircObj.Text;

    DSSText.Command = 'Export Buscoords';
    fid = fopen(DSSText.Result);
    busCoordStruct = textscan(fid, '%s %f %f', 'delimiter', ',');
    fclose(fid);
    busCoordNames = busCoordStruct{1};
    busCoordArray = [busCoordStruct{2} busCoordStruct{3}];

    delete(DSSText.Result);
end



