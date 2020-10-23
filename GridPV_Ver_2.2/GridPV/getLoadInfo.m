%% getLoadInfo
% Gets the information for all loads in the circuit
%
%% Syntax
%  Loads = getLoadInfo(DSSCircObj)
%  Loads = getLoadInfo(DSSCircObj, loadNames)
%
%% Description
% Function to get the information about the loads in the circuit and
% return a structure with the information. If the optional input of
% loadNames is filled, the function returns information for the specified
% subset of loads, excluding the miscellaneous parameters mentioned in the
% outputs below.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|loadNames|* - optional cell array of line names to get information for
%
%% Outputs
% *|Loads|* is a structure with all the parameters for the
% loads in the active circuit.  Fields are:
%
% * _name_ - Name of the load.
% * _busName_ - Name of the associated bus.
% * _numPhases_ - Number of phases associated with the load.
% * _enabled_ - {1|0} indicates whether this element is enabled in the simulation.
% * _coordinates_ - Coordinates for the load's bus, obtained from getBusInfo.
% * _distance_ - Line distance from the load's bus to the substation, obtained from getBusInfo.
% * _current_ - average phase current
% * _phaseVoltages_ - Value of voltage magnitudes calculated from
% the complex voltage returned by OpenDSS. Length is always 3,
% returning 0 for phases not on the bus
% * _phaseVoltagesPU_ - Per-unit value of voltage magnitudes calculated from
% the complex per-unit voltage returned by OpenDSS. Length is always 3,
% returning 0 for phases not on the bus.
% * _voltage_, _voltagePU_, _voltagePhasorPU_, _phaseVoltages_, _phaseVoltagePhasors_, ... 
% _phaseVoltagePhasorsPU_, _phaseVoltagesLL_, _phaseVoltagesLLPU_, _voltageLL_, _voltageLLPU_ - voltages and voltage phasors
% * _seqVoltages_, _cplxVoltages_, _seqCurrents_, _cplxSeqCurrents_ - zero, positive, and negative sequence voltages and currents magnitude or complex phasors
% * _phasePowerReal_ - 3-element array of the real components of each
% phase's complex power injected by generator. Phases that are not present will return 0.
% * _phasePowerReactive_ - 3-element array of the imaginary components of each
% phase's complex power injected by generator. Phases that are not present will return 0.
% * _powerReal_ - Total _phasePowerReal_.
% * _powerReactive_ - Total _phasePowerReactive_.
% * _losses_ - total real and imaginary power losses
% * _phaseLosses_ - real and imaginary power losses
% * _xfkVA_ - The kVA rating of the associated transformer.
% * _kW_, _kvar_, _kva_ - Rated power of the load.
% * _kV_ - Rated voltage.
% * _PF_ - Rate power factor of the load.
% * _Idx_, _pctMean_, _pctStdDev_, _allocationFactor_, _Cfactor_, _class_, _isDelta_, _CVRcurve_, CVRwatts_, _CVRvars_, _daily_, _duty_,
% _kwhdays_, _model_, _numCust_, _Rneut_, _spectrum_, _VmaxPU_, _VminEmerg_,
% _VminNorm_, _VminPU_, _Xneut_, _yearly_, _status_, _growth_ - OpenDSS load object properties
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
% Returns load information in the circuit
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% Loads = getLoadInfo(DSSCircObj) %Get information for all loads
% Loads = getLoadInfo(DSSCircObj,DSSCircObj.ActiveCircuit.Loads.AllNames) %Get information for all loads
% Loads = getLoadInfo(DSSCircObj, {'360667000'}) %Get information for one load
% Loads = getLoadInfo(DSSCircObj, [{'530877691_1'};{'331431200'}]); %Get information for two loads
%

function Loads = getLoadInfo(DSSCircObj, varargin)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addOptional('loadNames', 'noInput', @iscellstr);

p.parse(DSSCircObj, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

try
%% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

if strcmp(loadNames, 'noInput')
    loadNames = DSSCircuit.Loads.AllNames;
end

Loads = struct('name',loadNames);

% Return if there are no loads in the circuit
if strcmp(loadNames,'NONE')
    return;
end

% Get voltages bases
kVBases = DSSCircuit.Settings.VoltageBases;
kVBases = [kVBases kVBases/sqrt(3)]; % kvBases are only LL, adding LN


%% Get Load Info    
for ii=1:length(Loads)
    DSSCircuit.SetActiveElement(['load.' Loads(ii).name]);
    if ~strcmpi(['load.' loadNames{ii}], DSSCircuit.ActiveElement.Name)
        error('loadName:notfound',sprintf('Load ''%s'' is not found in the circuit.  Check that this is a load in the compiled circuit.', loadNames{ii}))
    end
    
    buses = DSSCircuit.ActiveElement.BusNames;
    buses=buses{1};
    if ~(contains(buses,'.'))
        Loads(ii).busName = buses;
    else
        buses=buses(1:regexp(buses,'\.','once')-1);
    end
    
    Loads(ii).busName = buses;
    
    Loads(ii).numPhases = DSSCircuit.ActiveCktElement.NumPhases;
    
    Loads(ii).enabled = DSSCircuit.ActiveElement.Enabled;
    phases = [~isempty(regexp(Loads(ii).busName,'.^?(\.1)','Tokens')),~isempty(regexp(Loads(ii).busName,'.^?(\.2)','Tokens')),~isempty(regexp(Loads(ii).busName,'.^?(\.3)','Tokens'))];
    if ~any(phases ~= 0) %If all phases are blank (ie they were not put into OpenDSS), then default to all 3 phases
        phases = [true true true];
    end
    Loads(ii).phases = phases;
    if ~Loads(ii).enabled % load is not enabled, so much of the active element properties will return errors
        continue;
    end
    
    numConductors = DSSCircuit.ActiveCktElement.NumConductors;
    
    nodes = DSSCircuit.ActiveElement.nodeOrder;
    nodes1 = nodes(1:numConductors);
    nodes1 = nodes1(nodes1~=0);
    
    Loads(ii).nodes = nodes1;
    
    % Currents
    currents = DSSCircuit.ActiveCktElement.Currents; %complex currents
    currents = reshape(currents,2,[]); %two rows for real and reactive
    currents = hypot(currents(1,:),currents(2,:)); %current magnitude
    
    Loads(ii).current = mean(currents(1:DSSCircuit.ActiveCktElement.NumPhases));
    
    % Voltages    
    DSSCircuit.SetActiveBus(Loads(ii).busName);
    if isempty(DSSCircuit.ActiveBus.Name)
        generalName = regexprep(Loads(ii).busName,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        DSSCircuit.SetActiveBus(generalName);
    end
    if isempty(DSSCircuit.ActiveBus.Name)
        error('busName:notfound',sprintf('Bus ''%s'' of Load ''%s'' is not found in the circuit.  Check that this is a bus in the compiled circuit.',Transformers(ii).bus1, Transformers(ii).name))
    end
    
    Loads(ii).coordinates = [DSSCircuit.ActiveBus.y, DSSCircuit.ActiveBus.x];
    Loads(ii).distance = DSSCircuit.ActiveBus.Distance;
    Loads(ii).Yprim = DSSCircuit.ActiveCktElement.Yprim;
    
    voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
    compVoltages = voltages(1:2:end) + 1j*voltages(2:2:end);
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    Loads(ii).voltage = mean(voltages(1:DSSCircuit.ActiveBus.NumNodes));
    
    voltagesPU = DSSCircuit.ActiveBus.puVoltages; %complex voltages
    compVoltagesPU = voltagesPU(1:2:end) + 1j*voltagesPU(2:2:end);
    voltagesPU = reshape(voltagesPU,2,[]); %two rows for real and reactive
    voltagesPU = hypot(voltagesPU(1,:),voltagesPU(2,:)); %voltage magnitude
    Loads(ii).voltagePU = mean(voltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    Loads(ii).voltagePhasorPU = mean(compVoltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    
    busPhaseVoltages = zeros(1,3);
    
    phaseVoltages = zeros(1,3);
    busPhaseVoltagesPU = zeros(1,3);
    phaseVoltagesPU = zeros(1,3);
    busPhaseVoltagePhasors = zeros(1,3);
    phaseVoltagePhasors = zeros(1,3);
    busPhaseVoltagePhasorsPU = zeros(1,3);
    phaseVoltagePhasorsPU = zeros(1,3);
    
    busPhaseVoltages(DSSCircuit.ActiveBus.Nodes) = voltages;
    phaseVoltages(nodes1) = busPhaseVoltages(nodes1);
    busPhaseVoltagesPU(DSSCircuit.ActiveBus.Nodes) = voltagesPU;
    phaseVoltagesPU(nodes1) = busPhaseVoltagesPU(nodes1);
    busPhaseVoltagePhasors(DSSCircuit.ActiveBus.Nodes) = compVoltages;
    phaseVoltagePhasors(nodes1) = busPhaseVoltagePhasors(nodes1);
    busPhaseVoltagePhasorsPU(DSSCircuit.ActiveBus.Nodes) = compVoltagesPU;
    phaseVoltagePhasorsPU(nodes1) = busPhaseVoltagePhasorsPU(nodes1);
    
    Loads(ii).phaseVoltages = phaseVoltages;
    Loads(ii).phaseVoltagesPU = phaseVoltagesPU;
    Loads(ii).phaseVoltagePhasors = phaseVoltagePhasors;
    Loads(ii).phaseVoltagePhasorsPU = phaseVoltagePhasorsPU;
    
    phaseVoltagesLN = abs(phaseVoltagePhasors);
    sngPhBus = sum(phaseVoltagesLN~=0, 2) == 1;
    
    phaseVoltagesLL = phaseVoltagesLN;
    if ~sngPhBus
        phaseVoltagesLL = abs([phaseVoltagePhasors(1) - phaseVoltagePhasors(2), ...
            phaseVoltagePhasors(2) - phaseVoltagePhasors(3), phaseVoltagePhasors(3) - phaseVoltagePhasors(1)] .* ...
            [phaseVoltagesLN(1) & phaseVoltagesLN(2), phaseVoltagesLN(2) & phaseVoltagesLN(3)...
            phaseVoltagesLN(3) & phaseVoltagesLN(1)]);
    end
    
    Loads(ii).phaseVoltagesLL = phaseVoltagesLL;
    
    % get pu
    phaseVoltagesLLAvg = sum(phaseVoltagesLL)./sum(phaseVoltagesLL~=0);
    baseDiff = kVBases - phaseVoltagesLLAvg/1000;
    [~, ind] = min(abs(baseDiff), [], 2);
    phaseVoltagesLLPU = phaseVoltagesLL./kVBases(ind)' / 1000;
    Loads(ii).phaseVoltagesLLPU = phaseVoltagesLLPU;
    % avg line to line voltages
    Loads(ii).voltageLL = phaseVoltagesLLAvg;
    Loads(ii).voltageLLPU = phaseVoltagesLLAvg/kVBases(ind)' / 1000;
    
    Loads(ii).seqVoltages = DSSCircuit.ActiveCktElement.SeqVoltages;
    Loads(ii).cplxSeqVoltages = DSSCircuit.ActiveCktElement.CplxSeqVoltages;
    Loads(ii).seqCurrents = DSSCircuit.ActiveCktElement.SeqCurrents;
    Loads(ii).cplxSeqCurrents = DSSCircuit.ActiveCktElement.CplxSeqCurrents;
    
    power = DSSCircuit.ActiveCktElement.Powers; %complex
    power = reshape(power,2,[]); %two rows for real and reactive
    
    Loads(ii).phasePowerReal = power(1,1:DSSCircuit.ActiveCktElement.NumPhases);
    Loads(ii).phasePowerReactive = power(2,1:DSSCircuit.ActiveCktElement.NumPhases);
    
    Loads(ii).powerReal = sum(power(1,1:DSSCircuit.ActiveCktElement.NumPhases));
    Loads(ii).powerReactive = sum(power(2,1:DSSCircuit.ActiveCktElement.NumPhases));
    
    losses = DSSCircuit.ActiveCktElement.Losses;
    Loads(ii).losses = losses(1)/1000 + 1i*losses(2)/1000;
    
    losses = DSSCircuit.ActiveCktElement.PhaseLosses;
    losses = reshape(losses,2,[]);
    Loads(ii).phaseLosses = losses(1,:) + 1i*losses(2,:);
    
end

%% Remove loads that are not enabled if no names were input to the function
condition = [Loads.enabled]==0;
if ~isempty(varargin) && any(condition) %if the user specified the load names, return warning for that load not being enabled
    warning(sprintf('Load %s is not enabled\n',Loads(condition).name));
else
    Loads = Loads(~condition);
end

%% Get load parameters

for ii=1:length(Loads)
    DSSCircuit.Loads.name = Loads(ii).name;
    Loads(ii).xfkVA = get(DSSCircuit.Loads,'xfkVA');
    Loads(ii).kW = get(DSSCircuit.Loads,'kW');
    Loads(ii).kvar = get(DSSCircuit.Loads,'kvar');
    Loads(ii).kva = get(DSSCircuit.Loads,'kva');
    Loads(ii).kV = get(DSSCircuit.Loads,'kV');
    Loads(ii).PF = get(DSSCircuit.Loads,'PF');
    
    %Get additional miscellaneous parameters that are present but are currently unused in the toolbox
    Loads(ii).Idx = get(DSSCircuit.Loads,'Idx');
    Loads(ii).pctMean = get(DSSCircuit.Loads,'PctMean');
    Loads(ii).pctStdDev = get(DSSCircuit.Loads,'PctStdDev');
    Loads(ii).allocationFactor = get(DSSCircuit.Loads,'AllocationFactor');
    Loads(ii).Cfactor = get(DSSCircuit.Loads,'Cfactor');
    Loads(ii).class = get(DSSCircuit.Loads,'Class');
    Loads(ii).isDelta = get(DSSCircuit.Loads,'IsDelta');
    Loads(ii).CVRcurve = get(DSSCircuit.Loads,'CVRcurve');
    Loads(ii).CVRwatts = get(DSSCircuit.Loads,'CVRwatts');
    Loads(ii).CVRvars = get(DSSCircuit.Loads,'CVRvars');
    Loads(ii).daily = get(DSSCircuit.Loads,'daily');
    Loads(ii).duty = get(DSSCircuit.Loads,'duty');
    Loads(ii).kwhdays = get(DSSCircuit.Loads,'kwhdays');
    Loads(ii).model = get(DSSCircuit.Loads,'Model');
    Loads(ii).numCust = get(DSSCircuit.Loads,'NumCust');
    Loads(ii).Rneut = get(DSSCircuit.Loads,'Rneut');
    Loads(ii).spectrum = get(DSSCircuit.Loads,'Spectrum');
    Loads(ii).VmaxPU = get(DSSCircuit.Loads,'Vmaxpu');
    Loads(ii).VminEmerg = get(DSSCircuit.Loads,'Vminemerg');
    Loads(ii).VminNorm = get(DSSCircuit.Loads,'Vminnorm');
    Loads(ii).VminPU = get(DSSCircuit.Loads,'Vminpu');
    Loads(ii).Xneut = get(DSSCircuit.Loads,'Xneut');
    Loads(ii).yearly = get(DSSCircuit.Loads,'Yearly');
    Loads(ii).status = get(DSSCircuit.Loads,'Status');
    Loads(ii).growth = get(DSSCircuit.Loads,'Growth');
    
end


%% As long as you are not in faultstudy mode, remove all loads that have zero volts on either side (not disabled but are isolated from the circuit)
if ~isempty(Loads) && isempty(varargin) && ~strcmp(DSSCircuit.Solution.ModeID,'Faultstudy')
    condition = [Loads.voltage]>100;
    Loads = Loads(condition);
end

catch err
    if ~strcmp(err.identifier,'loadName:notfound')
        allLines = [err.stack.line];
        allNames = {err.stack.name};
        fprintf(1, ['\nThere was an error in ' allNames{end} ' in line ' num2str(allLines(end)) ':\n'])
        fprintf(1, ['"' err.message '"' '\n\n'])
        fprintf(1, ['About to run circuitCheck.m to ensure the circuit is set up correctly in OpenDSS.\n\n'])
        fprintf(1, 'If the problem persists, change the MATLAB debug mode by entering in the command window:\n >> dbstop if caught error\n\n')
        fprintf(1, 'Running circuitCheck.m ....................\n')

        warnSt = circuitCheck(DSSCircObj, 'Warnings', 'on')
        assignin('base','warnSt',warnSt);
    end
    rethrow(err);
end
end