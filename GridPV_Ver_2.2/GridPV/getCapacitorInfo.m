%% getCapacitorInfo
% Gets the information for all capacitors in the circuit
%
%% Syntax
%  Capacitors = getCapacitorInfo(DSSCircObj);
%  Capacitors = getCapacitorInfo(DSSCircObj, capacitorNames);
%
%% Description
% Function to get the information about the capacitors in the circuit and
% return a structure with the information. If the optional input of
% capacitorNames is filled, the function returns information for the specified
% subset of capacitors, excluding the miscellaneous and additional 
% parameters mentioned in the outputs below.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|capacitorNames|* - optional cell array of capacitor names to get information for
%
%% Outputs
% *|Capacitors|* is a structure with all the parameters for the
% capacitors in the active circuit.  Fields are:
%
% * _name_ - The capacitor name.
% * _busName_ - Name of the associated bus.
% * _numPhases_ - Number of phases associated with the capacitor bank.
% * _enabled_ - {1|0} indicates whether this element is enabled in the simulation.
% * _coordinates_ - Coordinates for the capacitor's bus, obtained from getBusInfo.
% * _distance_ - Line distance from the capacitor's bus to the substation, obtained from getBusInfo.
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
% * _switching_ - {1|0} 1 if CapControl lists the capacitor as one of its elements, 0 otherwise.
% * _kvar_ - Total kvar, if one step, or ARRAY of kvar ratings for each step. Evenly divided among phases.
% * _isDelta_ - {1|0} 1 is it connected via delta connection, 0 otherwise.
% * _kV_ - For 2, 3-phase, kV phase-phase. Otherwise specify actual cap rating.
% * _capControl_ - Name of the CapControl element controlling the capacitor if the capacitor is being controlled.
% * _controlMode_ - Mode of control if the capacitor is being controlled.
% * _monitoredObj_, _monitoredTerm_, _CTratio_, _PTratio_, _onSetting_, _offSetting_, _Vmax_, _Vmin_,
% _useVoltOverride_, _delay_, _delayOff_, _deadTime_ - Settings from the CapControl
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
% Returns capacitor information in the circuit
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% Capacitors = getCapacitorInfo(DSSCircObj) %Get information for all capacitors
% Capacitors = getCapacitorInfo(DSSCircObj,DSSCircObj.ActiveCircuit.Capacitors.AllNames) %Get information for all capacitors
% Capacitors = getCapacitorInfo(DSSCircObj, {'cap_g2101ae7400'}) %Get information for one capacitor
% Capacitors = getCapacitorInfo(DSSCircObj, [{'cap_g2100pl6500'};{'cap_g2100fk7800'}]); %Get information for two capacitors
%

function Capacitors = getCapacitorInfo(DSSCircObj, varargin)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addOptional('capacitorNames', 'noInput', @iscellstr);

p.parse(DSSCircObj, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

try
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;

if strcmp(capacitorNames, 'noInput')
    capacitorNames = DSSCircuit.Capacitors.AllNames;
end

if length(capacitorNames) == 1 && strcmpi(capacitorNames{1},'none')
    Capacitors = [];
    return;
end

Capacitors = struct('name',capacitorNames);

% Return if there are no capacitors in the circuit
if strcmp(capacitorNames,'NONE')
    return;
end

% Get voltages bases
kVBases = DSSCircuit.Settings.VoltageBases;
kVBases = [kVBases kVBases/sqrt(3)]; % kvBases are only LL, adding LN


for ii=1:length(Capacitors)
    DSSCircuit.SetActiveElement(['capacitor.' Capacitors(ii).name]);
    
    if ~strcmpi(['capacitor.' capacitorNames{ii}], DSSCircuit.ActiveElement.Name)
        error('capName:notfound',sprintf('Capacitor ''%s'' is not found in the circuit.  Check that this is a capacitor in the compiled circuit.', capacitorNames{ii}))
    end
    
    capacitorBusNames = DSSCircuit.ActiveCktElement.BusNames;
    Capacitors(ii).busName = capacitorBusNames{1};
    phases = [~isempty(regexp(Capacitors(ii).busName,'.^?(\.1)','Tokens')),~isempty(regexp(Capacitors(ii).busName,'.^?(\.2)','Tokens')),~isempty(regexp(Capacitors(ii).busName,'.^?(\.3)','Tokens'))];
    if ~any(phases ~= 0) %If all phases are blank (ie they were not put into OpenDSS), then default to all 3 phases
        phases = [true true true];
    end
    Capacitors(ii).phases = phases;
    Capacitors(ii).numPhases = DSSCircuit.ActiveCktElement.NumPhases;
    
    Capacitors(ii).enabled = DSSCircuit.ActiveCktElement.Enabled;
    if ~Capacitors(ii).enabled % capacitor is not enabled, so much of the active element properties will return errors
        continue;
    end
    
%     Buses = getBusInfo(DSSCircObj,{Capacitors(ii).busName});
    DSSCircuit.SetActiveElement(['capacitor.' cell2mat(capacitorNames(ii))]); %reset active element after getBusInfo changes it
%     Capacitors(ii).coordinates = Buses.coordinates;
%     Capacitors(ii).distance = Buses.distance;
    
    % Currents
    Capacitors(ii).Yprim = DSSCircuit.ActiveCktElement.Yprim;
    currents = DSSCircuit.ActiveCktElement.Currents; %complex currents
    currents = reshape(currents,2,[]); %two rows for real and reactive
    currents = hypot(currents(1,:),currents(2,:)); %current magnitude
    
    Capacitors(ii).current = mean(currents(1:DSSCircuit.ActiveCktElement.NumPhases));
    
    % Voltages
    numConductors = DSSCircuit.ActiveCktElement.NumConductors;
    
    nodes = DSSCircuit.ActiveElement.nodeOrder;
    nodes1 = nodes(1:numConductors);
    nodes1 = nodes1(nodes1~=0);
    
    DSSCircuit.SetActiveBus(Capacitors(ii).busName);
    if isempty(DSSCircuit.ActiveBus.Name)
        generalName = regexprep(Capacitors(ii).busName,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        DSSCircuit.SetActiveBus(generalName);
    end
    if isempty(DSSCircuit.ActiveBus.Name)
        error('busName:notfound',sprintf('Bus ''%s'' of Capacitor ''%s'' is not found in the circuit.  Check that this is a bus in the compiled circuit.',Transformers(ii).bus1, Transformers(ii).name))
    end
    
    voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
    compVoltages = voltages(1:2:end) + 1j*voltages(2:2:end);
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    Capacitors(ii).voltage = mean(voltages(1:DSSCircuit.ActiveBus.NumNodes));
    
    voltagesPU = DSSCircuit.ActiveBus.puVoltages; %complex voltages
    compVoltagesPU = voltagesPU(1:2:end) + 1j*voltagesPU(2:2:end);
    voltagesPU = reshape(voltagesPU,2,[]); %two rows for real and reactive
    voltagesPU = hypot(voltagesPU(1,:),voltagesPU(2,:)); %voltage magnitude
    Capacitors(ii).voltagePU = mean(voltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    Capacitors(ii).voltagePhasorPU = mean(compVoltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    
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
    
    Capacitors(ii).phaseVoltages = phaseVoltages;
    Capacitors(ii).phaseVoltagesPU = phaseVoltagesPU;
    Capacitors(ii).phaseVoltagePhasors = phaseVoltagePhasors;
    Capacitors(ii).phaseVoltagePhasorsPU = phaseVoltagePhasorsPU;
    
    phaseVoltagesLN = abs(phaseVoltagePhasors);
    sngPhBus = sum(phaseVoltagesLN~=0, 2) == 1;
    
    phaseVoltagesLL = phaseVoltagesLN;
    if ~sngPhBus
        phaseVoltagesLL = abs([phaseVoltagePhasors(1) - phaseVoltagePhasors(2), ...
            phaseVoltagePhasors(2) - phaseVoltagePhasors(3), phaseVoltagePhasors(3) - phaseVoltagePhasors(1)] .* ...
            [phaseVoltagesLN(1) & phaseVoltagesLN(2), phaseVoltagesLN(2) & phaseVoltagesLN(3)...
            phaseVoltagesLN(3) & phaseVoltagesLN(1)]);
    end
    
    Capacitors(ii).phaseVoltagesLL = phaseVoltagesLL;
    
    % get pu
    phaseVoltagesLLAvg = sum(phaseVoltagesLL)./sum(phaseVoltagesLL~=0);
    baseDiff = kVBases - phaseVoltagesLLAvg/1000;
    [~, ind] = min(abs(baseDiff), [], 2);
    phaseVoltagesLLPU = phaseVoltagesLL./kVBases(ind)' / 1000;
    Capacitors(ii).phaseVoltagesLLPU = phaseVoltagesLLPU;
    % avg line to line voltages
    Capacitors(ii).voltageLL = phaseVoltagesLLAvg;
    Capacitors(ii).voltageLLPU = phaseVoltagesLLAvg/kVBases(ind)' / 1000;
    
    Capacitors(ii).seqVoltages = DSSCircuit.ActiveCktElement.SeqVoltages;
    Capacitors(ii).cplxSeqVoltages = DSSCircuit.ActiveCktElement.CplxSeqVoltages;
    Capacitors(ii).seqCurrents = DSSCircuit.ActiveCktElement.SeqCurrents;
    Capacitors(ii).cplxSeqCurrents = DSSCircuit.ActiveCktElement.CplxSeqCurrents;
    
    power = DSSCircuit.ActiveCktElement.Powers; %complex
    power = reshape(power,2,[]); %two rows for real and reactive
    
    Capacitors(ii).phasePowerReal = power(1,1:DSSCircuit.ActiveCktElement.NumPhases);
    Capacitors(ii).phasePowerReactive = power(2,1:DSSCircuit.ActiveCktElement.NumPhases);
    
    Capacitors(ii).powerReal = sum(power(1,1:DSSCircuit.ActiveCktElement.NumPhases));
    Capacitors(ii).powerReactive = sum(power(2,1:DSSCircuit.ActiveCktElement.NumPhases));
    
    losses = DSSCircuit.ActiveCktElement.Losses;
    Capacitors(ii).losses = losses(1)/1000 + 1i*losses(2)/1000;
    
    losses = DSSCircuit.ActiveCktElement.PhaseLosses;
    losses = reshape(losses,2,[]);
    Capacitors(ii).phaseLosses = losses(1,:) + 1i*losses(2,:);
    
    Capacitors(ii).switching = false;%default to not switching
    
end

%% Remove capacitor that are not enabled if no names were input to the function
condition = [Capacitors.enabled]==0;
if ~isempty(varargin) && any(condition) %if the user specified the cap names, return warning for that cap not being enabled
    warning(sprintf('Capacitor %s is not enabled\n',Capacitors(condition).name));
else
    Capacitors = Capacitors(~condition);
end

%% Get Capacitor parameters 
for ii=1:length(Capacitors)
    DSSCircuit.Capacitors.name = Capacitors(ii).name;
    Capacitors(ii).kvar = get(DSSCircuit.Capacitors,'kvar');
    Capacitors(ii).isDelta = get(DSSCircuit.Capacitors,'isDelta');
    Capacitors(ii).kV = get(DSSCircuit.Capacitors,'kV');
end

% get capacitor control information
CapControlNames = DSSCircuit.CapControls.AllNames;
if ~strcmp(CapControlNames,'NONE')
    for ii=1:length(CapControlNames)
        DSSCircuit.CapControls.Name = CapControlNames{ii};
        if DSSCircuit.CktElements(['CapControl.',CapControlNames{ii}]).Enabled
            CapControlsCapacitor = get(DSSCircuit.CapControls,'Capacitor');
            CapControlsMode = get(DSSCircuit.CapControls,'Mode');
            [tf, loc] = ismember({CapControlsCapacitor},{Capacitors.name});
            if tf
                Capacitors(loc).switching = true;
                Capacitors(loc).capControl = CapControlNames{ii};
                Capacitors(loc).controlMode = CapControlsMode;
                %Additional terms to retrieve for the control object
                Capacitors(loc).monitoredObj = get(DSSCircuit.CapControls,'MonitoredObj');
                Capacitors(loc).monitoredTerm = get(DSSCircuit.CapControls,'MonitoredTerm');
                Capacitors(loc).CTratio = get(DSSCircuit.CapControls,'CTratio');
                Capacitors(loc).PTratio = get(DSSCircuit.CapControls,'PTratio');
                Capacitors(loc).onSetting = get(DSSCircuit.CapControls,'ONSetting');
                Capacitors(loc).offSetting = get(DSSCircuit.CapControls,'OFFSetting');
                Capacitors(loc).Vmax = get(DSSCircuit.CapControls,'Vmax');
                Capacitors(loc).Vmin = get(DSSCircuit.CapControls,'Vmin');
                Capacitors(loc).useVoltOverride = get(DSSCircuit.CapControls,'UseVoltOverride');
                Capacitors(loc).delay = get(DSSCircuit.CapControls,'Delay');
                Capacitors(loc).delayOff = get(DSSCircuit.CapControls,'DelayOff');
                Capacitors(loc).deadTime = get(DSSCircuit.CapControls,'DeadTime');
            end
        end
    end
end

%% As long as you are not in faultstudy mode, remove all capacitors that have zero volts on either side (not disabled but are isolated from the circuit)
% if ~isempty(Capacitors) && isempty(varargin) && ~strcmp(DSSCircuit.Solution.ModeID,'Faultstudy')
%     condition = [Capacitors.voltage]>100;
%     Capacitors = Capacitors(condition);
% end

catch err
    
    if ~strcmp(err.identifier,'capName:notfound')
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




