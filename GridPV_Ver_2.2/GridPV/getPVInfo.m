%% getPVInfo
% Gets the information for all PV plants in the circuit
%
%% Syntax
%  PV = getPVInfo(DSSCircObj);
%  PV = getPVInfo(DSSCircObj, pvNames);
%
%% Description
% Function to get the information about the PV plants in the circuit and
% return a structure with the information. If the optional input of
% pvNames is filled, the function returns information for the specified
% subset of PV installations, excluding the miscellaneous parameters
% mentioned in the outputs below.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|pvNames|* - optional cell array of PV names to get information for
%
%% Outputs
% *|PV|* is a structure with all the parameters for the
% PV plants in the active circuit.  Fields are:
%
% * _name_ - Name of the PV source.
% * _numPhases_ - Number of phases associated with the PV.
% * _busName_ - Name of the associated bus.
% * _enabled_ - {1|0} indicates whether this element is enabled in the simulation.
% * _current_ - average phase current output
% * _coordinates_ - Coordinates for the PV bus
% * _distance_ - Line distance from the PV bus to the substation, obtained from getBusInfo.
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
% * _kW_, _kvar_, _kva_ - Rated power of the PV
% * _kV_ - Rated voltage.
% * _PF_ - Rated power factor of the PV.
% * _pmpp_ - DC power rating of the PV system.
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
% Returns PV information in the circuit
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_PV_Distributed_7_5.dss"'];
% DSSText.command = 'solve';
% PV = getPVInfo(DSSCircObj) %Get information for all PV
% PV = getPVInfo(DSSCircObj, {'PV05410_g2100nj9400'}) %Get information for one PV
% PV = getPVInfo(DSSCircObj, [{'PV05410_g2100nj9400'};{'PV05410_g2100sn7200'}]); %Get information for two PV
%

function PV = getPVInfo(DSSCircObj, varargin)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addOptional('pvNames', 'noInput', @iscellstr);

p.parse(DSSCircObj, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

try
%% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;

if strcmp(pvNames, 'noInput')
    DSSCircuit.SetActiveClass('PVSystem');
    pvNames = DSSCircuit.ActiveClass.AllNames;
    if isempty(pvNames)
        pvNames = 'NONE';
    end
end

PV = struct('name',pvNames);

% Return if there are no PV in the circuit
if strcmp(pvNames,'NONE')
    return;
end

% Get voltages bases
kVBases = DSSCircuit.Settings.VoltageBases;
kVBases = [kVBases kVBases/sqrt(3)]; % kvBases are only LL, adding LN

%% Get all info
for ii=1:length(PV)
    DSSCircuit.SetActiveElement(PV(ii).name);
    
    if ~strcmpi(['pvsystem.' pvNames{ii}], DSSCircuit.ActiveElement.Name)
        DSSCircuit.SetActiveElement(['pvsystem.', PV(ii).name]); %try with pvsystem in front of the name
        if ~strcmpi(['pvsystem.' pvNames{ii}], DSSCircuit.ActiveElement.Name)
            error('pvName:notfound',sprintf('PV ''%s'' is not found in the circuit.  Check that this is a PV system in the compiled circuit.', pvNames{ii}))
        end
    end
    
    PV(ii).numPhases = DSSCircuit.ActiveElement.numPhases;
    busNames = DSSCircuit.ActiveElement.BusNames;
    PV(ii).busName = busNames{1};
    PV(ii).enabled = DSSCircuit.ActiveElement.Enabled;
    
    if ~PV(ii).enabled % PV is not enabled, so much of the active element properties will return errors
        continue;
    end
    
    currents = DSSCircuit.ActiveCktElement.Currents; %complex voltages
    currents = reshape(currents,2,[]); %two rows for real and reactive
    currents = hypot(currents(1,:),currents(2,:)); %voltage magnitude
    
    PV(ii).current = mean(currents(1:DSSCircuit.ActiveCktElement.NumPhases));
    
    % get node info
    nodes = DSSCircuit.ActiveElement.nodeOrder;
    nodes = nodes(1:DSSCircuit.ActiveCktElement.NumConductors);
    nodes = nodes(nodes~=0);
    
    % set active bus
    generalBusName = regexprep(PV(ii).busName,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    DSSCircuit.SetActiveBus(generalBusName);
    
    if isempty(DSSCircuit.ActiveBus.Name)
        error('busName:notfound',sprintf('Bus ''%s'' of PV ''%s'' is not found in the circuit.  Check that this is a bus in the compiled circuit.',PV(ii).busName, PV(ii).name))
    end
    
    % get other information
    PV(ii).coordinates = [DSSCircuit.ActiveBus.y, DSSCircuit.ActiveBus.x];
    PV(ii).distance = DSSCircuit.ActiveBus.distance;
    buskVBase = DSSCircuit.ActiveBus.kVBase;
    
    % Begin getting voltages
    voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
    comPVoltages = voltages(1:2:end) + 1j*voltages(2:2:end);
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    PV(ii).voltage = mean(voltages(1:DSSCircuit.ActiveBus.NumNodes));
    
    voltagesPU = DSSCircuit.ActiveBus.puVoltages; %complex voltages
    comPVoltagesPU = voltagesPU(1:2:end) + 1j*voltagesPU(2:2:end);
    voltagesPU = reshape(voltagesPU,2,[]); %two rows for real and reactive
    voltagesPU = hypot(voltagesPU(1,:),voltagesPU(2,:)); %voltage magnitude
    PV(ii).voltagePU = mean(voltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    PV(ii).voltagePhasorPU = mean(comPVoltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    
    busPhaseVoltages = zeros(1,3);
    phaseVoltages = zeros(1,3);
    busPhaseVoltagesPU = zeros(1,3);
    phaseVoltagesPU = zeros(1,3);
    busPhaseVoltagePhasors = zeros(1,3);
    phaseVoltagePhasors = zeros(1,3);
    busPhaseVoltagePhasorsPU = zeros(1,3);
    phaseVoltagePhasorsPU = zeros(1,3);
    
    busPhaseVoltages(DSSCircuit.ActiveBus.Nodes) = voltages;
    phaseVoltages(nodes) = busPhaseVoltages(nodes);
    busPhaseVoltagesPU(DSSCircuit.ActiveBus.Nodes) = voltagesPU;
    phaseVoltagesPU(nodes) = busPhaseVoltagesPU(nodes);
    busPhaseVoltagePhasors(DSSCircuit.ActiveBus.Nodes) = comPVoltages;
    phaseVoltagePhasors(nodes) = busPhaseVoltagePhasors(nodes);
    busPhaseVoltagePhasorsPU(DSSCircuit.ActiveBus.Nodes) = comPVoltagesPU;
    phaseVoltagePhasorsPU(nodes) = busPhaseVoltagePhasorsPU(nodes);
    
    PV(ii).phaseVoltages = phaseVoltages;
    PV(ii).phaseVoltagesPU = phaseVoltagesPU;
    PV(ii).phaseVoltagePhasors = phaseVoltagePhasors;
    PV(ii).phaseVoltagePhasorsPU = phaseVoltagePhasorsPU;
    
    phaseVoltagesLN = abs(phaseVoltagePhasors);
    sngPhBus = sum(phaseVoltagesLN~=0, 2) == 1;
    
    phaseVoltagesLL = phaseVoltagesLN;
    if ~sngPhBus
        phaseVoltagesLL = abs([phaseVoltagePhasors(1) - phaseVoltagePhasors(2), ...
            phaseVoltagePhasors(2) - phaseVoltagePhasors(3), phaseVoltagePhasors(3) - phaseVoltagePhasors(1)] .* ...
            [phaseVoltagesLN(1) & phaseVoltagesLN(2), phaseVoltagesLN(2) & phaseVoltagesLN(3)...
            phaseVoltagesLN(3) & phaseVoltagesLN(1)]);
    end
    
    PV(ii).phaseVoltagesLL = phaseVoltagesLL;
    
    % get pu
    phaseVoltagesLLAvg = sum(phaseVoltagesLL)./sum(phaseVoltagesLL~=0);
    baseDiff = kVBases - phaseVoltagesLLAvg/1000;
    [~, ind] = min(abs(baseDiff), [], 2);
    phaseVoltagesLLPU = phaseVoltagesLL./kVBases(ind)' / 1000;
    PV(ii).phaseVoltagesLLPU = phaseVoltagesLLPU;
    
    % avg line to line voltages
    PV(ii).voltageLL = phaseVoltagesLLAvg;
    PV(ii).voltageLLPU = phaseVoltagesLLAvg/kVBases(ind)' / 1000;
    
    % re-set active element to avoid interferences of setting active bus
    DSSCircuit.SetActiveElement(PV(ii).name);
    
    power = DSSCircuit.ActiveCktElement.Powers; %complex
    power = reshape(power,2,[]); %two rows for real and reactive
    
    PV(ii).phasePowerReal = power(1,1:DSSCircuit.ActiveCktElement.NumPhases);
    PV(ii).phasePowerReactive = power(2,1:DSSCircuit.ActiveCktElement.NumPhases);
    
    PV(ii).powerReal = sum(power(1,1:DSSCircuit.ActiveCktElement.NumPhases));
    PV(ii).powerReactive = sum(power(2,1:DSSCircuit.ActiveCktElement.NumPhases));
    
    losses = DSSCircuit.ActiveCktElement.Losses;
    PV(ii).losses = losses(1)/1000 + 1i*losses(2)/1000;
    
    losses = DSSCircuit.ActiveCktElement.PhaseLosses;
    losses = reshape(losses,2,[]);
    PV(ii).phaseLosses = losses(1,:) + 1i*losses(2,:);
    

    PV(ii).seqVoltages = DSSCircuit.ActiveCktElement.SeqVoltages;
    PV(ii).cplxSeqVoltages = DSSCircuit.ActiveCktElement.CplxSeqVoltages;
    PV(ii).seqCurrents = DSSCircuit.ActiveCktElement.SeqCurrents;
    PV(ii).cplxSeqCurrents = DSSCircuit.ActiveCktElement.CplxSeqCurrents;
    PV(ii).seqPowers = DSSCircuit.ActiveCktElement.SeqPowers;
    
    DSSText.command = ['? ',PV(ii).name,'.kV'];
    PV(ii).kV = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.kva'];
    PV(ii).kVA = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.kvar'];
    PV(ii).kVAR = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.pf'];
    PV(ii).pf = str2num(DSSText.result);
    
    DSSText.command = ['? ',PV(ii).name,'.pmpp'];
    PV(ii).pmpp = str2num(DSSText.result);
end


%% Remove PV that are not enabled if no names were input to the function
condition = [PV.enabled]==0;
if ~isempty(varargin) && any(~condition) %if the user specified the PV names, return warning for that PV not being enabled
    warning(sprintf('PV %s is not enabled\n',PV(condition).name));
else
    PV = PV(~condition);
    if isempty(PV)
        PV = struct('name','NONE');
        return;
    end
end

%% Get coordinates
Buses = getBusInfo(DSSCircObj,{PV.busName},1);
busCoords = num2cell(reshape([Buses.coordinates],2,[])',2);
[PV.coordinates] = deal(busCoords{:});


%% As long as you are not in faultstudy mode, remove all PV that have zero volts on either side (not disabled but are isolated from the circuit)
if ~isempty(PV) && isempty(varargin) && ~strcmp(DSSCircuit.Solution.ModeID,'Faultstudy')
    condition = [PV.voltage]>100;
    PV = PV(condition);
end


catch err
    if ~strcmp(err.identifier,'pvName:notfound')
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