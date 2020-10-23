%% getGeneratorInfo
% Gets the information for all generators in the circuit
%
%% Syntax
%  Generators = getGeneratorInfo(DSSCircObj);
%  Generators = getGeneratorInfo(DSSCircObj, generatorNames);
%
%% Description
% Function to get the information about the generators in the circuit and
% return a structure with the information. If the optional input of
% generatorNames is filled, the function returns information for the specified
% subset of generators, excluding the miscellaneous parameters
% mentioned in the outputs below.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|generatorNames|* - optional cell array of generator names to get information for
%
%% Outputs
% *|Generators|* is a structure with all the parameters for the
% generators in the active circuit.  Fields are:
%
% * _name_ - Name of the generator.
% * _busName_ - Name of the associated bus.
% * _numPhases_ - Number of phases associated with the generator.
% * _enabled_ - {1|0} indicates whether this element is enabled in the simulation.
% * _nodes_ - the connection nodes at the bus
% * _current_ - average phase current output
% * _coordinates_ - Coordinates for the bus
% * _distance_ - Line distance from the bus to the substation
% * _phaseVoltages_ - Value of voltage magnitudes calculated from
% the complex voltage returned by OpenDSS. Length is always 3,
% returning 0 for phases not on the bus
% * _phaseVoltagesPU_ - Per-unit value of voltage magnitudes calculated from
% the complex per-unit voltage returned by OpenDSS. Length is always 3,
% returning 0 for phases not on the bus.
% * _voltage_, _voltagePU_, _voltagePhasorPU_, _phaseVoltages_, _phaseVoltagePhasors_, ... 
% _phaseVoltagePhasorsPU_, _phaseVoltagesLL_, _phaseVoltagesLLPU_, _voltageLL_, _voltageLLPU_ - voltages and voltage phasors
% * _phasePowerReal_ - 3-element array of the real components of each
% phase's complex power injected by generator. Phases that are not present will return 0.
% * _phasePowerReactive_ - 3-element array of the imaginary components of each
% phase's complex power injected by generator. Phases that are not present will return 0.
% * _powerReal_ - Total _phasePowerReal_.
% * _powerReactive_ - Total _phasePowerReactive_.
% * _seqVoltages_, _cplxVoltages_, _seqCurrents_, _cplxSeqCurrents_ - zero, positive, and negative sequence voltages and currents magnitude or complex phasors
% * _losses_ - total real and imaginary power losses
% * _phaseLosses_ - real and imaginary power losses
% * _kW_, _kvar_, _kva_ - Rated power of the generator
% * _kV_ - Rated voltage.
% * _PF_ - Rate power factor of the generator.
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
% Returns generator information in the circuit
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_Generators_Distributed_7_5.dss"'];
% DSSText.command = 'solve';
% Generators = getGeneratorInfo(DSSCircObj) %Get information for all generators
% Generators = getGeneratorInfo(DSSCircObj,DSSCircObj.ActiveCircuit.Generators.AllNames) %Get information for all generators
% Generators = getGeneratorInfo(DSSCircObj, {'pvn312429_1_2_3'}) %Get information for one generator
% Generators = getGeneratorInfo(DSSCircObj, [{'pvn300557_3'};{'pvn300587_2'}]); %Get information for two generators
%

function Generators = getGeneratorInfo(DSSCircObj, varargin)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addOptional('generatorNames', 'noInput', @iscellstr);

p.parse(DSSCircObj, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

try
% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;

% Get voltages bases
kVBases = DSSCircuit.Settings.VoltageBases;
kVBases = [kVBases kVBases/sqrt(3)]; % kvBases are only LL, adding LN

if strcmp(generatorNames, 'noInput')
    generatorNames = DSSCircuit.Generators.AllNames;
end

Generators = struct('name',generatorNames);

% Return if there are no generators in the circuit
if strcmp(generatorNames,'NONE')
    return;
end

% Get voltage bases
kVBases = DSSCircuit.Settings.VoltageBases;
kVBases = [kVBases kVBases/sqrt(3)]; % kvBases are only LL, adding LN

%% Get all info
for ii=1:length(generatorNames)
    DSSCircuit.SetActiveElement(['generator.' cell2mat(generatorNames(ii))]);
    
    if ~strcmpi(['generator.' generatorNames{ii}], DSSCircuit.ActiveElement.Name)
        error('GeneratorName:notfound',sprintf('Generator ''%s'' is not found in the circuit.  Check that this is a generator in the compiled circuit.', generatorNames{ii}))
    end
    
    generatorBusNames = DSSCircuit.ActiveElement.BusNames;
    Generators(ii).busName = generatorBusNames{1};
    
    Generators(ii).numPhases = DSSCircuit.ActiveElement.NumPhases;
    
    Generators(ii).enabled = DSSCircuit.ActiveCktElement.Enabled;
    
    if ~Generators(ii).enabled % generator is not enabled, so much of the active element properties will return errors
        continue;
    end
    
    nodes = DSSCircuit.ActiveElement.nodeOrder;
    Generators(ii).nodes = nodes(nodes~=0);
    
    currents = DSSCircuit.ActiveCktElement.Currents; %complex voltages
    currents = reshape(currents,2,[]); %two rows for real and reactive
    currents = hypot(currents(1,:),currents(2,:)); %voltage magnitude
    
    Generators(ii).current = mean(currents(1:DSSCircuit.ActiveCktElement.NumPhases));
    
    % get node info
    nodes = DSSCircuit.ActiveElement.nodeOrder;
    nodes = nodes(1:DSSCircuit.ActiveCktElement.NumConductors);
    nodes = nodes(nodes~=0);
    
    % set active bus
    generalBusName = regexprep(Generators(ii).busName,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
    DSSCircuit.SetActiveBus(generalBusName);
    
    if isempty(DSSCircuit.ActiveBus.Name)
        error('busName:notfound',sprintf('Bus ''%s'' of generator ''%s'' is not found in the circuit.  Check that this is a bus in the compiled circuit.',Generators(ii).busName, Generators(ii).name))
    end
    
    % get other information
    Generators(ii).coordinates = [DSSCircuit.ActiveBus.y, DSSCircuit.ActiveBus.x];
    Generators(ii).distance = DSSCircuit.ActiveBus.distance;
    buskVBase = DSSCircuit.ActiveBus.kVBase;
    
    % Begin getting voltages
    voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
    compVoltages = voltages(1:2:end) + 1j*voltages(2:2:end);
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    Generators(ii).voltage = mean(voltages(1:DSSCircuit.ActiveBus.NumNodes));
    
    voltagesPU = DSSCircuit.ActiveBus.puVoltages; %complex voltages
    compVoltagesPU = voltagesPU(1:2:end) + 1j*voltagesPU(2:2:end);
    voltagesPU = reshape(voltagesPU,2,[]); %two rows for real and reactive
    voltagesPU = hypot(voltagesPU(1,:),voltagesPU(2,:)); %voltage magnitude
    Generators(ii).voltagePU = mean(voltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    Generators(ii).voltagePhasorPU = mean(compVoltagesPU(1:DSSCircuit.ActiveBus.NumNodes));
    
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
    busPhaseVoltagePhasors(DSSCircuit.ActiveBus.Nodes) = compVoltages;
    phaseVoltagePhasors(nodes) = busPhaseVoltagePhasors(nodes);
    busPhaseVoltagePhasorsPU(DSSCircuit.ActiveBus.Nodes) = compVoltagesPU;
    phaseVoltagePhasorsPU(nodes) = busPhaseVoltagePhasorsPU(nodes);
    
    Generators(ii).phaseVoltages = phaseVoltages;
    Generators(ii).phaseVoltagesPU = phaseVoltagesPU;
    Generators(ii).phaseVoltagePhasors = phaseVoltagePhasors;
    Generators(ii).phaseVoltagePhasorsPU = phaseVoltagePhasorsPU;
    
    phaseVoltagesLN = abs(phaseVoltagePhasors);
    sngPhBus = sum(phaseVoltagesLN~=0, 2) == 1;
    
    phaseVoltagesLL = phaseVoltagesLN;
    if ~sngPhBus
        phaseVoltagesLL = abs([phaseVoltagePhasors(1) - phaseVoltagePhasors(2), ...
            phaseVoltagePhasors(2) - phaseVoltagePhasors(3), phaseVoltagePhasors(3) - phaseVoltagePhasors(1)] .* ...
            [phaseVoltagesLN(1) & phaseVoltagesLN(2), phaseVoltagesLN(2) & phaseVoltagesLN(3)...
            phaseVoltagesLN(3) & phaseVoltagesLN(1)]);
    end
    
    Generators(ii).phaseVoltagesLL = phaseVoltagesLL;
    
    % get pu
    phaseVoltagesLLAvg = sum(phaseVoltagesLL)./sum(phaseVoltagesLL~=0);
    baseDiff = kVBases - phaseVoltagesLLAvg/1000;
    [~, ind] = min(abs(baseDiff), [], 2);
    phaseVoltagesLLPU = phaseVoltagesLL./kVBases(ind)' / 1000;
    Generators(ii).phaseVoltagesLLPU = phaseVoltagesLLPU;
    
    % avg line to line voltages
    Generators(ii).voltageLL = phaseVoltagesLLAvg;
    Generators(ii).voltageLLPU = phaseVoltagesLLAvg/kVBases(ind)' / 1000;
    
    % re-set active element to avoid interferences of setting active bus
    DSSCircuit.SetActiveElement(['generator.' cell2mat(generatorNames(ii))]);
    
    power = DSSCircuit.ActiveCktElement.Powers; %complex
    power = reshape(power,2,[]); %two rows for real and reactive
    
    Generators(ii).phasePowerReal = power(1,1:DSSCircuit.ActiveCktElement.NumPhases);
    Generators(ii).phasePowerReactive = power(2,1:DSSCircuit.ActiveCktElement.NumPhases);
    
    Generators(ii).powerReal = sum(power(1,1:DSSCircuit.ActiveCktElement.NumPhases));
    Generators(ii).powerReactive = sum(power(2,1:DSSCircuit.ActiveCktElement.NumPhases));
    
    %Get additional miscellaneous parameters that are present but are currently unused in the toolbox
    Generators(ii).numTerminals = DSSCircuit.ActiveCktElement.NumTerminals;
    Generators(ii).numPhases = DSSCircuit.ActiveCktElement.NumPhases;
    
    losses = DSSCircuit.ActiveCktElement.Losses;
    Generators(ii).losses = losses(1)/1000 + 1i*losses(2)/1000;
    
    losses = DSSCircuit.ActiveCktElement.PhaseLosses;
    losses = reshape(losses,2,[]);
    Generators(ii).phaseLosses = losses(1,:) + 1i*losses(2,:);
    
    Generators(ii).seqVoltages = DSSCircuit.ActiveCktElement.SeqVoltages;
    Generators(ii).cplxSeqVoltages = DSSCircuit.ActiveCktElement.CplxSeqVoltages;
    Generators(ii).seqCurrents = DSSCircuit.ActiveCktElement.SeqCurrents;
    Generators(ii).cplxSeqCurrents = DSSCircuit.ActiveCktElement.CplxSeqCurrents;
    Generators(ii).seqPowers = DSSCircuit.ActiveCktElement.SeqPowers;
end



%% Remove generators that are not enabled if no names were input to the function
condition = [Generators.enabled]==0;
if ~isempty(varargin) && any(condition) %if the user specified the generator names, return warning for that generator not being enabled
    warning(sprintf('Generator %s is not enabled\n',Generators(condition).name));
else
    Generators = Generators(~condition);
    if isempty(Generators)
        Generators = struct('name','NONE');
        return;
    end
end


%% Get generator parameters

Buses = getBusInfo(DSSCircObj,{Generators.busName},1);
busCoords = num2cell(reshape([Buses.coordinates],2,[])',2);
[Generators.coordinates] = deal(busCoords{:});

for ii=1:length(Generators)
    DSSCircuit.Generators.name = Generators(ii).name;
    Generators(ii).kV = get(DSSCircuit.Generators,'kV');
    DSSText.command = sprintf('? Generator.%s.kw',Generators(ii).name);
    Generators(ii).kW = str2double(DSSText.result);
    Generators(ii).kvar = get(DSSCircuit.Generators,'kvar');
    Generators(ii).PF = get(DSSCircuit.Generators,'PF');
    Generators(ii).numPhases = get(DSSCircuit.Generators,'Phases');
end


%% As long as you are not in faultstudy mode, remove all generators that have zero volts on either side (not disabled but are isolated from the circuit)
if ~isempty(Generators) && isempty(varargin) && ~strcmp(DSSCircuit.Solution.ModeID,'Faultstudy')
    condition = [Generators.voltage]>100;
    Generators = Generators(condition);
end

catch err
    if ~strcmp(err.identifier,'GeneratorName:notfound')
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



