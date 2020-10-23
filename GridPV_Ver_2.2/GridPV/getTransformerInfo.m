%% getTransformerInfo
% Gets the information for all transformers in the circuit
%
%% Syntax
%  Transformers = getTransformerInfo(DSSCircObj)
%  Transformers = getTransformerInfo(DSSCircObj, transformerNames)
%
%% Description
% Function to get the information about the transformers in the circuit and
% return a structure with the information. If the optional input of
% transformerNames is filled, the function returns information for the specified
% subset of transformers, excluding the miscellaneous parameters mentioned
% in the outputs below.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|transformerNames|* - optional cell array of transformer names to get information for
%
%% Outputs
% *|Transformers|* is a structure with all the parameters for the
% transformers in the active circuit.  Fields are:
%
% * _name_ - Name of the transformer.
% * _bus1_ - Primary bus.
% * _bus2_ - Secondary bus.
% * _numPhases_ - Number of phases associated with the transformer.
% * _enabled_ - {1|0} indicates whether this element is enabled in the simulation.
% * _bus1PhasePowerReal_ - 3-element array of the real components of each
% phase's complex power at bus 1. Phases that are not present will return 0.
% * _bus1PhasePowerReactive_ - 3-element array of the imaginary components of each
% phase's complex power at bus 1. Phases that are not present will return 0.
% * _bus2PhasePowerReal_ - 3-element array of the real components of each
% phase's complex power at bus 2. Phases that are not present will return 0.
% * _bus2PhasePowerReactive_ - 3-element array of the imaginary components of each
% phase's complex power at bus 2. Phases that are not present will return 0.
% * _bus1PowerReal_ - Total real component at bus 1 of all present phases.
% * _bus1PowerReactive_ - Total imaginary component at bus 1 of all present phases.
% * _bus2PowerReal_ - Total real component at bus 2 of all present phases.
% * _bus2PowerReactive_ - Total imaginary component at bus 2 of all present phases.
% * _bus1NodeOrder_, _bus1Coordinates_, _bus1Distance_, _bus1PhaseVoltages_,
%  _bus1PhaseVoltagesPU_, _bus1Voltage_, _bus1VoltagePU_, _bus1VoltagePhasors_,
%  _bus1PhaseVoltagesLL_, _bus1PhaseVoltagesLLPU_, - Information
% regarding the starting bus. All obtained from the corresponding fields of
% the structure returned by getBusInfo when called with 'bus1' as an
% input.
% * _bus2NodeOrder_, _bus2Coordinates_, _bus2Distance_, _bus2PhaseVoltages_,
%  _bus2PhaseVoltagesPU_, _bus2Voltage_, _bus2VoltagePU_, _bus2VoltagePhasors_, 
% _bus2PhaseVoltagePhasorsPU_, _bus2PhaseVoltagePhasorsPU_, _bus2PhaseVoltagesLL_,
% _bus2PhaseVoltagesLLPU_ - Information
% regarding the ending bus. All obtained from the corresponding fields of
% the structure returned by getBusInfo when called with 'bus2' as an
% input. 
% * _losses_ - total real and imaginary power losses
% * _phaseLosses_ - real and imaginary power losses
% * _seqVoltages_, _cplxVoltages_, _seqCurrents_, _cplxSeqCurrents_ - zero, positive, and negative sequence voltages and currents magnitude or complex phasors
% * _inputkva_ - apparent power magnitude coming into the transformer
% * _controlled_ - Whether or not the transformer is tap-controlled.
% * _kva_ - Transformer power rating.
% * _XfmrCode_ - name of the transformer code if one is assigned
% * _wdg1R_, _wdg1Tap_, _wdg1minTap_, _wdg1maxTap_, _wdg1numTaps_, _bus1kV_, _isDelta_ - properties for the first winding of the transformer.  All properties exist for the wdg2 side as well
% * _Xneut_, _Rneut_, _Xhl_, _Xht_, _PCTnoLoadLoss_ - OpenDSS Transformer properties
% * _CTPrimary_, _delay_, _forwardBand_, _forwardR_, _forwardVreg_, _forwardX_, _isInverseTime_, __isReversible_, _maxTapChange_, _monitoredBus_,
% _PTratio_, _reverseBand_, _reverseR_, _reverseVreg_, _reverseX_, _tapDelay_, _tapWinding_, _voltageLimit_, _winding_ - OpenDSS properties for the RegControl object for transformers on which one is present
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
% Returns transformer information in the circuit
%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% Transformers = getTransformerInfo(DSSCircObj) %Get information for all transformers
% Transformers = getTransformerInfo(DSSCircObj,DSSCircObj.ActiveCircuit.Transformers.AllNames) %Get information for all transformers
% Transformers = getTransformerInfo(DSSCircObj, {'05410_g2101ak7700'}) %Get information for one transformer
% Transformers = getTransformerInfo(DSSCircObj, [{'05410_g2101ah4300'};{'05410_g2101ae2300'}]); %Get information for two transformers
%

function Transformers = getTransformerInfo(DSSCircObj, varargin)
%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);
p.addOptional('transformerNames', 'noInput', @iscellstr);

p.parse(DSSCircObj, varargin{:}); %parse inputs

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

try    
%% Define the circuit
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;

% Get voltages bases
kVBases = DSSCircuit.Settings.VoltageBases;
kVBases = [kVBases kVBases/sqrt(3)]; % kvBases are only LL, adding LN

%% Get bus connection of each transformer

if strcmp(transformerNames, 'noInput')
    transformerNames = DSSCircuit.Transformers.AllNames;
end

if length(transformerNames) == 1 && strcmpi(transformerNames{1},'none')
    Transformers = [];
    return;
end

Transformers = struct('name',transformerNames);

% Return if there are no transformers in the circuit
if strcmp(transformerNames,'NONE')
    return;
end

% Get voltage bases
kVBases = DSSCircuit.Settings.VoltageBases;
kVBases = [kVBases kVBases/sqrt(3)]; % kvBases are only LL, adding LN

%% Get all info
for ii=1:length(Transformers)
    DSSCircuit.SetActiveElement(['transformer.' Transformers(ii).name]);
    if ~strcmpi(['transformer.' transformerNames{ii}], DSSCircuit.ActiveElement.Name)
        error('xfmrName:notfound',sprintf('Transformer ''%s'' is not found in the circuit.  Check that this is a transformer in the compiled circuit.', transformerNames{ii}))
    end
    
    buses = DSSCircuit.ActiveElement.BusNames;
    Transformers(ii).bus1 = buses{1};
    Transformers(ii).bus2 = buses{2};
    numPhases = DSSCircuit.ActiveElement.NumPhases;
    Transformers(ii).numPhases = numPhases;
    Transformers(ii).enabled = DSSCircuit.ActiveElement.Enabled;
    
    phases = [~isempty(regexp(Transformers(ii).bus1,'.^?(\.1)','Tokens')),~isempty(regexp(Transformers(ii).bus1,'.^?(\.2)','Tokens')),~isempty(regexp(Transformers(ii).bus1,'.^?(\.3)','Tokens'))];
    if ~any(phases ~= 0) %If all phases are blank (ie they were not put into OpenDSS), then default to all 3 phases
        phases = [true true true];
    end
    Transformers(ii).phases = phases;
    
    if ~Transformers(ii).enabled % transformer is not enabled, so much of the active element properties will return errors
        continue;
    end
    
    numConductors = DSSCircuit.ActiveCktElement.NumConductors;
    
    nodes = DSSCircuit.ActiveElement.nodeOrder;
    nodes1 = nodes(1:numConductors);
    nodes1 = nodes1(nodes1>=1 & nodes1<=3);
    nodes2 = nodes(numConductors+1:numConductors*2);
    nodes2 = nodes2(nodes2>=1 & nodes2<=3);
    
    power = DSSCircuit.ActiveCktElement.Powers; %complex
    power = reshape(power,2,[]); %two rows for real and reactive
    
    bus1PhasePowerReal = zeros(1,3);
    bus1PhasePowerReactive = zeros(1,3);
    bus2PhasePowerReal = zeros(1,3);
    bus2PhasePowerReactive = zeros(1,3);
    bus1PhaseCurrent = zeros(1,3);
    bus2PhaseCurrent = zeros(1,3);
    
    bus1PhasePowerReal(nodes1) = power(1,1:length(nodes1));
    bus1PhasePowerReactive(nodes1) = power(2,1:length(nodes1));
    bus2PhasePowerReal(nodes2) = power(1,DSSCircuit.ActiveCktElement.NumConductors+1:DSSCircuit.ActiveCktElement.NumConductors+length(nodes2));
    bus2PhasePowerReactive(nodes2) = power(2,DSSCircuit.ActiveCktElement.NumConductors+1:DSSCircuit.ActiveCktElement.NumConductors+length(nodes2));
    
    Transformers(ii).bus1PhasePowerReal = bus1PhasePowerReal;
    Transformers(ii).bus1PhasePowerReactive = bus1PhasePowerReactive;
    Transformers(ii).bus2PhasePowerReal = bus2PhasePowerReal;
    Transformers(ii).bus2PhasePowerReactive = bus2PhasePowerReactive;
    
    Transformers(ii).bus1PowerReal = sum(power(1,1:length(nodes1)));
    Transformers(ii).bus1PowerReactive = sum(power(2,1:length(nodes1)));
    Transformers(ii).bus2PowerReal = sum(power(1,DSSCircuit.ActiveCktElement.NumConductors+1:DSSCircuit.ActiveCktElement.NumConductors+length(nodes2)));
    Transformers(ii).bus2PowerReactive = sum(power(2,DSSCircuit.ActiveCktElement.NumConductors+1:DSSCircuit.ActiveCktElement.NumConductors+length(nodes2)));
    
    
    % Get voltages
    numConductors = DSSCircuit.ActiveCktElement.NumConductors;
    
    % Voltage Phasors
    nodes = DSSCircuit.ActiveElement.nodeOrder;
    nodes1 = nodes(1:numConductors);
    nodes1 = nodes1(nodes1>=1 & nodes1<=3);
    nodes2 = nodes(numConductors+1:numConductors*2);
    nodes2 = nodes2(nodes2>=1 & nodes2<=3);
    Transformers(ii).bus1NodeOrder = nodes1;
    Transformers(ii).bus2NodeOrder = nodes2;
    
    % Voltages
    % bus1
    DSSCircuit.SetActiveBus(Transformers(ii).bus1);
    if isempty(DSSCircuit.ActiveBus.Name)
        generalName = regexprep(Transformers(ii).bus1,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        DSSCircuit.SetActiveBus(generalName);
    end
    if isempty(DSSCircuit.ActiveBus.Name)
        error('busName:notfound',sprintf('Bus ''%s'' of Transformer ''%s'' is not found in the circuit.  Check that this is a bus in the compiled circuit.',Transformers(ii).bus1, Transformers(ii).name))
    end
    
    voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
    busNodes = DSSCircuit.ActiveBus.Nodes;
    compVoltages = voltages(1:2:end) + 1j*voltages(2:2:end);
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    Transformers(ii).bus1Voltage = mean(voltages(busNodes>=1 & busNodes<=3));
    
    voltagesPU = DSSCircuit.ActiveBus.puVoltages; %complex voltages
    compVoltagesPU = voltagesPU(1:2:end) + 1j*voltagesPU(2:2:end);
    voltagesPU = reshape(voltagesPU,2,[]); %two rows for real and reactive
    voltagesPU = hypot(voltagesPU(1,:),voltagesPU(2,:)); %voltage magnitude
    Transformers(ii).bus1VoltagePU = mean(voltagesPU(busNodes>=1 & busNodes<=3));
    
    busPhaseVoltages = zeros(1,3);
    
    phaseVoltages = zeros(1,3);
    busPhaseVoltagesPU = zeros(1,3);
    phaseVoltagesPU = zeros(1,3);
    busPhaseVoltagePhasors = zeros(1,3);
    phaseVoltagePhasors = zeros(1,3);
    busPhaseVoltagePhasorsPU = zeros(1,3);
    phaseVoltagePhasorsPU = zeros(1,3);
    
    busPhaseVoltages(busNodes) = voltages;
    phaseVoltages(nodes1) = busPhaseVoltages(nodes1);
    busPhaseVoltagesPU(busNodes) = voltagesPU;
    phaseVoltagesPU(nodes1) = busPhaseVoltagesPU(nodes1);
    busPhaseVoltagePhasors(busNodes) = compVoltages;
    phaseVoltagePhasors(nodes1) = busPhaseVoltagePhasors(nodes1);
    busPhaseVoltagePhasorsPU(busNodes) = compVoltagesPU;
    phaseVoltagePhasorsPU(nodes1) = busPhaseVoltagePhasorsPU(nodes1);
    
    Transformers(ii).bus1PhaseVoltages = phaseVoltages;
    Transformers(ii).bus1PhaseVoltagesPU = phaseVoltagesPU;
    Transformers(ii).bus1PhaseVoltagePhasors = phaseVoltagePhasors;
    Transformers(ii).bus1PhaseVoltagePhasorsPU = phaseVoltagePhasorsPU;
    
    phaseVoltagesLN = abs(phaseVoltagePhasors);
    sngPhBus = sum(phaseVoltagesLN~=0, 2) == 1;
    
    phaseVoltagesLL = phaseVoltagesLN;
    if ~sngPhBus
        phaseVoltagesLL = abs([phaseVoltagePhasors(1) - phaseVoltagePhasors(2), ...
            phaseVoltagePhasors(2) - phaseVoltagePhasors(3), phaseVoltagePhasors(3) - phaseVoltagePhasors(1)] .* ...
            [phaseVoltagesLN(1) & phaseVoltagesLN(2), phaseVoltagesLN(2) & phaseVoltagesLN(3)...
            phaseVoltagesLN(3) & phaseVoltagesLN(1)]);
    end
    
    Transformers(ii).bus1PhaseVoltagesLL = phaseVoltagesLL;
    
    % get pu
    phaseVoltagesLLAvg = sum(phaseVoltagesLL)./sum(phaseVoltagesLL~=0);
    baseDiff = kVBases - phaseVoltagesLLAvg/1000;
    [~, ind] = min(abs(baseDiff), [], 2);
    phaseVoltagesLLPU = phaseVoltagesLL./kVBases(ind)' / 1000;
    Transformers(ii).bus1PhaseVoltagesLLPU = phaseVoltagesLLPU;
    % avg line to line voltages
    Transformers(ii).bus1VoltageLL = phaseVoltagesLLAvg;
    Transformers(ii).bus1VoltageLLPU = phaseVoltagesLLAvg/kVBases(ind)' / 1000;
    Transformers(ii).bus1kVBase = DSSCircuit.ActiveBus.kVBase;
    
    % bus2
    DSSCircuit.SetActiveBus(Transformers(ii).bus2);
    if isempty(DSSCircuit.ActiveBus.Name)
        generalName = regexprep(Transformers(ii).bus2,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        DSSCircuit.SetActiveBus(generalName);
    end
    if isempty(DSSCircuit.ActiveBus.Name)
        error('busName:notfound',sprintf('Bus ''%s'' of Transformer ''%s'' is not found in the circuit.  Check that this is a bus in the compiled circuit.',Transformers(ii).bus1, Transformers(ii).name))
    end
    voltages = DSSCircuit.ActiveBus.Voltages; %complex voltages
    busNodes = DSSCircuit.ActiveBus.Nodes;
    compVoltages = voltages(1:2:end) + 1j*voltages(2:2:end);
    voltages = reshape(voltages,2,[]); %two rows for real and reactive
    voltages = hypot(voltages(1,:),voltages(2,:)); %voltage magnitude
    Transformers(ii).bus2Voltage = mean(voltages(busNodes>=1 & busNodes<=3));
    
    voltagesPU = DSSCircuit.ActiveBus.puVoltages; %complex voltages
    compVoltagesPU = voltagesPU(1:2:end) + 1j*voltagesPU(2:2:end);
    voltagesPU = reshape(voltagesPU,2,[]); %two rows for real and reactive
    voltagesPU = hypot(voltagesPU(1,:),voltagesPU(2,:)); %voltage magnitude
    Transformers(ii).bus2VoltagePU = mean(voltagesPU(busNodes>=1 & busNodes<=3));
    
    busPhaseVoltages = zeros(1,3);
    phaseVoltages = zeros(1,3);
    busPhaseVoltagesPU = zeros(1,3);
    phaseVoltagesPU = zeros(1,3);
    busPhaseVoltagePhasors = zeros(1,3);
    phaseVoltagePhasors = zeros(1,3);
    busPhaseVoltagePhasorsPU = zeros(1,3);
    phaseVoltagePhasorsPU = zeros(1,3);
    
    busPhaseVoltages(busNodes) = voltages;
    phaseVoltages(nodes2) = busPhaseVoltages(nodes2);
    busPhaseVoltagesPU(busNodes) = voltagesPU;
    phaseVoltagesPU(nodes2) = busPhaseVoltagesPU(nodes2);
    busPhaseVoltagePhasors(busNodes) = compVoltages;
    phaseVoltagePhasors(nodes2) = busPhaseVoltagePhasors(nodes2);
    busPhaseVoltagePhasorsPU(busNodes) = compVoltagesPU;
    phaseVoltagePhasorsPU(nodes2) = busPhaseVoltagePhasorsPU(nodes2);
    
    Transformers(ii).bus2PhaseVoltages = phaseVoltages;
    Transformers(ii).bus2PhaseVoltagesPU = phaseVoltagesPU;
    Transformers(ii).bus2PhaseVoltagePhasors = phaseVoltagePhasors;
    Transformers(ii).bus2PhaseVoltagePhasorsPU = phaseVoltagePhasorsPU;
    
    phaseVoltagesLN = abs(phaseVoltagePhasors);
    sngPhBus = sum(phaseVoltagesLN~=0, 2) == 1;
    
    phaseVoltagesLL = phaseVoltagesLN;
    if ~sngPhBus
        phaseVoltagesLL = abs([phaseVoltagePhasors(1) - phaseVoltagePhasors(2), ...
            phaseVoltagePhasors(2) - phaseVoltagePhasors(3), phaseVoltagePhasors(3) - phaseVoltagePhasors(1)] .* ...
            [phaseVoltagesLN(1) & phaseVoltagesLN(2), phaseVoltagesLN(2) & phaseVoltagesLN(3)...
            phaseVoltagesLN(3) & phaseVoltagesLN(1)]);
    end
    
    Transformers(ii).bus2PhaseVoltagesLL = phaseVoltagesLL;
    
    % get pu
    phaseVoltagesLLAvg = sum(phaseVoltagesLL)./sum(phaseVoltagesLL~=0);
    baseDiff = kVBases - phaseVoltagesLLAvg/1000;
    [~, ind] = min(abs(baseDiff), [], 2);
    phaseVoltagesLLPU = phaseVoltagesLL./kVBases(ind)' / 1000;
    Transformers(ii).bus2PhaseVoltagesLLPU = phaseVoltagesLLPU;
    % avg line to line voltages
    Transformers(ii).bus2VoltageLL = phaseVoltagesLLAvg;
    Transformers(ii).bus2VoltageLLPU = phaseVoltagesLLAvg/kVBases(ind)' / 1000;
    Transformers(ii).bus2kVBase = DSSCircuit.ActiveBus.kVBase;
    
    % other fields
    Transformers(ii).bus1Distance = DSSCircuit.ActiveBus.Distance;
    DSSCircuit.SetActiveBus(Transformers(ii).bus2);
    if isempty(DSSCircuit.ActiveBus.Name)
        generalName = regexprep(Transformers(ii).bus2,'(\.[0-9]+)',''); %take out the phase numbers on buses if they have them
        DSSCircuit.SetActiveBus(generalName);
    end
    if isempty(DSSCircuit.ActiveBus.Name)
        error('busName:notfound',sprintf('Bus ''%s'' of Transformer ''%s'' is not found in the circuit.  Check that this is a bus in the compiled circuit.',Transformers(ii).bus2, Transformers(ii).name))
    end
    Transformers(ii).bus2Distance = DSSCircuit.ActiveBus.Distance;
    
    losses = DSSCircuit.ActiveCktElement.Losses;
    Transformers(ii).losses = losses(1)/1000 + 1i*losses(2)/1000;
    Transformers(ii).Yprim = DSSCircuit.ActiveCktElement.Yprim;
    losses = DSSCircuit.ActiveCktElement.PhaseLosses;
    losses = reshape(losses,2,[]);
    Transformers(ii).phaseLosses = losses(1,:) + 1i*losses(2,:);
    

    Transformers(ii).seqVoltages = DSSCircuit.ActiveElement.SeqVoltages;
    Transformers(ii).cplxSeqVoltages = DSSCircuit.ActiveElement.CplxSeqVoltages;
    Transformers(ii).seqCurrents = DSSCircuit.ActiveElement.SeqCurrents;
    Transformers(ii).cplxSeqCurrents = DSSCircuit.ActiveElement.CplxSeqCurrents;
    Transformers(ii).seqPower = DSSCircuit.ActiveElement.SeqPowers;
    
    power = reshape(DSSCircuit.ActiveElement.Powers,2,[]); %two rows for real and reactive
    Transformers(ii).inputkva = sqrt(sum(power(1,1:DSSCircuit.ActiveCktElement.NumTerminals))^2 + sum(power(2,1:DSSCircuit.ActiveCktElement.NumTerminals))^2);

    Transformers(ii).controlled = false; % Default is 0
    
end


%% Remove transformers that are not enabled if no names were input to the function
condition = [Transformers.enabled]==0;
if ~isempty(varargin) && any(condition) %if the user specified the load names, return warning for that load not being enabled
    warning(sprintf('Transformer %s is not enabled\n',Transformers(condition).name));
else
    Transformers = Transformers(~condition);
end


%% Get transformer parameters
for ii=1:length(Transformers)
    normhkva = DSSCircuit.CktElements(['transformer.',Transformers(ii).name]).Properties('normhkva').val;
    DSSCircuit.Transformers.name = Transformers(ii).name;
    DSSCircuit.Transformers.Wdg=1;
    Transformers(ii).kva = get(DSSCircuit.Transformers,'kva');
    if ~isempty(normhkva)
        Transformers(ii).normhkva = str2double(normhkva);
    else
        Transformers(ii).normhkva = Transformers(ii).kva*1.1;
    end
    Transformers(ii).XfmrCode = get(DSSCircuit.Transformers,'XfmrCode');
    Transformers(ii).wdg1R = get(DSSCircuit.Transformers,'R');
    Transformers(ii).wdg1Tap = get(DSSCircuit.Transformers,'Tap');
    Transformers(ii).wdg1minTap = get(DSSCircuit.Transformers,'MinTap');
    Transformers(ii).wdg1maxTap = get(DSSCircuit.Transformers,'MaxTap');
    Transformers(ii).wdg1numTaps = get(DSSCircuit.Transformers,'NumTaps');
    Transformers(ii).bus1kV = get(DSSCircuit.Transformers,'kV');
    Transformers(ii).Xneut = get(DSSCircuit.Transformers,'Xneut');
    Transformers(ii).Rneut = get(DSSCircuit.Transformers,'Rneut');
    Transformers(ii).wdg1IsDelta = get(DSSCircuit.Transformers,'IsDelta');
    Transformers(ii).Xhl = get(DSSCircuit.Transformers,'Xhl');
    Transformers(ii).Xht = get(DSSCircuit.Transformers,'Xlt');
    
    DSSCircuit.Transformers.Wdg=2;
    Transformers(ii).bus2kV = get(DSSCircuit.Transformers,'kV');
    Transformers(ii).wdg2R = get(DSSCircuit.Transformers,'R');
    Transformers(ii).wdg2Tap = get(DSSCircuit.Transformers,'Tap');
    Transformers(ii).wdg2minTap = get(DSSCircuit.Transformers,'MinTap');
    Transformers(ii).wdg2maxTap = get(DSSCircuit.Transformers,'MaxTap');
    Transformers(ii).wdg2numTaps = get(DSSCircuit.Transformers,'NumTaps');
    Transformers(ii).wdg2IsDelta = get(DSSCircuit.Transformers,'IsDelta');
    
    DSSText.command = sprintf('? transformer.%s.%%loadloss',Transformers(ii).name);
    Transformers(ii).PCTloadLoss = str2double(DSSText.result);
    DSSText.command = sprintf('? transformer.%s.%%noloadloss',Transformers(ii).name);
    Transformers(ii).PCTnoLoadLoss = str2double(DSSText.result);
    
    DSSCircuit.Transformers.next;
end

% get transformer control information
XfmControlNames = DSSCircuit.RegControls.AllNames;
DSSCircuit.RegControls.first;
for ii=1:length(XfmControlNames)
    regControlXfmr = get(DSSCircuit.RegControls,'Transformer');
    %XfmControls(ii).Mode = get(DSSCircuit.SwtControls,'Mode');
    [tf, loc] = ismember(regControlXfmr,{Transformers.name});
    if tf
        Transformers(loc).controlled = true;
        Transformers(loc).controller = XfmControlNames{ii};
        %Additional terms to retrieve for the control object
        Transformers(loc).CTPrimary = get(DSSCircuit.RegControls,'CTPrimary');
        Transformers(loc).delay = get(DSSCircuit.RegControls,'Delay');
        Transformers(loc).forwardBand = get(DSSCircuit.RegControls,'ForwardBand');
        Transformers(loc).forwardR = get(DSSCircuit.RegControls,'ForwardR');
        Transformers(loc).forwardVreg = get(DSSCircuit.RegControls,'ForwardVreg');
        Transformers(loc).forwardX = get(DSSCircuit.RegControls,'ForwardX');
        Transformers(loc).isInverseTime = get(DSSCircuit.RegControls,'IsInverseTime');
        Transformers(loc).isReversible = get(DSSCircuit.RegControls,'IsReversible');
        Transformers(loc).maxTapChange = get(DSSCircuit.RegControls,'MaxTapChange');
        Transformers(loc).monitoredBus = get(DSSCircuit.RegControls,'MonitoredBus');
        Transformers(loc).PTratio = get(DSSCircuit.RegControls,'PTRatio');
        Transformers(loc).reverseBand = get(DSSCircuit.RegControls,'ReverseBand');
        Transformers(loc).reverseR = get(DSSCircuit.RegControls,'ReverseR');
        Transformers(loc).reverseVreg = get(DSSCircuit.RegControls,'ReverseVreg');
        Transformers(loc).reverseX = get(DSSCircuit.RegControls,'ReverseX');
        Transformers(loc).tapDelay = get(DSSCircuit.RegControls,'TapDelay');
        Transformers(loc).tapWinding = get(DSSCircuit.RegControls,'TapWinding');
        Transformers(loc).voltageLimit = get(DSSCircuit.RegControls,'VoltageLimit');
        Transformers(loc).winding = get(DSSCircuit.RegControls,'Winding');
    end
    DSSCircuit.RegControls.next;
end


%% As long as you are not in faultstudy mode, remove all transformers that have zero volts on either side (not disabled but are isolated from the circuit)
if ~isempty(Transformers) && isempty(varargin) && ~strcmp(DSSCircuit.Solution.ModeID,'Faultstudy')
    condition = [Transformers.bus1Voltage]>100;
    Transformers = Transformers(condition);
end

catch err
    if ~strcmp(err.identifier,'xfmrName:notfound')
        allTransformers = [err.stack.line];
        allNames = {err.stack.name};
        fprintf(1, ['\nThere was an error in ' allNames{end} ' in line ' num2str(allTransformers(end)) ':\n'])
        fprintf(1, ['"' err.message '"' '\n\n'])
        fprintf(1, ['About to run circuitCheck.m to ensure the circuit is set up correctly in OpenDSS.\n\n'])
        fprintf(1, 'If the problem persists, change the MATLAB debug mode by entering in the command window:\n >> dbstop if caught error\n\n')
        fprintf(1, 'Running circuitCheck.m ....................\n')

        circuitCheck(DSSCircObj, 'Warnings', 'on')
    end
    rethrow(err);
end

end