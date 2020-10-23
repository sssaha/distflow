%% OpenDSS and GridPV Toolbox Tutorial

%% Initiating the COM Interface

    [DSSCircObj, DSSText, gridpvPath] = DSSStartup;

    % DSSCircObj is passed into the toolbox functions
    %   The Active Circuit and the Text Interface are extracted within the
    %   toolbox functions
    DSSCircuit = DSSCircObj.ActiveCircuit;
    DSSText = DSSCircObj.Text;

    % Note that the circuit is currently empty
    DSSCircuit.Lines.get

%% Getting data out of OpenDSS

    % Load the circuit into OpenDSS
    DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];

    % Now you can start interfacing via the COM interface or start pulling data 
    % out of the circuit structure
    % (Note that we did not have to redefine DSSCircuit after loading)

    %Show the line info in the command interface
    DSSCircuit.Lines.get
    %Get individual fields
    DSSCircuit.Lines.name
    get(DSSCircuit.Lines, 'name')

    DSSCircuit.Lines.next;
    DSSCircuit.Lines.name

    DSSCircuit.Lines.first
    firstLineNm = DSSCircuit.Lines.name
    
 %% Difference b/w the ActiveElement data and the specific object (eg 'Lines') data
 
    DSSCircuit.SetActiveElement(['line.' firstLineNm]);
    DSSCircuit.ActiveElement.get
    DSSCircuit.Lines.get
    DSSCircuit.ActiveElement.currents

    DSSCircuit.SetActiveElement('line.05410_103405ug');
    DSSCircuit.ActiveElement.get
    DSSCircuit.ActiveElement.voltages

 %% Other commands
    
    % Usefule for iterating
    lineNames = DSSCircuit.Lines.AllNames
    DSSCircuit.Lines.Count
    
    % Many generic calls work for most objects
    xfmrNames = DSSCircuit.Transformers.AllNames
    DSSCircuit.Transformers.Count
    DSSCircuit.Capacitors.Count
    
    %Get all elements in the circuit
    DSSCircuit.AllElementNames
    
    % List the Methods associated with a certain object
    DSSCircuit.methods
    DSSCircuit.Lines.methods
    
%% Buses
    %Buses are not in the structure explicitly, but we can still easily
    %access their names, etc
    DSSText.Command = 'Export Buscoords';
    busStruct = importdata(DSSText.Result);
    busNames = busStruct.rowheaders;
    %This is all done in the toolbox using getBusCoordinates(DSSCircObj)
    
%% Adding/Editing Elements in OpenDSS 
    
    % Note that there are currently no generators
    DSSCircuit.Generators.get;

    % Add PV in the form of a generator object
    DSSText.command = 'new generator.PV bus1= n292757 phases=3 	kv=34.5 kw=500 pf=1 enabled=true';

    % You can now see the generator that was added
    DSSCircuit.Generators.get;

    % Set it as the active element and view its bus information
    DSSCircuit.SetActiveElement('generator.pv');
    DSSCircuit.ActiveElement.BusNames

    % Now change it to another bus and observe the change
    DSSText.command = 'edit generator.PV bus1=n1325391 kv=13.2';
    DSSCircuit.ActiveElement.BusNames

    
    DSSText.command = 'Set controlmode=static';
    DSSText.command = 'Set mode=snapshot number=1  hour=0  h=1 sec=0';
    DSSText.command = 'solve';
    % Calling 'solve' again will solve for the next timestep, 'h', making
    % timeseries analysis simpler.
    % However, if you're interested in the same snapshot, be sure to reset
    % the time to the timestep you are interested in before each solve.
    
    % This method of adding only exists for the lifetime of the COM server
    % There are other,  more permanent routes such as editing the .DSS file
    
%% Toolbox

    %% Get Functions
        Buses = getBusInfo(DSSCircObj);
        % {Buses(1:4).name}  <-- Makes a cell array of the first 4 names
        coordinates = getCoordinates(DSSCircObj,{Buses(1:4).name})
        [busCoordNames busCoordArray] = getBusCoordinatesArray(DSSCircObj);
        Lines = getLineInfo(DSSCircObj);
        Transformers = getTransformerInfo(DSSCircObj);
        Loads = getLoadInfo(DSSCircObj);
        PV = getGeneratorInfo(DSSCircObj);
    
    %% Circuit Check
        % The example circuit has been tested with the circuit checker
        % and all warnings have been fixed.
        % Therefore, there will be no warnings returned.
        warnSt = circuitCheck(DSSCircObj);
        
    %% Circuit Analysis Functions
        downstreamBuses = findDownstreamBuses(DSSCircObj,'N292303');
        %FindHighestImpedanceBus is beging replaced by a built-in OpenDSS function
        [longestDistance toBus] = findLongestDistanceBus(DSSCircObj, 'perPhase');
        upstreamBuses = findUpstreamBuses(DSSCircObj,'N292303');
    
    %% Plotting Functions
    % See examples in the function headers for more details on all the
    % options
    
    % Figures are interactive and left and right click funcitons are enabled
    
        % plotCircuitLines
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
        figure; plotCircuitLines(DSSCircObj,'Coloring','voltage','Thickness','current','SubstationMarker','off')
        figure; plotCircuitLines(DSSCircObj,'CapacitorMarker','on')
        figure; plotCircuitLines(DSSCircObj);
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_PV_Central_7_5.dss"'];
        DSSText.command = 'Set mode=duty number=10  hour=13  h=1 sec=1800';
        DSSText.command = 'Set controlmode = static';
        DSSText.command = 'solve';
        figure; plotCircuitLinesOptions(DSSCircObj);
        
        
        % plotVoltageProfile
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
        figure; plotVoltageProfile(DSSCircObj,'BusName','N1311915','Downstream','on');
        figure; plotVoltageProfile(DSSCircObj)
        figure; plotVoltageProfile(DSSCircObj,'Only3Phase','on')
        figure; plotVoltageProfile(DSSCircObj,'SecondarySystem','off','AveragePhase','addition','Only3Phase','on')
                DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_PV_Distributed_7_5.dss"'];
        DSSText.command = 'Set mode=duty number=1  hour=12  h=1 sec=0';
        DSSText.command = 'Set controlmode=static';
        DSSText.command = 'solve';
        figure; plotVoltageProfile(DSSCircObj)
        
        % plotKWProfile
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
        figure; plotKWProfile(DSSCircObj,'AveragePhase','addition','BusName','N1311915');
        figure; plotKWProfile(DSSCircObj,'AveragePhase','on');
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_PV_Central_7_5.dss"'];
        DSSText.command = 'Set mode=duty number=1  hour=12  h=1 sec=0';
        DSSText.command = 'Set controlmode=static';
        DSSText.command = 'solve';
        figure; plotKWProfile(DSSCircObj,'BusName','N300552')
        
        % plotKVARProfile
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
        figure; plotKVARProfile(DSSCircObj,'AveragePhase','addition','BusName','N1311915');
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\Ckt24_PV_Central_7_5.dss"'];
        DSSText.command = 'Set mode=duty number=1  hour=12  h=1 sec=0';
        DSSText.command = 'Set controlmode=static';
        DSSText.command = 'solve';
        figure; plotKVARProfile(DSSCircObj,'BusName','N300552')

        % plotMonitor
        DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
        DSSText.command = 'Set mode=duty number=8760  hour=0  h=1h sec=0';
        DSSText.command = 'Set controlmode = time';
        DSSText.command = 'solve';
        plotMonitor(DSSCircObj,'fdr_05410_Mon_PQ')
        
    %% Coordinate Conversion
        initCoordConversion();
        
    %% Solar Modeling Functions
        placePVplant();
        createPVscenarioFiles();

        

    
    
