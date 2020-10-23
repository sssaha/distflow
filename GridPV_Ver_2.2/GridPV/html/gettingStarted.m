%% Getting Started
%
%% Introduction
% The power industry is seeing large amounts of distributed generation being added onto the electric power distribution system.  This presents a new set of issues, especially for renewable generation with variability.  This presents a new set of issues, especially for renewable generation with variable intermittent power output.  It is important to precisely model the impact of solar energy on the grid and to help distribution planners perform the necessary interconnection impact studies.  The variability in the load, throughout the day and year, and the variability of solar, throughout the year and because of clouds, makes the analysis increasingly complex.  Both accurate data and timeseries simulations are required to fully understand the impact of variability on distribution system operations and reliability.  
% 
% This help manual describes the functionality and use of a MATLAB toolbox for using OpenDSS to model the variable nature of the distribution system load and solar energy.  OpenDSS is used to model the distribution system with MATLAB providing the frontend user interface through a COM interface.  OpenDSS is designed for distribution system analysis and is very good at timeseries analysis with changing variables and dynamic control.  OpenDSS is command based and has limited visualization capabilities.  By bringing control of OpenDSS to MATLAB, the functionality of OpenDSS is utilized while adding the looping, advanced analysis, and visualization abilities of MATLAB.
%
% The functions in the toolbox are categorized into five main sections: OpenDSS functions, Solar Modeling functions, Plotting functions, Geographic Mapping functions, and Example Simulations.  Each function is documented with the function use syntax, full description, function input list, function output list, and an example use.  The function example also includes an example output of the function.

%% Objectives
% The GridPV Toolbox for Matlab provides a set of well-documented functions for simulating the performance of photovoltaic energy systems.
% 
% The toolbox was developed at Georgia Institute of Technology and Sandia National Laboratories.  It implements many of the models and methods developed at the Labs.  Future versions are planned that will add more functions and capability.

%% What is OpenDSS?
% The distribution system electrical modeling is done in the open source software OpenDSS from the Electric Power Research Institute (EPRI). OpenDSS is commonly used to model solar on the grid because of its high-resolution time series analysis capabilities. It is a 3-phase distribution system analysis power flow solver that can handle unbalanced phases.

%% Initial Steps
% Each toolbox function has its own example contained in the header file. These examples will run on their own using the example circuit and may be useful for becoming familiar with the toolbox.
%
% The basic process for getting started with the toolbox is:

% 1. Start the OpenDSS COM. Needs to be done each time MATLAB is opened
    [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% 2. Compiling the circuit
    DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% 3. Solve the circuit. Call anytime you want the circuit to resolve
    DSSText.command = 'solve';
% 4. Run circuitCheck function to double-check for any errors in the circuit before using the toolbox
    warnSt = circuitCheck(DSSCircObj);
