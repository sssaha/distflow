%% Functions by Category
% GridPV Toolbox
%
%
%% OpenDSS Functions
% * <DSSStartup_help.html |DSSStartup|> - Function for starting up OpenDSS and linking to MATLAB
% * <getBusCoordinatesArray_help.html |getBusCoordinatesArray|> - Gets the coordinates for all buses that have a location in OpenDSS
% * <getBusInfo_help.html |getBusInfo|> - Gets the information for all Bus in busNames
% * <getCapacitorInfo_help.html |getCapacitorInfo|> - Gets the information for all capacitors in the circuit
% * <getCoordinates_help.html |getCoordinates|> - Gets the coordinates for the buses in busNames
% * <getGeneratorInfo_help.html |getGeneratorInfo|> - Gets the information for all generators in the circuit
% * <getLineInfo_help.html |getLineInfo|> - Gets the information for all lines in the circuit
% * <getLoadInfo_help.html |getLoadInfo|> - Gets the information for all loads in the circuit
% * <getPVInfo_help.html |getPVInfo|> - Gets the information for all PV plants in the circuit
% * <getTransformerInfo_help.html |getTransformerInfo|> - Gets the information for all transformers in the circuit
% * <isinterfaceOpenDSS_help.html |isinterfaceOpenDSS|> - Used to check for a valid interface input.
%
%% Circuit Analysis Functions
% * <circuitCheck_help.html |circuitCheck|> - Used to error-check the circuit for any obvious abnormalities
% * <findDownstreamBuses_help.html |findDownstreamBuses|> - Finds all buses downstream of the busName
% * <findHighestImpedanceBus_help.html |findHighestImpedanceBus|> - Finds the highest impedance bus for each phase to the source bus
% * <findLongestDistanceBus_help.html |findLongestDistanceBus|> - Finds the bus for each phase that is farthest distance away
% * <findSubstationLocation_help.html |findSubstationLocation|> - Locates the substation coordinates
% * <findUpstreamBuses_help.html |findUpstreamBuses|> - Finds all buses upstream of the busName
%
%% Plotting Functions
% * <plotAmpProfile_help.html |plotAmpProfile|> - Plots the line currents profile and line rating vs. distance
% * <plotCircuitLines_help.html |plotCircuitLines|> - Plots the feeder circuit diagram
% * <plotCircuitLinesOptions_help.html |plotCircuitLinesOptions|> - GUI for providing options for how to plot the feeder circuit diagram
% * <plotKVARProfile_help.html |plotKVARProfile|> - Plots the feeder profile for the kVAR power flow on the lines
% * <plotKWProfile_help.html |plotKWProfile|> - Plots the feeder profile for the kW power flow on the lines
% * <plotMonitor_help.html |plotMonitor|> - Plots a monitor from the simulation
% * <plotVoltageProfile_help.html |plotVoltageProfile|> - Plots the voltage profile for the feeder (spider plot)
%
%% Geographic Mapping Functions
% * <initCoordConversion_help.html |initCoordConversion|> - Function to initialize the coordinate conversion process
% * <createCircuitCoordConversion_help.html |createCircuitCoordConversion|> - Function to create conversion of circuit coordinates to GPS coordinates
% * <createCircuitCoordConversionUTM_help.html |createCircuitCoordConversionUTM|> - Function to create conversion of circuit coordinates in UTM to GPS coordinates
% * <plotGoogleMap_help.html |plotGoogleMap|> - Plots a Google map on the current axes using the Google Static Maps API
%
%% Solar Modeling Functions
% * <placePVplant_help.html |placePVplant|> - Draw PV on the circuit diagram and save plant info for WVM input
% * <createPVscenarioFiles_help.html |createPVscenarioFiles|> - Runs the WVM model and puts out the OpenDSS PV scenario files
% * <distributePV_help.html |distributePV - Allocates|> PV based off of the load transformer size (kva)
% * <findMaxPenetrationTime_help.html |findMaxPenetrationTime|> - Finds the max penetration time
% * <IneichenClearSkyModel_help.html |IneichenClearSkyModel|> - Generates the clear sky irradiance using Ineichen and Perez model
% * <makePFoutputFunction_help.html |makePFoutputFunction|> - GUI for creating power factor as a function of PV power output
% * <makePFprofile_help.html |makePFprofile|> - Creates varying Power Factor profile by schedule or PV output
% * <makePFschedule_help.html |makePFschedule|> - GUI for creating a power factor daily schedule
% * <makeVVCcurve_help.html |makeVVCcurve|> - GUI for setting up the OpenDSS VVControl function parameters
% * <pvl_WVM_help.html |pvl_WVM|> - WVM Wavelet Variability Model
%
%% Example Simulations
% * <examplePeakTimeAnalysis_help.html |examplePeakTimeAnalysis|> - Runs simulation during peak penetration time and generates plots
% * <exampleTimeseriesAnalyses_help.html |exampleTimeseriesAnalyses|> - Timeseries analysis and plots monitor values from the simulation
% * <exampleVoltageAnalysis_help.html |exampleVoltageAnalysis|> - Example analysis of maximum and minimum feeder voltages through time
%


