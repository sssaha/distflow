import numpy as np
import opendssdirect as dss
from utils.device.Inverter import Inverter
from utils.controller.AdaptiveInvController import AdaptiveInvController
from utils.controller.FixedInvController import FixedInvController
import matplotlib.pyplot as plt
from math import tan, acos
import copy
import pandas as pd
import random
from scipy.interpolate import interp1d
import os 

# Global variable initialization and error checking

mva_base = 1 # mva_base is set to 1 as per unit values are not going to be used, rather kw and kvar are going to be used
load_scaling_factor = 1.5 # scaling factor to tune the loading values
generation_scaling_factor = 2.0 # scaling factor to tune the generation values 
slack_bus_voltage = 1.04 # slack bus voltage, tune this parameter to get a different voltage profile for the whole network
noise_multiplier = 0 # If you want to add noise to the signal, put a positive value
start_time = 42900  # Set simulation analysis period - the simulation is from StartTime to EndTime
end_time = 44000
end_time += 1  # creating a list, last element does not count, so we increase EndTime by 1
# Set hack parameters
time_step_of_hack = 500 # When you want to initiate the hacking time 
# how much hacking we are going to allow, 0.5 means 50% of the invereter capacity can not be changed after the attack is being initiated, Make sure size of this array >= number of inverters
# Choose the time settings accordingly 


# Set initial VBP parameters for un-compromised inverters
VQ_start = 0.98
VQ_end = 1.01
VP_start = 1.02
VP_end = 1.05
hacked_settings = np.array([1.0, 1.001, 1.001, 1.01])
# Set delays for each node

percent_hacked = [0.5, .5, 0.5, 0, .5, .5, .5, 0, .5, 0, 0, .5, 0.5]
Delay_VBPCurveShift = [60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60]

# Set observer voltage threshold
threshold_vqvp = 0.25
power_factor = 0.9
pf_converted = tan(acos(power_factor))
number_of_inverters = 13  # even feeder is 34Bus, we only have 13 inverters, which is chosen randomly for now

# File directory
FileDirectoryBase = '../../../Data Files/testpvnum10/'  # Get the data from the Testpvnum folder
network_model_directory = 'feeder/feeder34_B_NR/feeder34_B_NR.dss'

# Error checking of the global variable
if end_time < start_time or end_time < 0 or start_time < 0:
    print('Setup Simulation Times Inappropriately.')
if noise_multiplier < 0:
    noise_multiplier = 0
    print('Setup Noise Multiplier Correctly.')

# Error checking for percent hacked and the delay value
if len(percent_hacked) < number_of_inverters:
    print('Adjust the percent hacked array.')
    raise SystemError
if len(Delay_VBPCurveShift) < number_of_inverters:
    print('Adjust the delay list accordingly..')
    raise SystemError    
# Global variable initialization done

dss.run_command('Redirect ' + network_model_directory)  # redirecting to the model
# print (os.getcwd())
dss.Vsources.PU(slack_bus_voltage)  # setting up the slack bus voltage
# Setting up the solution parameters, check OpenDSS documentation for details
dss.Monitors.ResetAll()
dss.Solution.Mode(1)
dss.Solution.Number(1)
dss.Solution.StepSize(1)
dss.Solution.ControlMode(-1)
dss.Solution.MaxControlIterations(1000000)
dss.Solution.MaxIterations(30000)
dss.Solution.Solve()  # solve commands execute the power flow
if not dss.Solution.Converged():
    print('Initial Solution Not Converged. Check Model for Convergence')
    raise SystemError
else:
    print('Initial Model Converged. Proceeding to Next Step.')
    total_loads = dss.Loads.Count()
    print ('Load count',total_loads)
    all_load_names = dss.Loads.AllNames()
    print('OpenDSS Model Compilation Done.')


# print ('load Name and KW before change:', dss.Loads.Name(), dss.Loads.kW())
dss.Solution.Solve() 
dss.Loads.First()
dss.Circuit.SetActiveElement('load.' + dss.Loads.Name()) # set active element
dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[0]) # grab the bus for the active element
voltage = dss.Bus.puVmagAngle()[::2] # get the pu information directly
print (f'load Name: {dss.Loads.Name()} KW : {dss.Loads.kW()} Voltage at the load bus: {voltage} before change')
# double the load kw value
dss.Loads.First()
dss.Loads.kW(dss.Loads.kW()*2)
dss.Solution.Solve() 
dss.Circuit.SetActiveElement('load.' + dss.Loads.Name()) # set active element
dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[0]) # grab the bus for the active element
voltage = dss.Bus.puVmagAngle()[::2] # get the pu information directly
print (f'load Name: {dss.Loads.Name()} KW : {dss.Loads.kW()} Voltage at the load bus: {voltage} after change')

# dss.Basic.ClearAll()
# del dss
# import opendssdirect as dss
dss.run_command('Redirect ' + network_model_directory)  # redirecting to the model
dss.Loads.First()
print ('load Name and KW After Compiling:', dss.Loads.Name(), dss.Loads.kW())
dss.Solution.Solve() 
dss.Loads.First()
dss.Circuit.SetActiveElement('load.' + dss.Loads.Name()) # set active element
dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[0]) # grab the bus for the active element
voltage = dss.Bus.puVmagAngle()[::2] # get the pu information directly
print (f'load Name: {dss.Loads.Name()} KW : {dss.Loads.kW()} Voltage at the load bus: {voltage} after recompiling')