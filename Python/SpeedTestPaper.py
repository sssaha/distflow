import numpy as np
import opendssdirect as dss
import matplotlib.pyplot as plt
import copy
import time
# Global variable initialization and error checking



slack_bus_voltage = 1.02
dss.run_command('Compile ' + 'C:/feeders/MultiPhase/13Bus/IEEE13Nodeckt.dss')  # redirecting to the model
# print (os.getcwd())
dss.Vsources.PU(slack_bus_voltage)  # setting up the slack bus voltage
# Setting up the solution parameters, check OpenDSS documentation for details
start_time = time.time()
for i in range(200):
    dss.Solution.Solve()  # solve commands execute the power flow
    # print('Iterations:',  dss.Solution.Iterations())
end_time = time.time()
if not dss.Solution.Converged():
    print('Initial Solution Not Converged. Check Model for Convergence')
    raise SystemError
else:
    print('Model Converged.')
    print(f'Time elapsed is {(end_time-start_time)} Power Flow with networksize of 13 Bus')
    total_loads = dss.Loads.Count()
    print ('Load count',total_loads)
    print('OpenDSS Model Compilation Done.')
    print ('\n')




slack_bus_voltage = 1.02
dss.run_command('Compile ' + 'C:/feeders/MultiPhase/34Bus/ieee34Mod1.dss')  # redirecting to the model
# print (os.getcwd())
dss.Vsources.PU(slack_bus_voltage)  # setting up the slack bus voltage
# Setting up the solution parameters, check OpenDSS documentation for details
start_time = time.time()
for i in range(200):
    dss.Solution.Solve()  # solve commands execute the power flow
    # print('Iterations:',  dss.Solution.Iterations())
end_time = time.time()
if not dss.Solution.Converged():
    print('Initial Solution Not Converged. Check Model for Convergence')
    raise SystemError
else:
    print('Model Converged.')
    print(f'Time elapsed is {(end_time-start_time)} Power Flow with networksize of 34 Bus')
    total_loads = dss.Loads.Count()
    print ('Load count',total_loads)
    print('OpenDSS Model Compilation Done.')
    print ('\n')


slack_bus_voltage = 1.02
dss.run_command('Compile ' + 'C:/feeders/MultiPhase/37Bus/ieee37.dss')  # redirecting to the model
# print (os.getcwd())
dss.Vsources.PU(slack_bus_voltage)  # setting up the slack bus voltage
# Setting up the solution parameters, check OpenDSS documentation for details
start_time = time.time()
for i in range(200):
    dss.Solution.Solve()  # solve commands execute the power flow
    # print('Iterations:',  dss.Solution.Iterations())
end_time = time.time()
if not dss.Solution.Converged():
    print('Initial Solution Not Converged. Check Model for Convergence')
    raise SystemError
else:
    print('Model Converged.')
    print(f'Time elapsed is {(end_time-start_time)} Power Flow with networksize of 37 Bus')
    total_loads = dss.Loads.Count()
    print ('Load count',total_loads)
    all_load_names = dss.Loads.AllNames()
    print('OpenDSS Model Compilation Done.')
    print ('\n')


slack_bus_voltage = 1.02
dss.run_command('Compile ' + 'C:/feeders/MultiPhase/123Bus/IEEE123Master.dss')  # redirecting to the model
# print (os.getcwd())
dss.Vsources.PU(slack_bus_voltage)  # setting up the slack bus voltage
# Setting up the solution parameters, check OpenDSS documentation for details
start_time = time.time()
for i in range(200):
    dss.Solution.Solve()  # solve commands execute the power flow
    # print('Iterations:',  dss.Solution.Iterations())
end_time = time.time()
if not dss.Solution.Converged():
    print('Initial Solution Not Converged. Check Model for Convergence')
    raise SystemError
else:
    print('Model Converged.')
    print(f'Time elapsed is {(end_time-start_time)} Power Flow with networksize of 123 Bus')
    total_loads = dss.Loads.Count()
    print ('Load count',total_loads)
    all_load_names = dss.Loads.AllNames()
    print('OpenDSS Model Compilation Done.')
    print ('\n')