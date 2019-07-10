import numpy as np
import opendssdirect as dss
import matplotlib.pyplot as plt
import copy
import time
# Global variable initialization and error checking

#######################################################################################

slack_bus_voltage = 1.02
dss.run_command('Compile ' + 'C:/feeders/MultiPhase/123BusQSTS/IEEE123Master.dss')  # redirecting to the model
# print (os.getcwd())
dss.Vsources.PU(slack_bus_voltage)  # setting up the slack bus voltage
# Setting up the solution parameters, check OpenDSS documentation for details
all_voltage = []
start_time = time.time()
for i in range(8760):
    dss.Solution.Solve()  # solve commands execute the power flow
    all_voltage.append(dss.Circuit.AllBusMagPu())
end_time = time.time()
if not dss.Solution.Converged():
    print('Initial Solution Not Converged. Check Model for Convergence')
    raise SystemError
else:
    print('Model Converged.')
    print(f'Time elapsed is {(end_time-start_time)} QSTS Power Flow with networksize of 123 Bus')
    total_loads = dss.Loads.Count()
    print ('Load count',total_loads)
    print('OpenDSS Model Compilation Done.')
    print ('\n')


####################################################################    