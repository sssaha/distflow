import numpy as np
import opendssdirect as dss
import matplotlib.pyplot as plt
import copy
import time
import datetime
import os 
from DSSStartup import DSSStartup


result= DSSStartup()
dssText=result['dsstext']
dssSolution=result['dsssolution']
dssCircuit=result['dsscircuit']
DSSObj=result['dssobj']
# Global variable initialization and error checking



slack_bus_voltage = 1.02
# print(os.getcwd())
# os.chdir('C:/ceds-cigar-external/PyCIGAR')
# print(os.getcwd())
dssText.Command='Compile IEEE13Nodeckt.dss'  # redirecting to the model
# Setting up the solution parameters, check OpenDSS documentation for details
start_time = time.time()
for i in range(200):
    dssSolution.Solve()  # solve commands execute the power flow
    
    # print('Iterations:',  dss.Solution.Iterations())
end_time = time.time()

if not dssSolution.Converged:
    print('Initial Solution Not Converged. Check Model for Convergence')
    raise SystemError
else:
    print('Model Converged.')
    print(f'Time elapsed is {(end_time-start_time)} Power Flow with networksize of 13 Bus')
    total_loads = dss.Loads.Count()
    print ('Load count',total_loads)
    print('OpenDSS Model Compilation Done.')
    print ('\n')

# x = np.random.rand(13,13)
# y = np.random.rand(13,1)
# start_time =  datetime.datetime.now()
# print(start_time)
# z = np.dot(x,y)
# # print (z)
# end_time =  datetime.datetime.now()
# print (end_time)
# # c = end_time-start_time
# # print (c)

# # print(x)

# start_time = time.time()
# dss.run_command('Compile ' + 'C:/feeders/MultiPhase/13Bus/IEEE13Nodeckt.dss') 
# dss.Solution.Solve()
# end_time = time.time()
# print(f'Time elapsed is {(end_time-start_time)} Power Flow with networksize of 13 Bus')