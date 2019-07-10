import numpy as np
import opendssdirect as dss
import matplotlib.pyplot as plt
import copy
import time
import datetime
# Global variable initialization and error checking



slack_bus_voltage = 1.02
dss.run_command('Compile ' + 'C:/feeders/MultiPhase/13Bus/IEEE13Nodeckt.dss')  # redirecting to the model
# print (os.getcwd())
dss.Vsources.PU(slack_bus_voltage)  # setting up the slack bus voltage
# Setting up the solution parameters, check OpenDSS documentation for details
start_time = time.time()
for i in range(10):
    dss.Solution.Solve()  # solve commands execute the power flow
    dss.Circuit.SetActiveBus('671')
    x=dss.Bus.PuVoltage()
    # print (x)
    y=[]
    for i in range(0,len(x),2):
        y.append(dss.CmathLib.cabs(x[i],x[i+1]))
    print(y)
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

x = np.random.rand(13,13)
y = np.random.rand(13,1)
start_time =  datetime.datetime.now()
print(start_time)
z = np.dot(x,y)
# print (z)
end_time =  datetime.datetime.now()
print (end_time)
# c = end_time-start_time
# print (c)

# print(x)

start_time = time.time()
dss.run_command('Compile ' + 'C:/feeders/MultiPhase/13Bus/IEEE13Nodeckt.dss') 
dss.Solution.Solve()
end_time = time.time()
print(f'Time elapsed is {(end_time-start_time)} Power Flow with networksize of 13 Bus')