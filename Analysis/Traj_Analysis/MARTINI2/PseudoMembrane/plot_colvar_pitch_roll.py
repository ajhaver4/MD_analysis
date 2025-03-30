import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# Read the CSV file using pandas
data = pd.read_csv(sys.argv[1],delimiter='\s+',comment='#',header=None,names=['time', 'pitchA', 'pitchB', 'rollA', 'rollB', 
                                                                              'res_pitchA.bias', 'res_pitchB.bias', 'res_rollA.bias', 
                                                                              'res_rollB.bias'])
fig,ax=plt.subplots(2,1)
# Plot the data using matplotlib
ax[0].plot(data['time'],data['pitchA'],color='crimson',label='chain 1')
ax[0].plot(data['time'],data['pitchB'],color='steelblue',label='chain 2')

ax[1].plot(data['time'],data['rollA'],color='crimson',label='chain 1')
ax[1].plot(data['time'],data['rollB'],color='steelblue',label='chain 2')

ax[0].set_xlabel('Time')
ax[0].set_ylabel('Radians ')
ax[0].set_title('Plot of Pitch Angle')
ax[0].legend()


ax[1].set_xlabel('Time')
ax[1].set_ylabel('Radians ')
ax[1].set_title('Plot of Roll Angle')
ax[1].legend()



fig,ax = plt.subplots()

ax.hist(data['pitchA'],bins=100,color='crimson',alpha=0.5,label='chain 1',density=True)
ax.hist(data['pitchB'],bins=100,color='steelblue',alpha=0.5,label='chain 2',density=True)

ax.set_xlabel('Radians')
ax.set_ylabel('Frequency')
ax.set_title('Histogram of Pitch')
ax.legend()

fig,ax = plt.subplots()

ax.hist(data['rollA'],bins=100,color='crimson',alpha=0.5,label='chain 1',density=True)
ax.hist(data['rollB'],bins=100,color='steelblue',alpha=0.5,label='chain 2',density=True)

ax.set_xlabel('Radians')
ax.set_ylabel('Frequency')
ax.set_title('Histogram of Roll')
ax.legend()
plt.show()




