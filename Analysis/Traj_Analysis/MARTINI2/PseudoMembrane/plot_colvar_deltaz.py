import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# Read the CSV file using pandas
data = pd.read_csv(sys.argv[1],delimiter='\s+',comment='#',header=None,names=['time','p1.z','p2.z','r1.bias','r2.bias'])

fig,ax=plt.subplots()
# Plot the data using matplotlib
ax.plot(data['time'],data['p1.z'],color='crimson',label='chain 1')
ax.plot(data['time'],data['p2.z'],color='steelblue',label='chain 2')

ax.set_xlabel('Time')
ax.set_ylabel('Position (nm)')
ax.set_title('Plot of Delta Z COLVAR')
ax.legend()
plt.show()


fig,ax = plt.subplots()

ax.hist(data['p1.z'],bins=100,color='crimson',alpha=0.5,label='chain 1')
ax.hist(data['p2.z'],bins=100,color='steelblue',alpha=0.5,label='chain 2')

ax.set_xlabel('Position (nm)')
ax.set_ylabel('Frequency')
ax.set_title('Histogram of Delta Z COLVAR')
ax.legend()
plt.show()




