import sys
import numpy as np
import pandas

print("Reading file...", sys.argv[1])
col_names=['Timestep','Bond','HP','G96','PDih','ImPDih','LJ','Col','Pot','Col-Pro-Pro','LJ-Pro-Pro','Col-Pro-Mem','LJ-Pro-Mem','Col-Pro-W','LJ-Pro-W','Col-Mem-Mem','LJ-Mem-Mem','Col-Mem-W','LJ-Mem-W','Col-W-W','LJ-W-W']
ener_data = pandas.read_csv(sys.argv[1],comment='#',delimiter="\s+",names=col_names)

print("Reading file...", sys.argv[2])
rw_data = pandas.read_csv(sys.argv[2],comment='#',delimiter="\s+",names=['Time','d1','d2','bias','rct','rbias'])

print("Reading file..", sys.argv[3])
vol_data = pandas.read_csv(sys.argv[3],comment='#',delimiter="\s+",names=['Timestep','Volume'])

print("Total length: ",len(rw_data['Time']))
print("Total length: ",len(ener_data['Timestep']))


energy_ts = 200   # 20ps interval of energies
colvar_ts = 200  #20 ps time interval for COLVARS

filtered_ener = ener_data.drop_duplicates(subset='Timestep',ignore_index=True)
filtered_vol = vol_data.drop_duplicates(subset='Timestep',ignore_index=True)
print("Non duplicates length Energy : ",len(filtered_ener['Timestep']))
print("Non duplicates length Volume: ",len(filtered_vol['Timestep']))

join_data = filtered_ener.join([rw_data,filtered_vol['Volume']])

print(join_data)


count=300001

ntimes=int(len(rw_data['Time'])/count)


print("Lenght of FINAL array: ",len(filtered_ener['Timestep']))

print("Writing to a file...")
join_data.to_csv("Reweigthed_Data_FINAL",sep='\t')
