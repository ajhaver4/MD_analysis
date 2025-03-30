import sys
import pandas
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle as markerstyle
from scipy.interpolate import griddata
import numpy as np
import matplotlib.cm as cm

clust_file = sys.argv[1]
data_file = sys.argv[2]
plot_on_fes = False

if sys.argv[3] is not None:
    plot_on_fes = sys.argv[3]
    fes_file = sys.argv[4]
clust_data = pandas.read_csv(clust_file,comment='#',sep='\t')
master_data = pandas.read_csv(data_file,comment='#',sep='\t')


cols = clust_data.columns
for col in cols:
    clust_data[col] = clust_data[col].astype(float)

remapped_clust = clust_data.copy()

count=0
print(type(clust_data['d1'][0]),type(clust_data['d2']))
clust_dict_dist = {}
flag=False
with open(clust_file,'r') as fl:
    data = fl.read()
    for line in data.split("\n"):
        if line:
            # print(line)
            if line[0]=='#':
                if 'cluster:' in line:
                    stable_state = int(line.split(":")[1])
                    print(stable_state)
                continue


label_colours = {0:'forestgreen', 1:'teal', 2:'crimson', 3:'gold', 4:'orchid', 5:'peru', 6:'mediumpurple', 7:'darkorange',8:'steelblue',-1:'olivedrab'}
fig,ax = plt.subplots()
time_dict={}
for i in range(len(clust_data['Cluster'])):
    if clust_data['Cluster'][i] == 0 or clust_data['Cluster'][i] == 1 or clust_data['Cluster'][i] == 2 \
        or clust_data['Cluster'][i] == 4 or clust_data['Cluster'][i]==6 or clust_data['Cluster'][i] == 7:
        time_dict[clust_data['Timestep'][i]] = clust_data['Cluster'][i]
        
    elif clust_data['Cluster'][i] == 3:
        if clust_data['d1'][i] < 2.5 and clust_data['d2'][i] < 2.5:
            time_dict[clust_data['Timestep'][i]] = -1
            
            remapped_clust['Cluster'][i] = -1

        else:
            time_dict[clust_data['Timestep'][i]] = 4
            remapped_clust['Cluster'][i] = 4
    elif clust_data['Cluster'][i] == 5:
        if clust_data['d2'][i] >= 3.0 :
            time_dict[clust_data['Timestep'][i]] = 2
            remapped_clust['Cluster'][i] = 2
        else:
            time_dict[clust_data['Timestep'][i]] = 0
            remapped_clust['Cluster'][i] = 0


n_values = len(master_data['d1'])

cry_timesteps=[]
# sys.exit()
for i in range(0,n_values):
    if master_data['RMSD-B'][i] < 1.5:
    # if master_data['d1'][i] <= 2.0 or master_data['d2'][i]<=2.0:
        cry_timesteps.append(master_data['Timestep'][i])

cry_data = master_data.loc[master_data['Timestep'].isin(cry_timesteps)]
cry_data.insert(0,"Cluster",-1)
print(cry_data)

#Adding the Crystal state timesteps from the master data
for i in range(len(cry_timesteps)):
    time_dict[cry_timesteps[i]] = -1

remapped_clust = remapped_clust.append(cry_data,ignore_index=True)

clus_values = list(set(remapped_clust['Cluster'].values))
num_clusters=len(clus_values)
print("CLUSTER VALUES: ",clus_values)

if plot_on_fes:
    print("Reading data...")
    data = [line.split() for line in open(fes_file, "r")]
    data2 = [x for x in data if not x == []]  # strip out headers
    d1, d2, free, dd1, dd2 = [], [], [], [], []
    for elem in data2[9:]:
        if elem[2] == 'inf' or elem[2] == '-inf':
            free.append(0.0)
        else:
            free.append(float(elem[2]))

        d1.append(float(elem[0]))
        d2.append(float(elem[1]))
    #    dd1.append(float(elem[3]))
    #    dd2.append(float(elem[4]))

    X = np.linspace(min(d1), max(d1), 1318)
    Y = np.linspace(min(d2), max(d2), 1322)

    print("Creating data grid. This may take a while...")
    D1, D2 = np.meshgrid(X, Y)

    #Normalize unbound state to zero
    free = np.array(free)
    # free = free+110.6

    #Zero value
    d1_arr = np.array(d1)
    d2_arr = np.array(d2)
    mask1 = (d1_arr >= 18.0 ) & (d1_arr < 19.0 )
    mask2 = (d2_arr >= 18.0 ) & (d2_arr < 19.0 )
    mask3 = mask1 * mask2

    correction = np.mean(free[mask3])
    free = free - correction

    #Shift max value to a constant. To shifht the color map
    max_val = 50.0
    mask4 = (free>=max_val)
    free[mask4]=max_val

    ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)
    print(np.max(free))
    print(np.min(free))
    levels = np.arange(np.min(free), np.max(free), (np.max(free)-np.min(free))/25)
    levels = levels.tolist()
    levels = levels
    levels = np.array(levels)
    print(levels)
    print(len(levels))
    contour = plt.contour(D1, D2, ENER, colors='k', linewidths=0.3, levels=levels)
    contourf = plt.contourf(D1, D2, ENER, cmap=cm.Greys,levels=levels,alpha=0.98)
    plt.colorbar(label="Free Energy (kJ/mol)")
    plt.scatter([0.8525], [0.7766], marker='x', c='k')
    plt.xlabel(r"$d_1$ (nm)",fontdict={'FontSize':20})
    plt.ylabel(r"$d_2$ (nm)",fontdict={'FontSize':20})
    ax=plt.gca()
    ax.set_xlim(left=0.5,right=19.0)
    ax.set_ylim(bottom=0.5,top=19.0)


    ax.tick_params(axis='both',labelsize='xx-large')



for clus in clus_values:
    plot_data = remapped_clust.loc[remapped_clust['Cluster']==clus]
    if clus==-1:
        lab='Crystal'
        ax.scatter(plot_data['d1'],plot_data['d2'],c=label_colours[int(clus)],label=lab,s=5)
    else:
        lab=clus
        ax.scatter(plot_data['d1'],plot_data['d2'],c=label_colours[int(clus)],label='C-'+str(int(lab)),s=5)

ax.legend(fancybox=True,framealpha=0.4,fontsize='large',handletextpad=0.1,markerscale=2.5,loc='best',ncol=num_clusters+1,columnspacing=0.5)
ax.set_title("Remapped Clustering")
ax.set_xlabel('d1 (nm)',fontsize='x-large')
ax.set_ylabel('d2 (nm)',fontsize='x-large')



plt.show()




with open("Clust_times.py",'w') as fl1:
    fl1.write("time_dict = {}")
    fl1.write("\n")
    for time,st in time_dict.items():
        count+=1
        fl1.write("time_dict[{:.1f}] = {}".format(time,st))
        fl1.write("\n")




