import numpy as np 
import pandas
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.collections import LineCollection

file=sys.argv[1]
data = pandas.read_csv(file,comment='#',sep='\t')

timesteps_array = data['Timestep']/1e6

#Timepoints to ignore 
 #33.0 - 66.0 us 
#95.2 - 97.8 us
#101.95 to 103.3 us 
#127.6 to 128.5 us

mask1 = (timesteps_array > 33.0) & (timesteps_array < 66.0)
mask2 = (timesteps_array > 95.2) & (timesteps_array < 97.8)
mask3 = (timesteps_array > 101.95) & (timesteps_array < 103.3)
mask4 = (timesteps_array > 127.6) & (timesteps_array < 128.5)
mask5 = (timesteps_array > 139) & (timesteps_array < 142.5)


mask = mask1 | mask2 | mask3 | mask4 | mask5

timesteps_array = timesteps_array[~mask]
data = data[~mask]

disp_ts = 100   #ns
interval = int(disp_ts*1e3 / 200)  #200 ps time interval between data points in the data file
filter_data = data[::interval] 


bd_unb_mask = (filter_data['d1'] < 12.0) & (filter_data['d2'] < 12.0)
mask_bound = filter_data['RMSD-BtoA'] < 2.0

def get_segments(xdata,ydata):

    points = np.array([xdata, ydata]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

t_size = 45
fsize= 48
lw = 2.5
msize= 17
fig_size=(10,8)


# top=0.945
# bottom=0.16
# left=0.09
# right=0.97
# hspace=0.2
# wspace=0.05

top=0.965
bottom=0.175
left=0.105
right=0.965
hspace=0.1
wspace=0.1

def plot_d1d2(lc,bd_series,bd_bool,col, clr, lbl,axes,y_lbl, single_plot=False):
    ax = axes[0]
    ax1 = axes[1]
    ax2 = axes[2]
    ax3 = axes[3]
    ax4 = axes[4]
    ax5 = axes[5]

    
    



    #Plotting as 1 series
    if single_plot:
        series1 = lc
        ax.plot(series1['Timestep']/1e6,series1[col],color=clr,linewidth=lw,label=lbl,linestyle='-')
        ax1.plot(series1['Timestep']/1e6,series1[col],color=clr,linewidth=lw,label=lbl,linestyle='-')
        ax2.plot(series1['Timestep']/1e6,series1[col],color=clr,linewidth=lw,label=lbl,linestyle='-')
        ax3.plot(series1['Timestep']/1e6,series1[col],color=clr,linewidth=lw,label=lbl,linestyle='-')
        ax4.plot(series1['Timestep']/1e6,series1[col],color=clr,linewidth=lw,label=lbl,linestyle='-')
        ax5.plot(series1['Timestep']/1e6,series1[col],color=clr,linewidth=lw,label=lbl,linestyle='-')
    
    else:
        #Plotting as line collection of segements
        ax.add_collection(lc[0])
        ax1.add_collection(lc[1])
        ax2.add_collection(lc[2])
        ax3.add_collection(lc[3])
        ax4.add_collection(lc[4])
        ax5.add_collection(lc[5])

    if bd_bool:
        ax.plot(bd_series['Timestep']/1e6,bd_series[col],'black',linewidth=lw,linestyle='',marker='*',markersize=msize)
        ax1.plot(bd_series['Timestep']/1e6,bd_series[col],'black',linewidth=lw,linestyle='',marker='*',markersize=msize)
        ax2.plot(bd_series['Timestep']/1e6,bd_series[col],'black',linewidth=lw,linestyle='',marker='*',markersize=msize)
        ax3.plot(bd_series['Timestep']/1e6,bd_series[col],'black',linewidth=lw,linestyle='',marker='*',markersize=msize)
        ax4.plot(bd_series['Timestep']/1e6,bd_series[col],'black',linewidth=lw,linestyle='',marker='*',markersize=msize)
        
        ax5.plot(bd_series['Timestep']/1e6,bd_series[col],'black',linewidth=lw,linestyle='',marker='*',markersize=msize)
    
    ax.set_xlim(0, 33)
    ax1.set_xlim(66, 95.2)
    ax2.set_xlim(97.8,101.95)
    ax3.set_xlim(103.3,127.6)
    ax4.set_xlim(127.6,139)
    ax5.set_xlim(142.5,165)


    # hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)

    ax4.spines['right'].set_visible(False)
    ax4.spines['left'].set_visible(False)

    ax5.spines['left'].set_visible(False)


    ax.yaxis.tick_left()
    # ax.tick_params(labelright='off')
    ax5.yaxis.tick_right()

    # This looks pretty good, and was fairly painless, but you can get that
    # cut-out diagonal lines look with just a bit more work. The important
    # thing to know here is that in axes coordinates, which are always
    # between 0-1, spine endpoints are at these locations (0, 0), (0, 1),
    # (1, 0), and (1, 1).  Thus, we just need to put the diagonals in the
    # appropriate corners of each of our axes, and so long as we use the
    # right transform and disable clipping.
    4
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False) # switch to the bottom axes
    ax1.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax1.plot((-d, +d), (-d, +d), **kwargs)
    ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    ax.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)

    ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)
    ax2.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs = dict(transform=ax3.transAxes, color='k', clip_on=False)
    ax3.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax3.plot((-d, +d), (-d, +d), **kwargs)
    ax3.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax3.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs = dict(transform=ax4.transAxes, color='k', clip_on=False)
    ax4.plot((-d, +d), (1-d, 1+d), **kwargs)     #Left side top axis
    ax4.plot((-d, +d), (-d, +d), **kwargs)       #Left side bottom axis
    ax4.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax4.plot((1-d, 1+d), (1-d, 1+d), **kwargs)


    kwargs = dict(transform=ax5.transAxes, color='k', clip_on=False)
    ax5.plot((-d, +d), (1-d, 1+d), **kwargs)     #Left side top axis
    ax5.plot((-d, +d), (-d, +d), **kwargs)    

    # What's cool about this is that now if we vary the distance between
    # ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
    # the diagonal lines will move accordingly, and stay right at the tips
    # of the spines they are 'breaking'


    plt.subplots_adjust(wspace=0.05)

    # # ax.plot(unbound_data['Timestep']/1e6,unbound_data['d2'],'royalblue',linewidth=2.0,label='d2',linestyle='-')

    # # ax.plot(bound_data['Timestep']/1e6,bound_data['d1'],'gold',linewidth=2.5,label='d1',linestyle='',marker='.')
    # # ax.plot(bound_data['Timestep']/1e6,bound_data['d2'],'gold',linewidth=2.5,label='d2',linestyle='',marker='.')

    # ax.legend(fontsize=30)
    # ax.set_xlabel(r'Time ($ \mu s$)',fontsize=40)
    ylbl = y_lbl + (' (nm)')
    ax.set_ylabel(ylbl,fontsize=fsize)
    # #ax_rmsd.set_title("RMSD")
    ax.tick_params(labelsize=t_size)
    ax1.tick_params(labelsize=t_size)
    ax2.tick_params(labelsize=t_size)
    ax3.tick_params(labelsize=t_size)
    ax4.tick_params(labelsize=t_size)
    ax5.tick_params(labelsize=t_size)

    ax1.tick_params(axis='y',left=False)
    ax2.tick_params(axis='y',left=False)
    ax3.tick_params(axis='y',left=False)
    ax4.tick_params(axis='y',left=False)

    # ticks= [0,40,80,120,160,200,240]
    ax.set_xticks([0,30])
    ax1.set_xticks([80])
    ax2.set_xticks([100])
    ax3.set_xticks([120])
    ax4.set_xticks([132])
    ax5.set_xticks([150,250])

    ax.set_yticks([0,10,20,30])


    ax3.set_xlabel(r'Time ($ \mu s$)',fontsize=fsize)
    # ax4.legend(fontsize=25)


bd_series = filter_data[mask_bound]

"""
Plot d1
"""
color_map = np.where(bd_unb_mask , 'firebrick', 'silver')
line_maps = np.where(bd_unb_mask , 'solid', 'solid')
alpha_map = np.where(bd_unb_mask , 0.8, 0.8)
segments_d1 = get_segments(filter_data['Timestep']/1e6,filter_data['d1'])
lc_d1 = [LineCollection(segments_d1, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),LineCollection(segments_d1, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps)
         ,LineCollection(segments_d1, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),
         LineCollection(segments_d1, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),
         LineCollection(segments_d1, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),
         LineCollection(segments_d1, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps)]
bd_series = filter_data[mask_bound]

fig,ax_d1 = plt.subplots(1,6,sharey=True,figsize=fig_size)

plot_d1d2(lc_d1,bd_series,True,'d1', 'crimson', r'$d_1$',ax_d1, r'$d_1$')
plt.subplots_adjust(top=top,bottom=bottom,left=left,right=right,hspace=hspace,wspace=wspace)
plt.show()
fig.savefig("d1_FINAL.svg",dpi=600,format="svg",transparent=True)
fig.savefig("d1_FINAL.png",dpi=600,format="png",transparent=True)

"""
Plot d1
"""
color_map = np.where(bd_unb_mask , 'firebrick', 'silver')
line_maps = np.where(bd_unb_mask , 'solid', 'solid')
alpha_map = np.where(bd_unb_mask , 0.8, 0.8)
segments_d2 = get_segments(filter_data['Timestep']/1e6,filter_data['d2'])
lc_d2 = [LineCollection(segments_d2, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),LineCollection(segments_d2, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),
         LineCollection(segments_d2, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),
         LineCollection(segments_d2, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),
         LineCollection(segments_d2, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps),
         LineCollection(segments_d2, colors=color_map, linewidth=lw,alpha=0.8,linestyle=line_maps)]
fig,ax_d1 = plt.subplots(1,6,sharey=True,figsize=fig_size)
plot_d1d2(lc_d2,bd_series,True,'d2', 'blue', r'$d_2$',ax_d1,r'$d_2$')
plt.subplots_adjust(top=top,bottom=bottom,left=left,right=right,hspace=hspace,wspace=wspace)
plt.show()
fig.savefig("d2_FINAL.svg",dpi=600,format="svg",transparent=True)
fig.savefig("d2_FINAL.png",dpi=600,format="png",transparent=True)


# fig1,ax1 = plt.subplots(1,6,sharey=True,figsize=(16,8))
# plot_d1d2(series1,bd_series,False,'comA', 'darkorange', r'$chainA$',ax1,'COM-z')
# plot_d1d2(series1,bd_series,False,'comB', 'forestgreen', r'$chainB$',ax1, 'COM-z')
# ax1[-1].legend(fontsize=25)
# fig.savefig("COM_FINAL.svg",dpi=600,format="svg",transparent=True)
# plt.show()




"""
Plot RMSD
"""
# fig,ax = plt.subplots(1,6,sharey=True,figsize=fig_size)

# plot_d1d2(filter_data,bd_series,False,'RMSD-B', 'darkcyan', r'$RMSD$',ax,'RMSD',single_plot=True)


# # ax.plot(bd_data['Timestep']/1e6,bd_data['d2'],color='black',linewidth=2.0,linestyle='',marker='*',markersize=10)
# # ax.legend(fontsize=30)

# for ax1 in ax:
#     ax1.axhline(y=2.0,linestyle='-.',linewidth=2.0,color='black')

# # plt.tight_layout(w_pad=0.05)
# plt.subplots_adjust(top=top,bottom=bottom,left=left,right=right,hspace=hspace,wspace=wspace)
# figlbl = "RMSD_FINAL_" + str(disp_ts) + "ns.svg"
# plt.savefig(figlbl,dpi=600,format="svg",transparent=True)
# plt.show()


"""
Plot RMSD
"""
fig,ax = plt.subplots(1,6,sharey=True,figsize=fig_size)

plot_d1d2(filter_data,bd_series,False,'comB', 'sienna', r'$COM_{z}$',ax,r'$COM_{z}$',single_plot=True)


# ax.plot(bd_data['Timestep']/1e6,bd_data['d2'],color='black',linewidth=2.0,linestyle='',marker='*',markersize=10)
# ax.legend(fontsize=30)

# for ax1 in ax:
#     ax1.axhline(y=2.0,linestyle='-.',linewidth=2.0,color='black')
for ax1 in ax:
    ax1.set_ylim(0,36.0)
    ax1.set_yticks([0,9,18,27,36])
# plt.tight_layout(w_pad=0.05)
plt.subplots_adjust(top=top,bottom=bottom,left=left,right=right,hspace=hspace,wspace=wspace)
figlbl = "COM_FINAL_" + str(disp_ts)
plt.savefig(figlbl + "ns.svg",dpi=600,format="svg",transparent=True)
plt.savefig(figlbl+"ns.png",dpi=600,format="png",transparent=True)
plt.show()
