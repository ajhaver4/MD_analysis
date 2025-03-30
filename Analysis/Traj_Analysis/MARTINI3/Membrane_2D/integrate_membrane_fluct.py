import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")
pitch_data = pd.read_csv(sys.argv[2],comment='#',delimiter="\t")
deltaz_data = pd.read_csv(sys.argv[3],comment='#',delimiter="\t")

print("Total Number of Configurations: ",len(data['Timestep']))

#Selecting bound and unbound timesteps
mask_bound = (data['d1'] < 12.0) & (data['d2'] < 12.0)
mask_unbound = (data['d1'] > 12.0) | (data['d2'] > 12.0)
#Also selecting mask to ensure proteins are bound on membrane
# mask_mem_bdA = data['min_dist_A'] < 1.2
# mask_mem_bdB = data['min_dist_B'] < 1.2

mask_mem_bdA = (data['comA'] > 10) & (data['comA'] < 15)
mask_mem_bdB = (data['comB'] > 10) & (data['comB'] < 15)


#Mask for bound and unbound protein states including membrane bound chains only
bound_timesteps= (mask_bound & mask_mem_bdA & mask_mem_bdB)       #Using 0.8 nm as cutoff
unbound_timesteps = (mask_unbound & mask_mem_bdA & mask_mem_bdB)

"""
!!!!!!!! - This line is only for debugging. Remove it for actual itnegration
"""
# unbound_timesteps = (mask_bound & mask_mem_bdA & mask_mem_bdB)    

print("Number of Bound Configurations: ",len(data['Timestep'][bound_timesteps]))
print("Number of Unbound Configurations: ",len(data['Timestep'][unbound_timesteps]))

#Another mask which is only for membrane bound chains. 
memb_bd_mask = (mask_mem_bdA & mask_mem_bdB)


n_frames = len(data['Timestep'])
timesteps_array = pitch_data['Timestep']/1e6


figp,[axp4,axp5] = plt.subplots(2,1)

p4A=axp4.hist(pitch_data['PitchA3'][bound_timesteps],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p4B=axp4.hist(pitch_data['PitchB3'][bound_timesteps],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp4.set_title("Distribution of Pitch in Bound State")
axp4.legend()

p5A=axp5.hist(pitch_data['PitchA3'][unbound_timesteps],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p5B=axp5.hist(pitch_data['PitchB3'][unbound_timesteps],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp5.set_title("Distribution of Pitch in Unbound State")
axp5.legend()

figr,[axr4,axr5] = plt.subplots(2,1)
r4A=axr4.hist(pitch_data['RollA4'][bound_timesteps],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r4B=axr4.hist(pitch_data['RollB4'][bound_timesteps],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr4.set_title("Distribution of Roll in Bound State")
axr4.legend()

r5A=axr5.hist(pitch_data['RollA4'][unbound_timesteps],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r5B=axr5.hist(pitch_data['RollB4'][unbound_timesteps],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr5.set_title("Distribution of Roll in Unbound State")
axr5.legend()


figA,[axA,axA1,axA2] = plt.subplots(3,1,sharex=True)

axA.plot(timesteps_array[memb_bd_mask],pitch_data['PitchA3'][memb_bd_mask],label='Pitch Chain A',linewidth=0.7)
axA1.plot(timesteps_array[memb_bd_mask],pitch_data['PitchB3'][memb_bd_mask],label='Pitch Chain B',linewidth=0.7)
axA2.plot(data['Timestep']/1e6,data['comA'],label='COM of A (z)',linewidth=0.7)
axA2.plot(data['Timestep']/1e6,data['comB'],label='COM of B (z)',linewidth=0.7)

axA.legend()
axA1.legend()
axA.set_xlabel("Time")
axA.set_ylabel("Angle ")
axA1.set_ylabel("Angle ")
axA.set_title("Chain A")
axA1.set_title("Chain B")

figB,[axB,axB1,axB2] = plt.subplots(3,1,sharex=True)

axB.plot(timesteps_array[memb_bd_mask],pitch_data['RollA4'][memb_bd_mask],label="Roll Chain A",linewidth=0.7)
axB1.plot(timesteps_array[memb_bd_mask],pitch_data['RollB4'][memb_bd_mask],label="Roll Chain B",linewidth=0.7)
axB2.plot(data['Timestep'][memb_bd_mask]/1e6,data['comA'][memb_bd_mask],label='COM of A (z)',linewidth=0.7)
axB2.plot(data['Timestep']/1e6,data['comB'],label='COM of B (z)',linewidth=0.7)

axB.legend()
axB1.legend()
axB.set_xlabel("Time")
axB.set_ylabel("Angle ")
axB1.set_ylabel("Angle ")
axB.set_title("Chain A")
axB1.set_title("Chain B")
plt.show()


def integrate_pitch(xvalues,ax,low_lim,upp_lim,plot=False):

    radian_val = np.radians(xvalues)
    low_lim = np.radians(low_lim)
    upp_lim = np.radians(upp_lim)
    xhisto = np.histogram(radian_val,bins=len(radian_val))
    frac = xhisto[0]  #Fraction of max value
    if plot:
        ax.hist(radian_val,bins=500,alpha=0.5,color='gold',density=True)
        ax.hist(radian_val,bins=500,weights=frac*np.sin(radian_val),alpha=0.5,color='crimson',density=True)
    
    
    sin_histo = np.histogram(radian_val,bins=500,weights=frac*np.sin(radian_val))
    dtheta = sin_histo[1][1]-sin_histo[1][0]    #bin size
    mask2 = (sin_histo[1][1:] > low_lim)  & (sin_histo[1][1:] <= upp_lim)     #Integration limits
    integral = np.sum(dtheta*sin_histo[0][mask2])/np.max(sin_histo[0])    # Sum of all values in the integration area = Ni*bin_size
    return(integral)

def intetgrate_roll(xvalues,ax,low_lim,upp_lim):
    radian_val = np.radians(xvalues)
    low_lim = np.radians(low_lim)
    upp_lim = np.radians(upp_lim)
    xhisto = np.histogram(radian_val,bins=500)
    ax.hist(radian_val,bins=500,alpha=0.5,color='crimson')
    dtheta = xhisto[1][1]-xhisto[1][0]    #bin size
    mask2 = (xhisto[1][1:] > low_lim)  & (xhisto[1][1:] <= upp_lim)     #Integration limits
    integral = np.sum(dtheta*xhisto[0][mask2])/np.max(xhisto[0][mask2])    # Sum of all values in the integration area = Ni*bin_size
    return(integral)

#Integrate Pitch for chain A
fig1,[ax1,ax2] = plt.subplots(2,1)


total_len = len(pitch_data['PitchA3'][bound_timesteps])+len(pitch_data['PitchA3'][unbound_timesteps])
pitch_bd_intA = integrate_pitch(pitch_data['PitchA3'][bound_timesteps],ax1,70,108,plot=True)

pitch_unb_intA = integrate_pitch(pitch_data['PitchA3'][unbound_timesteps],ax2,70,108,plot=True)

plt.show()



#Integrate Roll for chain A 
roll_bd_intA = intetgrate_roll(pitch_data['RollA4'][bound_timesteps],axr4,0,180)
roll_unb_intA = intetgrate_roll(pitch_data['RollA4'][unbound_timesteps],axr5,0,180)

#Integrate Pitch for chain B
pitch_bd_intB = integrate_pitch(pitch_data['PitchB3'][bound_timesteps],axp4,0,180)
pitch_unb_intB = integrate_pitch(pitch_data['PitchB3'][unbound_timesteps],axp5,0,180)

#Integrate Roll for chain B
roll_bd_intB = intetgrate_roll(pitch_data['RollB4'][bound_timesteps],axr4,0,180)
roll_unb_intB = intetgrate_roll(pitch_data['RollB4'][unbound_timesteps],axr5,0,180)


# #FItting
# def harm2(theta,k):
#     v = k*(1-np.cos(1*(np.radians(theta))))
#     return(v)

def harm(theta,k,theta0):
    v = 0.5*k*(np.radians(theta-theta0))**2
    return(v)

# # """Chain A Pitch"""
# #Bound
# hist_pitchA = p4A
# #Adding additional mask to remove outliers that are observed when proteins are not bound to membrane
# mask1 = (p4A[1][:-1] < 108) & (p4A[1][:-1] > 70)
# x_bd_data_pA = hist_pitchA[1][:-1][mask1]
# y_bd_data_pA = -2.5*np.log(hist_pitchA[0][mask1])+2.5*np.log(np.max(hist_pitchA[0][mask1]))
# params_pitchA_bd = curve_fit(harm,x_bd_data_pA,y_bd_data_pA,[100,90])
# #Unbound
# hist_pitchA2 = p5A
# mask2 = (p5A[1][:-1] < 108) & (p5A[1][:-1] > 70)
# x_unb_data_pA = hist_pitchA2[1][:-1][mask2]
# y_unb_data_pA = -2.5*np.log(hist_pitchA2[0][mask2])+2.5*np.log(np.max(hist_pitchA2[0][mask2]))
# params_pitchA_unb = curve_fit(harm,x_unb_data_pA,y_unb_data_pA,[100,90])

# # """Chain B Pitch"""
# #Bound
# hist_pitchB = p4B
# maskB = (p4B[1][:-1] < 108) & (p4B[1][:-1] > 70)
# x_bd_data_pB = hist_pitchB[1][:-1][maskB]
# y_bd_data_pB = -2.5*np.log(hist_pitchB[0][maskB]) + 2.5*np.log(np.max(hist_pitchB[0][maskB]))
# params_pitchB_bd = curve_fit(harm,x_bd_data_pB,y_bd_data_pB,[100,90],maxfev=800,bounds=(1.0,2000))
# #Unbound
# hist_pitchB2 = p5B
# maskB2 = (p5B[1][:-1] < 108) & (p5B[1][:-1] > 70)
# x_unb_data_pB = hist_pitchB2[1][:-1][maskB2]
# y_unb_data_pB = -2.5*np.log(hist_pitchB2[0][maskB2]) + 2.5*np.log(np.max(hist_pitchB2[0][maskB2]))
# params_pitchB_unb = curve_fit(harm,x_unb_data_pB,y_unb_data_pB,[100,90],maxfev=800,bounds=(1.0,2000))

# #Chain A Pitch plotting
# fig_fit,[ax_fit,ax_fit2] = plt.subplots(2,1)
# ax_fit.plot(x_bd_data_pA,y_bd_data_pA,label='V(x)=-2.5 * log(P(x))')
# ax_fit.plot(x_bd_data_pA,harm(x_bd_data_pA,*params_pitchA_bd[0]),label='V(x) = k*(1-cos(theta))')
# ax_fit.legend()

# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_pitchA_bd[0][0])
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_pitchA_bd[0][1]))
# # s2 = 'n = {:.2f}'.format(params_02[0][1])
# # s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
# ax_fit.text(np.mean(x_bd_data_pA),10,s1)
# ax_fit.text(np.mean(x_bd_data_pA),8,s2)
# ax_fit.set_title("Bound Configurations")

# ax_fit2.plot(x_unb_data_pA,y_unb_data_pA,label='V(x)=-2.5 * log(P(x))')
# ax_fit2.plot(x_unb_data_pA,harm(x_unb_data_pA,*params_pitchA_unb[0]),label='V(x) = k*(1-cos(theta))')
# ax_fit2.legend()

# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_pitchA_unb[0][0])
# # s2 = 'n = {:.2f}'.format(params_02[0][1])
# # s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_pitchA_unb[0][1]))
# ax_fit2.text(np.mean(x_unb_data_pA),10,s1)
# ax_fit2.text(np.mean(x_unb_data_pA),8,s2)
# ax_fit2.set_title("Unound Configurations")
# fig_fit.suptitle('Chain A Pitch')

# #Chain B Pitch plotting
# fig_fit,[ax_fit,ax_fit2] = plt.subplots(2,1)
# ax_fit.plot(x_bd_data_pB,y_bd_data_pB,label='V(x)=-2.5 * log(P(x))')
# ax_fit.plot(x_bd_data_pB,harm(x_bd_data_pB,*params_pitchB_bd[0]),label='V(x) = k*(1-cos(theta))')
# ax_fit.legend()

# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_pitchB_bd[0][0])
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_pitchB_bd[0][1]))
# # s2 = 'n = {:.2f}'.format(params_02[0][1])
# # s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
# ax_fit.text(np.mean(x_bd_data_pB),10,s1)
# ax_fit.text(np.mean(x_bd_data_pB),8,s2)
# ax_fit.set_title("Bound Configurations")

# ax_fit2.plot(x_unb_data_pB,y_unb_data_pB,label='V(x)=-2.5 * log(P(x))')
# ax_fit2.plot(x_unb_data_pB,harm(x_unb_data_pB,*params_pitchB_unb[0]),label='V(x) = k*(1-cos(theta))')
# ax_fit2.legend()
# ax_fit2.set_title('PitchA3')
# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_pitchB_unb[0][0])
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_pitchB_unb[0][1]))
# # s2 = 'n = {:.2f}'.format(params_02[0][1])
# # s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
# ax_fit2.text(np.mean(x_unb_data_pB),10,s1)
# ax_fit2.text(np.mean(x_unb_data_pB),8,s2)
# ax_fit2.set_title("Unbound Configurations")
# fig_fit.suptitle('Chain B Pitch')



# """Chain A Roll"""
# #Bound States
# mask1 = r4A[1][:-1] < 50
# # mask1 = np.nonzero(r4A[0])
# x_bd_data_rA = r4A[1][:-1][mask1]
# y_bd_data_rA = -2.5*np.log(r4A[0][mask1])+2.5*np.log(np.max(r4A[0][mask1]))
# params_bd_rA = curve_fit(harm,x_bd_data_rA,y_bd_data_rA,[100,10])

# #Unbound States
# mask2 = r5A[1][:-1] < 50
# # mask2 = np.nonzero(r5A[0])
# x_unb_data_rA = r5A[1][:-1][mask2]
# y_unb_data_rA = -2.5*np.log(r5A[0][mask2])+2.5*np.log(np.max(r5A[0][mask2]))
# params_unb_rA = curve_fit(harm,x_unb_data_rA,y_unb_data_rA,[100,20])

# """Chain B Roll"""
# #Bound States
# mask1 = r4B[1][:-1] < 50
# # mask1 = np.nonzero(r4B[0])
# x_bd_data_rB = r4B[1][:-1][mask1]
# y_bd_data_rB = -2.5*np.log(r4B[0][mask1])+2.5*np.log(np.max(r4B[0][mask1]))
# params_bd_rB = curve_fit(harm,x_bd_data_rB,y_bd_data_rB,[100,15])

# #Unbound States
# mask2 = r5B[1][:-1] < 50
# # mask2 = np.nonzero(r5B[0])
# x_unb_data_rB = r5B[1][:-1][mask2]
# y_unb_data_rB = -2.5*np.log(r5B[0][mask2])+2.5*np.log(np.max(r5B[0][mask2]))
# params_unb_rB = curve_fit(harm,x_unb_data_rB,y_unb_data_rB,[100,20])

# #Chain A Roll plotting
# fig_rfit,[ax_rfit,ax_rfit2] = plt.subplots(2,1)
# ax_rfit.plot(x_bd_data_rA,y_bd_data_rA,label='V(x)=-2.5 * log(P(x))')
# ax_rfit.plot(x_bd_data_rA,harm(x_bd_data_rA,*params_bd_rA[0]),label='V(x) = k*(1-cos(theta))')
# ax_rfit.legend()
# ax_rfit.set_title('Bound Configurations')
# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_bd_rA[0][0])
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_bd_rA[0][1]))
# # s2 = 'n = {:.2f}'.format(params_02A[0][1])
# ax_rfit.text(np.mean(x_bd_data_rA),10,s1)
# ax_fit.text(np.mean(x_bd_data_rA),8,s2)

# ax_rfit2.plot(x_unb_data_rA,y_unb_data_rA,label='V(x)=-2.5 * log(P(x))')
# ax_rfit2.plot(x_unb_data_rA,harm(x_unb_data_rA,*params_unb_rA[0]),label='V(x) = k*(1-cos(theta))')
# ax_rfit2.legend()
# ax_rfit2.set_title('Unbound Configurations')
# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_unb_rA[0][0])
# # s2 = 'n = {:.2f}'.format(params_02A[0][1])
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_unb_rA[0][1]))
# ax_rfit2.text(np.mean(x_unb_data_rA),10,s1)
# ax_rfit2.text(np.mean(x_unb_data_rA),8,s2)
# fig_rfit.suptitle('Chain A Roll')

# #Chain B Roll plotting
# fig_rfit,[ax_rfit,ax_rfit2] = plt.subplots(2,1)
# ax_rfit.plot(x_bd_data_rB,y_bd_data_rB,label='V(x)=-2.5 * log(P(x))')
# ax_rfit.plot(x_bd_data_rB,harm(x_bd_data_rB,*params_bd_rB[0]),label='V(x) = k*(1-cos(theta))')
# ax_rfit.legend()
# ax_rfit.set_title('Bound Configurations')
# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_bd_rB[0][0])
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_bd_rB[0][1]))
# # s2 = 'n = {:.2f}'.format(params_02A[0][1])
# ax_rfit.text(np.mean(x_bd_data_rB),10,s1)
# ax_rfit.text(np.mean(x_bd_data_rB),8,s2)

# ax_rfit2.plot(x_unb_data_rB,y_unb_data_rB,label='V(x)=-2.5 * log(P(x))')
# ax_rfit2.plot(x_unb_data_rB,harm(x_unb_data_rB,*params_unb_rB[0]),label='V(x) = k*(1-cos(theta))')
# ax_rfit2.legend()
# ax_rfit2.set_title('Unbound Configurations')
# s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_unb_rB[0][0])
# s2 = r'$\mu = {:.2f}$ rad'.format(np.radians(params_unb_rB[0][1]))
# # s2 = 'n = {:.2f}'.format(params_02A[0][1])
# ax_rfit2.text(np.mean(x_unb_data_rB),10,s1)
# ax_rfit2.text(np.mean(x_unb_data_rB),8,s2)
# fig_rfit.suptitle('Chain B Roll')

# plt.show()


#Printing out the integrals
print("Pitch Integral Bound Chain A: ",pitch_bd_intA)
print("Pitch Integral Unbound Chain A: ",pitch_unb_intA)
print("Roll Integral Bound Chain A: ",roll_bd_intA)
print("Roll Integral Unbound Chain A: ",roll_unb_intA)

print("Pitch Integral Bound Chain B: ",pitch_bd_intB)
print("Pitch Integral Unbound Chain B: ",pitch_unb_intB)
print("Roll Integral Bound Chain B: ",roll_bd_intB)
print("Roll Integral Unbound Chain B: ",roll_unb_intB)

"""
DeltaZ Calculations and Integration
"""
figp,[axp4,axp5] = plt.subplots(2,1)

p4A=axp4.hist(deltaz_data['z_A_membres'][bound_timesteps],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p4B=axp4.hist(deltaz_data['z_B_membres'][bound_timesteps],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp4.set_title("Distribution of deltaZ in Bound State")
axp4.legend()

p5A=axp5.hist(deltaz_data['z_A_membres'][unbound_timesteps],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p5B=axp5.hist(deltaz_data['z_B_membres'][unbound_timesteps],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp5.set_title("Distribution of deltaZ in Unbound State")
axp5.legend()

plt.show()

def integrate_deltaz(xvalues):
    mask = (xvalues >0.5) & (xvalues < 2.0)
    xhisto = np.histogram(xvalues[mask],bins=100)
    dz = xhisto[1][1]-xhisto[1][0]    #bin size
    integral = np.sum(dz*xhisto[0])/np.max(xhisto[0])    # Sum of all values in the integration area = Ni*bin_size
    return(integral)

#Integrate DeltaZ for chain A
z_bd_intA = integrate_deltaz(deltaz_data['z_A_membres'][bound_timesteps])
z_unb_intA = integrate_deltaz(deltaz_data['z_A_membres'][unbound_timesteps])

#Integrate DeltaZ for chain B
z_bd_intB = integrate_deltaz(deltaz_data['z_B_membres'][bound_timesteps])
z_unb_intB = integrate_deltaz(deltaz_data['z_B_membres'][unbound_timesteps])

#Fitting deltaz to harmonic potential
def harm2(x,k,r0):
    c = (2*3.14*2.5/k)**0.5
    v = 0.5*k*(x-r0)**2 + 2.5*np.log(c)
    return(v)

# #Fitting deltaz distributions
# def fit_deltaz(hist_data,low_lim,upp_lim,ax,lbl):
#     mask0 = (hist_data[1][:-1] > low_lim)  & (hist_data[1][:-1] < upp_lim)      #Selecting appropriate range for fitting based on previous histogram
#     mask1 = np.nonzero(hist_data[0][mask0])                          #Sometimes additional filter is needed to remove non-zero values
#     x_bd_data_z = hist_data[1][:-1][mask0][mask1]
#     y_bd_data_z = -2.5*np.log(hist_data[0][mask0][mask1])
#     params_bd_dz = curve_fit(harm2,x_bd_data_z,y_bd_data_z,[12,1.0])

#     ax.plot(x_bd_data_z,y_bd_data_z,label='V(x)=-2.5 * log(P(x))')
#     ax.plot(x_bd_data_z,harm2(x_bd_data_z,*params_bd_dz[0]),label='V(x) = 0.5*k*x^2')
#     s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_bd_dz[0][0])
#     s2 = 'r0 = {:.2f} nm'.format(params_bd_dz[0][1])
#     ax.text(np.mean(hist_data[1]),10,s1)
#     ax.text(np.mean(hist_data[1]),8,s2)
#     ax.legend()
#     ax.set_title(lbl)

#     return(params_bd_dz)

# fig_fit,[ax_fit2,ax_fit3] = plt.subplots(2,1)
# params_bd_dzA = fit_deltaz(p4A,0.5,1.8,ax_fit2,'Bound Configurations')
# params_unb_dzA = fit_deltaz(p5A,0.5,1.8,ax_fit3,'Unbound Configurations')
# fig_fit.suptitle("Chain A Fit DeltaZ")

# fig_fit,[ax_fit2,ax_fit3] = plt.subplots(2,1)
# params_bd_dzB = fit_deltaz(p4B,0.5,1.8,ax_fit2,'Bound Configurations')
# params_unb_dzB = fit_deltaz(p5B,0.5,1.8,ax_fit3,'Unbound Configurations')
# fig_fit.suptitle("Chain B Fit DeltaZ")

# plt.show()



"""Printing out the integrals"""
print("DeltaZ Integral Bound Chain A: ",z_bd_intA)
print("DeltaZ Integral Unbound Chain A: ",z_unb_intA)
print("DeltaZ Integral Bound Chain B: ",z_bd_intB)
print("DeltaZ Integral Unbound Chain B: ",z_unb_intB)


# #Plotting the distributions
# fsize=(6,8)
# t_size=18
# leg_size=15
# f_size=20

# fig,ax = plt.subplots(figsize=fsize)

# ax.hist(deltaz_data['z_A_membres'][bound_timesteps],bins=100,label='Chain A',alpha=0.7,density=True,color='darkturquoise')
# ax.hist(deltaz_data['z_B_membres'][bound_timesteps],bins=100,label='Chain B',alpha=0.5,density=True,color='cyan')
# ax.set_xlabel(r'$dz$ (nm)',fontsize=f_size)
# ax.set_ylabel(r'$P(dz)$',fontsize=f_size)
# fig.tight_layout()
# axA.tick_params(labelsize=t_size)
# axA.legend(fontsize=leg_size)
# # axA.set_title("UnBound State Distribution - Roll (Chain A & Chain B)")

# fig,ax = plt.subplots(figsize=fsize)

# ax.hist(pitch_data['PitchA3'][bound_timesteps],bins=100,label='Chain A',alpha=0.7,density=True,color='limegreen')
# ax.hist(pitch_data['PitchB3'][bound_timesteps],bins=100,label='Chain B',alpha=0.5,density=True,color='seagreen')
# ax.set_xlabel(r'$ Pitch (\thetha)$ (deg)',fontsize=f_size)
# ax.set_ylabel(r'$P(\thetha)$',fontsize=f_size)
# fig.tight_layout()
# axA.tick_params(labelsize=t_size)
# axA.legend(fontsize=leg_size)
# # axA.set_title("UnBound State Distribution - Roll (Chain A & Chain B)")

# fig,ax = plt.subplots(figsize=fsize)

# ax.hist(pitch_data['RollA4'][bound_timesteps],bins=100,label='Chain A',alpha=0.7,density=True,color='darkorange')
# ax.hist(pitch_data['RollB4'][bound_timesteps],bins=100,label='Chain B',alpha=0.5,density=True,color='sandybrown')
# ax.set_xlabel(r'$ Roll (\phi)$ (deg)',fontsize=f_size)
# ax.set_ylabel(r'$P(\phi)$',fontsize=f_size)
# fig.tight_layout()
# axA.tick_params(labelsize=t_size)
# axA.legend(fontsize=leg_size)

# plt.show()
