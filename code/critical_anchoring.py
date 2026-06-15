import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
from matplotlib import patches
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
from reader import readQtensor
from matplotlib.patches import Circle
from tools import detect_plateau, structure_factor, correlation_plot, calculate_charge_density, correlation_oriented, peak_k, linear_fit, linear_fit_2
import os
import glob
import scipy.interpolate as spint
from scipy.signal import find_peaks

# formatting for nice colors
my_cmap=sns.color_palette("rocket", as_cmap=True)
my_cmap_inv=sns.color_palette("rocket_r", as_cmap=True)
my_blue_cmap=sns.color_palette("mako", as_cmap=True)
my_blue_cmap_inv=sns.color_palette("mako_r", as_cmap=True)
my_cmap.set_bad(color='lightgray')   # NaNs appear gray

my_brown = "#bb9d7b"

# for latex font
plt.rcParams.update({
     "text.usetex": True,
     "font.family": "Computer Modern Serif"
})

# makes the font recognizable for inkscape:
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams["axes.grid"] = "False"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["axes.linewidth"]='1.25'
#plt.rcParams["axes.titlepad"]='20'
plt.rcParams["axes.labelpad"]='5'
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['ytick.major.width']='1.25'
plt.rcParams['xtick.major.width']='1.25'

# ------------ critical Ws ------------ # 

Lx=50
Ly=Lx
Lz=1

maindir = "/Users/s2862303/PhD/NLBF/Wc"

Ws_neg=np.arange(-0.045,0,0.001)
Ws_pos=np.arange(0.02,0.086, 0.001)     

kphis = np.arange(0.001,0.015,0.001)
Ls = np.arange(0.005,0.035,0.001) 

# Wcs_normal = np.full((len(kphis),len(Ls)),np.nan)
# Wcs_planar = np.full((len(kphis),len(Ls)),np.nan)

# for i,kphi in enumerate(kphis):
#     for j,L1 in enumerate(Ls):
#         for W in Ws_neg[::-1]:
#             try:            
#                 filename=maindir+"/kphi_{:.3f}_L1_{:.3f}".format(kphi,L1)+"/Qtensor_phi0_1.000_gamma0_2.000_k_{:.3f}_W_{:.3f}.txt".format(kphi,W)
#                 v, d, phi, d_avg, b, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)
                
#                 Sk, magnitude = peak_k(phi,Lx,Ly)           
#                 # if there is a clear laminar structure
#                 if (magnitude>1e4): # and 2*np.pi/Sk < 10 and np.max(phi)<3):
#                     Wcs_normal[i][j] = W        
                                       
                      # plotting individual figures
#                     # fig = plt.figure(figsize=(4,4))
#                     # plt.imshow(np.transpose(phi),cmap=my_cmap,interpolation="gaussian",origin="lower")
#                     # plt.quiver(np.transpose(v[0,:,:]*3*d),
#                     #        np.transpose(v[1,:,:]*3*d),
#                     #        np.transpose(d),
#                     #        linewidth=4,
#                     #        width=.006,
#                     #        headlength=0,       # No head
#                     #        headaxislength=0,   # No axis part of head
#                     #        headwidth=0,        # No head width
#                     #        minlength=0,        # Draw even short arrows
#                     #        pivot='middle',     # Centered arrows (optional)
#                     #        cmap=my_cmap_inv,
#                     #        scale=100,
#                     #        clim=(0,1))
#                     # # print(np.max(phi))
#                     # plt.savefig(maindir+"/kphi_{:.3f}_L1_{:.3f}_N.pdf".format(kphi,L1))
                    
#                     break
                
#             except Exception as e:
#                 print(f"Skipping {filename}: {e}")
#                 continue
            
# print("----- \\ -----")

# for i,kphi in enumerate(kphis):
#     for j,L1 in enumerate(Ls): 
#         for W in Ws_pos:
#             try:    
#                 filename=maindir+"/kphi_{:.3f}_L1_{:.3f}".format(kphi,L1)+"/Qtensor_phi0_1.000_gamma0_2.000_k_{:.3f}_W_{:.3f}.txt".format(kphi,W)
#                 v, d, phi, d_avg, b, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)
                
#                 Sk, magnitude = peak_k(phi,Lx,Ly)    
                
#                 if (magnitude>1e4): # and 2*np.pi/Sk < 10 and np.max(phi)<3):
#                     Wcs_planar[i][j] = W  
                      
                      # plotting individual figures
#                     # fig = plt.figure(figsize=(4,4))
#                     # plt.imshow(np.transpose(phi),cmap=my_blue_cmap,interpolation="gaussian",origin="lower")
#                     # plt.quiver(np.transpose(v[0,:,:]*3*d),
#                     #        np.transpose(v[1,:,:]*3*d),
#                     #        np.transpose(d),
#                     #        linewidth=4,
#                     #        width=.006,
#                     #        headlength=0,       # No head
#                     #        headaxislength=0,   # No axis part of head
#                     #        headwidth=0,        # No head width
#                     #        minlength=0,        # Draw even short arrows
#                     #        pivot='middle',     # Centered arrows (optional)
#                     #        cmap=my_blue_cmap_inv,
#                     #        scale=100,
#                     #        clim=(0,1))
#                     # # plt.savefig(maindir+"/kphi_{:.3f}_L1_{:.3f}_P.pdf".format(kphi,L1))
#                     # plt.show()
                    
#                     break                  
                
#             except Exception as e:
#                 #print(f"Skipping {filename}: {e}")
#                 continue

# # save data
# np.savetxt("Wcs_normal_total.txt", Wcs_normal)
# np.savetxt("Wcs_planar_total.txt", Wcs_planar)

# load data
Wcs_normal = np.loadtxt("../Wc/Wcs_normal_total.txt")
Wcs_planar = np.loadtxt("../Wc/Wcs_planar_total.txt")


plt.rcParams['xtick.labelbottom']=True
plt.rcParams['ytick.labelleft']=True

# X, Y = np.meshgrid(kphis,Ls)     
# fig, ax = plt.subplots(1,2,figsize=(6,3))
# c1 = ax[0].pcolormesh(X, Y, np.transpose(Wcs_normal), cmap=my_cmap, shading='auto')
# c2 = ax[1].pcolormesh(X, Y, np.transpose(Wcs_planar), cmap=my_cmap, shading='auto')
# fig.colorbar(c1, ax=ax[0])
# fig.colorbar(c2, ax=ax[1])
# ax[0].set_title(r'Critical homeotropic anchoring, $W^c_\perp$')
# ax[0].set_xlabel(r'surface tension, $k_\phi$')
# ax[0].set_ylabel(r'elasticity, $L_1$')
# ax[1].set_title(r'Critical planar anchoring, $W^c_\parallel$')
# ax[1].set_xlabel(r'surface tension, $k_\phi$')
# ax[1].set_ylabel(r'elasticity, $L_1$')
# plt.tight_layout()
# plt.show()

fig = plt.figure(figsize=(4.5,4))

colors_hot = [my_cmap(i / (len(Ls) - 1)) for i in range(len(Ls))]
colors_cool = [my_blue_cmap(i / (len(Ls) - 1)) for i in range(len(Ls))]
for i,L1 in enumerate(Ls):    
    plt.plot(kphis[1:],np.abs(Wcs_normal[1:,i]),color=colors_hot[i],marker=".")
    plt.plot(kphis[1:],Wcs_planar[1:,i],color=colors_cool[i],marker=".")

plt.plot(kphis[1:], 4*kphis[1:]+0.025, linestyle="dashed",color="black",linewidth=1.0,label=r"$4\kappa$")
plt.plot(kphis[1:], 2*kphis[1:]+0.015, linestyle="dashed",color="black",linewidth=1.0,label=r"$2\kappa$")
plt.text(0.01,0.05,r"$\uparrow K$")
plt.legend()
plt.ylabel(r"critical anchoring, $|W_c|$")
plt.xlabel(r"surface tension coefficient, $\kappa$")
plt.savefig(maindir+"/Wc_lines.pdf")
#plt.show()

# # # ---- fitting ------# 
linear_fit_2(kphis, Ls, Wcs_normal,"normal_anchoring",planar=False)
linear_fit_2(kphis, Ls, Wcs_planar,"planar_anchoring",planar=True)