import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import shendrukGroupFormat as ed
#plt.style.use("shendrukGroupStyle")
import seaborn as sns
import matplotlib
from matplotlib import patches
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
from reader import readQtensor
from matplotlib.patches import Circle
from tools import structure_factor, correlation_plot, peak_k, correlation_oriented
import os
import glob
import scipy.interpolate as spint
from scipy.signal import find_peaks

# formatting for nice colors
my_cmap=sns.color_palette("rocket", as_cmap=True)
my_cmap_inv=sns.color_palette("rocket_r", as_cmap=True)
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
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['ytick.major.width']='1.25'
plt.rcParams['xtick.major.width']='1.25'

# plt.rcParams['xtick.major.pad']='20'
# plt.rcParams['ytick.major.pad']='20'
# plt.rcParams['xtick.minor.pad']='25'
# plt.rcParams['ytick.minor.pad']='25'

#plt.rcParams['ytick.major.width']='1.2'
#plt.rcParams['xtick.major.width']='1.2'
#plt.rcParams['ytick.major.size']='0'
#plt.rcParams['xtick.major.size']='0'
#plt.rcParams['xtick.labelbottom']=False
#plt.rcParams['ytick.labelleft']=False

# -------------- header ------------------ #

Lx, Ly = 256, 256
Lz = 1
filename3="/Users/s2862303/PhD/NLBF/changing_anchoring/256x256_9585/Qtensor_phi0_1.000_gamma0_2.000_k_0.010_W_0.057.txt"
filename1="/Users/s2862303/PhD/NLBF/changing_anchoring/256x256_9585/Qtensor_phi0_1.000_gamma0_2.000_k_0.010_W_-0.030.txt"
filename2="/Users/s2862303/PhD/NLBF/changing_anchoring/256x256_9585/Qtensor_phi0_1.000_gamma0_2.000_k_0.010_W_0.000.txt"


# # ------------ structure factor ------------ # 
# filename="/Users/s2862303/PhD/NLBF/Wc_kphi_L1L/kphi_0.001_L1_0.025/Qtensor_phi0_1.000_gamma0_2.000_k_0.001_W_-0.020.txt"
# Lx, Ly = 50,50

# v, d, phi, d_avg, U, Qxx, Qxy = readQtensor(Lx,Ly,1,filename1)
# kx, ky, C_k = structure_factor(phi,Lx,Ly)

# # find the peak
# max_index = np.unravel_index(np.argmax(C_k), C_k.shape)
# diff = max(np.abs(ky[max_index]-np.pi/4)/np.pi/4,np.abs(kx[max_index]-np.pi/4)/np.pi/4)
# kmod = np.sqrt(np.abs(ky[max_index])**2+np.abs(kx[max_index])**2)

# # plot
# fig = plt.figure(figsize=(5,4))
# ax = fig.add_axes([0.1, 0.2, 0.8, 0.7]) 
# im = ax.pcolormesh(kx,ky,C_k,cmap=my_cmap)

# cax1 = fig.add_axes([0.8, 0.2, 0.04, 0.7])  # same height as plot
# cbar = fig.colorbar(im, cax=cax1, orientation='vertical')
# cbar.set_label(r"structure factor, $S(k)$", labelpad=7.5)
# cbar.ax.tick_params(axis='x',length=0,width=0,direction="out",bottom=False,labelbottom=True, pad=0.1)
# ax.set_aspect('equal', adjustable='box') 
# ax.set_xlabel(r"$k_x$")
# ax.set_ylabel(r"$k_y$")
# ax.set_title(r"Normal Anchoring",fontsize=13)
# # ax.axvline(kx[max_index],color="red",linestyle="dashed",linewidth=0.5)
# # ax.axvline(-kx[max_index],color="red",linestyle="dashed",linewidth=0.5)
# # ax.axhline(-ky[max_index],color="red",linestyle="dashed",linewidth=0.5)
# # ax.axhline(ky[max_index],color="red",linestyle="dashed",linewidth=0.5)
# # ax.text(1.9, 2.7, r"$\approx \pm\pi/4 (\pm {}\%)$".format(np.round(diff*100)), color='red', fontsize=11, ha='center', va='center',
# #         bbox=dict(facecolor='black', alpha=0., boxstyle='round,pad=0.3'))
# #ax.text(-3,-3,r"$\lambda^\ast \approx {:.1f}$".format(2*np.pi/kmod),fontsize=12,color="white")
# #plt.savefig("/Users/s2862303/PhD/NLBF/steady_state_diagnosis/{}x{}/structure_factor.pdf".format(Lx,Ly))
# plt.savefig("/Users/s2862303/PhD/NLBF/SI/correlation_SF/SF_normal_anchoring.png",dpi=300)
# plt.tight_layout()
# plt.show()
 
# ------------ G of r - directional ------------ # 
plt.figure(figsize=(4,4))

v, d, phi, d_avg, binder, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename1)
rs, corr = correlation_oriented(phi, v, planar=False)
plt.plot(rs,corr)

v, d, phi, d_avg, binder, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename3)
rs, corr = correlation_oriented(phi, v, planar=True)
plt.plot(rs,corr)

plt.ylabel(r"$G(r_\perp)$")
plt.xlabel(r"$r$")
#plt.savefig("directional_correlation_no_anchoring.pdf")
plt.show()

plt.figure(figsize=(4,4))
plt.plot(rs,1-corr)
plt.plot(rs,-(1-corr),linestyle="dashed")
plt.yscale("log")
plt.xlim([0,50])
plt.ylabel(r"$|1-G(r_\perp)|$")
plt.xlabel(r"$r$")
plt.savefig("directional_correlation_no_anchoring_semilog.pdf")
plt.show()

# ------------ G of r ------------ # 

def compute_Gor(filename):

    v, d, phi, d_avg, binder, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)
    corr = correlation_plot(phi,Lx,Ly)

    # fig = plt.figure(figsize=(4,4))
    # im = plt.imshow(np.transpose(corr),cmap=my_cmap, interpolation="spline16")
    # cbar = fig.colorbar(im,shrink=0.7)
    # plt.xlabel(r"$x$")
    # plt.ylabel(r"$y$")
    # plt.title(r"G(r) - {}x{}".format(Lx,Ly))
    # plt.savefig("/Users/s2862303/PhD/NLBF/steady_state_diagnosis/{}x{}/G_r.pdf".format(Lx,Ly))
    # plt.savefig("/Users/s2862303/PhD/NLBF/steady_state_diagnosis/{}x{}/G_r.png".format(Lx,Ly),dpi=300)
    # plt.tight_layout()
    # plt.show()

    # interpolating the correlation from the 2D grid
    xs = np.arange(-Lx//2,Lx//2,1)
    ys = np.arange(-Ly//2,Ly//2,1)
    interpolator = spint.RectBivariateSpline(xs, ys, corr)

    return interpolator

def decaying_oscillation(t, A, lam, f, phi, C):
    return A * np.exp(-lam * t) * np.cos(2 * np.pi * f * t + phi) + C

# output="/Users/s2862303/PhD/NLBF/paper_draft/figures"
# lim, small_lim = 256, 256
# phi0=1.0
# gamma0=2.0

# int1 = compute_Gor(filename1)
# int2 = compute_Gor(filename2)
# int3 = compute_Gor(filename3)
# ints = [int1,int2,int3]
# Us=[[0,1],[0,1],[0,1]]

# # Interpolate all points
# fig = plt.figure(figsize=(3,3))
# ax = fig.add_axes([0.15,0.15,0.65,0.75])
# plt.xlabel(r"$r$")
# plt.ylabel(r"$G_{\hat{e}_y}(r)$")
# #plt.title(r"$G_u(r), u=({:.2f},{:.2f})$".format(U[0],U[1]))
# plt.axhline(1,linestyle="dashed",color="black")

# colors=["#bb9d7b", "#ff78e4", "#00c2ab"]
# Ws = [-0.030, 0.00, 0.057]

# for i,inter in enumerate(ints):
    
#     # defining the direction vector
#     ux, uy = Us[i]
#     U = np.array([ux,uy])/np.sqrt(ux**2+uy**2)
    
#     # defining the grid
#     grid_dens = 8
#     points = np.zeros((Lx*grid_dens,2))
#     points[:,0] = U[0]*np.arange(0,Lx,1/grid_dens)
#     points[:,1] = U[1]*np.arange(0,Ly,1/grid_dens)
    
#     values = [inter(x, y)[0, 0] for x, y in points]
#     rs = np.sqrt(points[:,0]**2+points[:,1]**2)


#     ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1)) 
#     ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05)) 
#     ax.xaxis.set_major_locator(ticker.MultipleLocator(30)) 
# #    plt.xaxis.set_minor_locator(ticker.MultipleLocator(0.1)) 
#     ax.plot(rs,values,color=colors[i],label=r"$W={:.3f}$".format(Ws[i]))

#     # # fitting the curves
#     # # --- Fit the model ---
#     # initial_guess = [0.9, 0.1, 1.0, 0.0, 0.0]
#     # popt, pcov = curve_fit(decaying_oscillation, rs, values, p0=initial_guess)
#     # A_fit, lam_fit, f_fit, phi_fit, C_fit = popt

#     # # --- Compute derived quantities ---
#     # T_osc = 1 / f_fit       # oscillation period
#     # tau_env = 1 / lam_fit    # envelope decay time constant

#     # # --- Display results ---
#     # print(f"Fitted parameters:")
#     # print(f"  A = {A_fit:.3f}")
#     # print(f"  δ = {lam_fit:.3f}  -> envelope time constant τ = {tau_env:.3f}")
#     # print(f"  f = {f_fit:.3f} -> oscillation period λ = {T_osc:.3f} -> leg width w = {T_osc/2:.3f}")
#     # print(f"  φ = {phi_fit:.3f}")
#     # print(f"  C = {C_fit:.3f}")

#     #plt.plot(rs, decaying_oscillation(rs, *popt), 'r-', label='Fit')

# plt.legend()
# plt.xlim([0,128])
# #plt.tight_layout()
# # plt.savefig("/Users/s2862303/PhD/NLBF/paper_draft/figures/Gor_256x256.pdf")
# plt.show()

# ---- specific structure factor ------# 

# Lx, Ly = 50, 50
# v, d, phi, d_avg, b = readQtensor(Lx,Ly,Lz,filename)
# print(peak_k_2(phi))