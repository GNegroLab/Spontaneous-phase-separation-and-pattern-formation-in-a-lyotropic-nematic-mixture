import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import patches
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
from reader import readQtensor
from matplotlib.patches import Circle
from tools import calculate_charge_density, calculate_charge_density_2, smooth_data_log_savgol
from scipy.ndimage import zoom, label
from matplotlib.ticker import ScalarFormatter

# formatting for nice colors
my_cmap=sns.color_palette("rocket", as_cmap=True)
my_cmap_inv=sns.color_palette("rocket_r", as_cmap=True)
my_blue_cmap=sns.color_palette("mako", as_cmap=True)
my_blue_cmap_inv=sns.color_palette("mako_r", as_cmap=True)
my_cmap.set_bad(color='lightgray')   # NaNs appear gray

my_purple = "#830a83"
my_brown_2 = "#9d6f3a"
my_brown = "#bb9d7b"
my_blue = "#3491d3"

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
plt.rcParams["axes.linewidth"]='1.5'
#plt.rcParams["axes.titlepad"]='20'
plt.rcParams["axes.labelpad"]='5'
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.major.width']='1.5'
plt.rcParams['xtick.major.width']='1.5'
plt.rcParams['ytick.major.size']='4.0'
plt.rcParams['xtick.major.size']='4.0'
plt.rcParams['ytick.minor.width']='1.'
plt.rcParams['xtick.minor.width']='1.'
plt.rcParams['ytick.minor.size']='3.0'
plt.rcParams['xtick.minor.size']='3.0'

main_dir="/Users/s2862303/PhD/NLBF/domain_growth_fig3/256x256"

# prepocessing compute number of defects and average distance between defects
n_seeds = 3
Nmax = 100000
step = 1000

Lx, Ly = 256, 256
Lz = 1

time = np.arange(1000,1076000,1000)

n_defects = np.zeros((2,n_seeds,len(time)))
charge_density = np.zeros((2,n_seeds,len(time)))

for i in range(n_seeds):
    seed_dir = main_dir+"/seed{}".format(i+1)
    
    # for j,t in enumerate(time):
        
    #     try:
    #         # no surface tension
    #         filename = seed_dir+"/Qtensor_phi0_1.00_gamma0_2.00_k_0.00_W_0.00_{}.txt".format(t)
    #         v, d, phi, d_avg, b, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)
            
    #         # number of defects
    #         fine_director = zoom(v, (1, 2, 2), order=1) # interpolating the director field to a finer mesh
    #         q1 = calculate_charge_density(v)
        
    #         defect_mask = np.abs(q1) > 0.4    
    #         labeled_array, num_defects = label(defect_mask, structure=np.ones((3, 3), dtype=int) ) # count in clusters for defects between gridpoints
            
    #         one_half_id = np.where(q1>0.4)
    #         minus_one_half_id = np.where(q1<-0.4)
    #         n_defects[0][i][j] = np.sqrt((len(one_half_id[0])+len(one_half_id[1]))/(Lx*Ly))
            
    #         # local charge density
    #         q2 = calculate_charge_density_2(Qxx,Qxy)
    #         charge_density[0][i][j] = np.sum(np.abs(q2))/(Lx*Ly)
            
    #         # surface tension
    #         filename = seed_dir+"/Qtensor_phi0_1.00_gamma0_2.00_k_0.01_W_0.00_{}.txt".format(t)
    #         v, d, phi, d_avg, b, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)
            
    #         # number of defects
    #         q1 = calculate_charge_density(v)
    #         one_half_id = np.where(q1>0.4)
    #         minus_one_half_id = np.where(q1<-0.4)
    #         n_defects[1][i][j] = np.sqrt((len(one_half_id[0])+len(one_half_id[1]))/(Lx*Ly))
            
    #         # local charge density
    #         q2 = calculate_charge_density_2(Qxx,Qxy)
    #         charge_density[1][i][j] = np.sum(np.abs(q2))/(Lx*Ly)
            
    #     except Exception as e:
    #         print(f"Skipping {filename}: {e}")
    #         continue
        
    # np.savetxt("../domain_growth_fig3/charge_data{}.txt".format(i+1),charge_density[:,i,:])
    # np.savetxt("../domain_growth_fig3/n_defects{}.txt".format(i+1),n_defects[:,i,:])
    charge_density[:,i,:] = np.loadtxt("../domain_growth_fig3/charge_data{}.txt".format(i+1))
    n_defects[:,i,:] = np.loadtxt("../domain_growth_fig3/n_defects{}.txt".format(i+1))

# dir = "/Users/s2862303/PhD/NLBF/defect_distance"
# n_defects_st = np.loadtxt(dir+"/kphi=0.01/dd1.dat")
# n_defects_nst = np.loadtxt(dir+"/kphi=0/dd1.dat")
        
# fig, axs = plt.subplots(2,1,figsize=(4,6),sharex=True)
        
# for i in range(n_seeds):
    
    
#     # no surface tension
#     mask = n_defects[0,i,:]!=0.
#     time_m = time[mask]
#     defect_distance = 1/n_defects[0,i,mask]
#     charge_density_m = charge_density[0,i,mask]
    
#     x, y = smooth_data_log_savgol(time_m, defect_distance,window_length=10,polyorder=1,noise_par=0.4)
#     axs[0].plot(x, y,color=my_brown_2,alpha=0.5)
#     axs[1].plot(time_m, charge_density_m,color=my_brown_2,alpha=0.5)
    
#     # surface tension
#     defect_distance = 1/n_defects[1,i,mask]
#     charge_density_m = charge_density[1,i,mask]
#     x, y = smooth_data_log_savgol(time_m, defect_distance,window_length=10,polyorder=1,noise_par=0.4)
#     axs[0].plot(x, y,color=my_purple,alpha=0.5)
#     axs[1].plot(time_m, charge_density_m,color=my_purple,alpha=0.5)

# axs[0].set_xscale("log")
# axs[0].set_yscale("log")
# axs[1].set_xscale("log")
# axs[1].set_yscale("log")
    
# # for the ticks
# formatter = ScalarFormatter()
# formatter.set_scientific(False)
# axs[0].yaxis.set_major_formatter(formatter)
# axs[0].yaxis.set_minor_formatter(formatter)    

# axs[0].set_ylabel(r"average defect distance, $\xi_q$")
# axs[1].set_ylabel(r"average charge density, $\langle q_t \rangle$")
# axs[1].set_xlabel(r"time")
# plt.subplots_adjust()
# plt.tight_layout()
# plt.show()
###

# remove the missing entries
mask = np.all(n_defects != 0, axis=(0, 1))
time = time[mask]
n_defects = n_defects[:,:,mask]
charge_density = charge_density[:,:,mask]

mean_def_noST  = np.mean(1/n_defects[0, :, :], axis=0)
std_def_noST   = np.std(1/n_defects[0, :, :], axis=0)

mean_def_ST    = np.mean(1/n_defects[1, :, :], axis=0)
std_def_ST     = np.std(1/n_defects[1, :, :], axis=0)

mean_charge_noST = np.mean(charge_density[0, :, :], axis=0)
std_charge_noST  = np.std(charge_density[0, :, :], axis=0)

mean_charge_ST   = np.mean(charge_density[1, :, :], axis=0)
std_charge_ST    = np.std(charge_density[1, :, :], axis=0)

fig, axs = plt.subplots(2,1,figsize=(4,6),sharex=True)

axs[0].plot(time, mean_def_noST, color=my_brown_2)
axs[0].fill_between(time, mean_def_noST - std_def_noST,
                           mean_def_noST + std_def_noST,
                           color=my_brown_2, alpha=0.2)

axs[0].plot(time, mean_def_ST, color=my_purple)
axs[0].fill_between(time, mean_def_ST - std_def_ST,
                           mean_def_ST + std_def_ST,
                           color=my_purple, alpha=0.2)


axs[1].plot(time, mean_charge_noST, color=my_brown_2)
axs[1].fill_between(time, mean_charge_noST - std_charge_noST,
                           mean_charge_noST + std_charge_noST,
                           color=my_brown_2, alpha=0.2)

axs[1].plot(time, mean_charge_ST, color=my_purple)
axs[1].fill_between(time, mean_charge_ST - std_charge_ST,
                           mean_charge_ST + std_charge_ST,
                           color=my_purple, alpha=0.2)

# Data for the largest system
trend = np.arange(np.min(time),4000,10)
# exp = 0.35
# axs[0].plot(trend,trend**exp,color="black",linestyle="dashed")
ex1
axs[0].plot(trend,0.4*trend**exp,color="black",linestyle="dashed",label=r"$t^{1/2}$")
# axs[0].plot(n_defects_st[:-1,0],np.sqrt(1/n_defects_st[:-1,1]),color="cyan",label=r"w surface tension $\kappa=0$")
# axs[0].plot(n_defects_nst[:-1,0],np.sqrt(1/n_defects_nst[:-1,1]),color="purple",label=r"wo surface tension $\kappa\neq 0$")

axs[0].set_xscale("log")
axs[0].set_yscale("log")
axs[1].set_xscale("log")
axs[1].set_yscale("log")

# for the ticks
formatter = ScalarFormatter()
formatter.set_scientific(False)
axs[0].yaxis.set_major_formatter(formatter)
axs[0].yaxis.set_minor_formatter(formatter)

axs[0].set_ylabel(r"average defect distance, $\xi_q$")
axs[1].set_ylabel(r"average charge density, $\langle q_t \rangle$")
axs[1].set_xlabel(r"time")
axs[0].legend()
plt.tight_layout()
plt.show()
            
        
        
