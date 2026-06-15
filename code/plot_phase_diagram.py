import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import patches
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
from reader import readQtensor

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


def plot_phase_space(dir):
    

    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(7,4),constrained_layout=True)
    
    gamma0s = np.arange(1.6,2.600,0.1)
    phi0s = np.arange(0.6,2.025,0.1)
    avg_q = np.zeros((len(phi0s),len(gamma0s)),dtype=float)
    binder = np.zeros((len(phi0s),len(gamma0s)),dtype=float)

    for i,phi0 in enumerate(phi0s):
        for j,gamma0 in enumerate(gamma0s):

            filename = dir+"/Qtensor_phi0_{:.3f}_gamma0_{:.3f}.txt".format(phi0, gamma0)
            try:
                v, d, phi, d_avg, b = readQtensor(Lx,Ly,Lz,filename)
                avg_q[i][j] = d_avg
                binder[i][j] = b

            except Exception as e:
                print(f"An error occurred: {e}")
                continue

    axs[0].set_xticks(np.arange(1.6, 2.6, 0.2))     
    axs[0].set_yticks(np.arange(0.8, 2.2, 0.2))             
    im = axs[0].imshow(avg_q,origin="lower",cmap=my_cmap,aspect="equal",interpolation="gaussian",extent=[gamma0s[0], gamma0s[-1], phi0s[0], phi0s[-1]])
    axs[0].set_xlabel(r"bare coupling coefficient, $\gamma_0$")
    axs[0].set_ylabel(r"nematic composition, $\phi_0$ ")
    axs[0].set_title(r"average nematic order $\langle q\rangle$")
    # cbar1 = fig.colorbar(im, ax=axs[0], orientation='vertical', fraction=0.03, pad=0.04)

    im = axs[1].imshow(binder,origin="lower",cmap=my_cmap,aspect="equal",interpolation="spline16",extent=[gamma0s[0], gamma0s[-1], phi0s[0], phi0s[-1]])
    axs[1].set_xticks(np.arange(1.6, 2.6, 0.2))
    axs[1].set_yticks(np.arange(0.8, 2.2, 0.2))
    axs[1].set_xlabel(r"bare coupling coefficient, $\gamma_0$")
    axs[1].set_ylabel(r"$\phi_0$ nematic composition")
    
    axs[1].set_title(r"binder cumulant $U_\phi (\kappa=0,W=0)$")
    cbar1 = fig.colorbar(im, ax=axs[1], orientation='vertical', fraction=0.03, pad=0.04)
    
    plt.savefig(dir+"/phase_space_1.pdf")
    plt.savefig(dir+"/phase_space_.png",dpi=300)

    plt.show()

    return

def plot_phase_space_inv(dir): # Added defaults and directory arg

    gamma0s = np.arange(1.6,2.600,0.1)#0.025)
    phi0s = np.arange(0.6,2.025,0.1)#0.025)
    avg_q = np.zeros((len(phi0s),len(gamma0s)),dtype=float)
    binder = np.zeros((len(phi0s),len(gamma0s)),dtype=float)

    # Read data
    for i,phi0 in enumerate(phi0s):
        for j,gamma0 in enumerate(gamma0s):
            filename = dir+"/Qtensor_phi0_{:.3f}_gamma0_{:.3f}.txt".format(phi0, gamma0)
            try:
                output = readQtensor(Lx,Ly,Lz,filename)
                avg_q[i][j] = output[3]
                binder[i][j] = output[4]

            except Exception as e:
                print(f"An error occurred: {e}")
                continue
    
    phi_min, phi_max = np.min(phi0s), np.max(phi0s)
    gamma_min, gamma_max = np.min(gamma0s), np.max(gamma0s)
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(6,4), gridspec_kw={'wspace': 0})
    
    # --- Plot 1: Average Nematic Order ---
    data_to_plot_0 = np.transpose(avg_q)
    im0 = axs[0].imshow(data_to_plot_0[::-1],
                       origin="upper",
                       cmap=my_cmap,
                       aspect="auto", 
                       interpolation="gaussian",
                       vmin=0,
                       extent=[phi_min, phi_max, gamma_min, gamma_max]) 

    axs[0].set_xticks(np.arange(0.6, 2.2, 0.2)) 
    axs[0].set_yticks(np.arange(1.6, 2.6, 0.2)) 
    axs[0].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    axs[0].xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    axs[0].set_xlabel(r"nematic composition, $\phi_0$ ")
    axs[0].set_ylabel(r"bare coupling coefficient, $\gamma_0$")
    axs[0].invert_yaxis()
    axs[0].set_title(r"average nematic order, $\langle q\rangle$")
    
    cax = fig.add_axes([0.2,0.25,0.25,0.75])
    cbar = fig.colorbar(im0, ax=cax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.ax.tick_params(length=0)
    print(np.max(avg_q))
    ticks = np.linspace(0,np.max(avg_q),5)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([str(round(x,2)) for x in ticks])
    cax.axis('off')

    # --- Plot 2: Binder Cumulant ---
    data_to_plot_1 = np.transpose(binder) 
    im1 = axs[1].imshow(data_to_plot_1[::-1],
                       origin="upper",
                       cmap=my_cmap,
                       aspect="auto", 
                       interpolation="gaussian", 
                       vmin=0,
                       extent=[phi_min, phi_max, gamma_min, gamma_max]) 

    axs[1].set_xticks(np.arange(0.6, 2.2, 0.2)) # Ticks for phi on x-axis
    axs[1].set_yticks(np.arange(1.6, 2.6, 0.2)) # Ticks for gamma on y-axis
    axs[1].yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    axs[1].xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    axs[1].set_xlabel(r"nematic composition, $\phi_0$ ")
    # axs[1].set_ylabel(r"bare coupling coefficient, $\gamma_0$") 
    axs[1].invert_yaxis()
    axs[1].set_title(r"binder cumulant, $U_\phi$")
    
    cax = fig.add_axes([0.5,0.25,0.25,0.75])
    cbar = fig.colorbar(im1, ax=cax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.ax.tick_params(length=0)
    ticks = np.linspace(0,np.max(binder),5)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([str(round(x,2)) for x in ticks])
    cax.axis('off')

    #plt.savefig(dir+"/phase_space_inv.pdf")
    #plt.savefig(dir+"/phase_space_inv.png",dpi=300)

    plt.show()

    return

def experimental_phase_space():
    
    # Isotropic (blue) data points
    isotropic_data = []
    isotropic_data.extend([(26, y) for y in [20, 22, 24, 26, 28, 30, 32, 34, 35, 36, 37, 38, 39, 45]])
    isotropic_data.extend([(27, y) for y in [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]])
    isotropic_data.extend([(28, y) for y in [37, 38, 39, 45]])
    isotropic_data.extend([(29, y) for y in [43, 44, 45]])

    # Isotropic/nematic (orange) data points
    isotropic_nematic_data = []
    isotropic_nematic_data.extend([(27, y) for y in [20, 21, 22, 23, 24, 25, 26, 27]])
    isotropic_nematic_data.extend([(28, y) for y in [24, 26, 28, 30, 32, 34, 35, 36]])
    isotropic_nematic_data.extend([(29, y) for y in [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]])
    isotropic_nematic_data.extend([(30, y) for y in [36, 37, 38, 39, 40, 41, 42, 43, 44, 45]])
    isotropic_nematic_data.extend([(31, y) for y in [43, 44, 45]])

    # Nematic (pink) data points
    nematic_data = []
    nematic_data.extend([(28, y) for y in [20, 22]])
    nematic_data.extend([(29, y) for y in [20, 21, 22, 23, 24, 25, 26, 27]])
    nematic_data.extend([(30, y) for y in [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]])
    nematic_data.extend([(31, y) for y in [20, 24, 28, 32, 36, 40, 42]])
    nematic_data.extend([(32, y) for y in [20, 24, 28, 32, 36, 40, 42, 43, 44, 45]])
 
    iso_x, iso_y = zip(*isotropic_data)
    iso_nem_x, iso_nem_y = zip(*isotropic_nematic_data)
    nem_x, nem_y = zip(*nematic_data)

    plt.figure(figsize=(3, 4))
    plt.plot(iso_x, iso_y, 's', color='cornflowerblue', label='Isotropic', markersize=5)
    plt.plot(iso_nem_x, iso_nem_y, 's', color='orange', label='Isotropic/nematic', markersize=5)
    plt.plot(nem_x, nem_y, 's', color='deeppink', label='Nematic', markersize=5)

    ax = plt.gca()
    ax.set_xticks(range(26, 33))
    ax.set_yticks(range(20, 46, 5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))

    concentration = np.arange(25.5,33,0.5)
    ax.fill_between(concentration,concentration*3,concentration*8-188,color="cornflowerblue",alpha=0.3,edgecolor="none")
    ax.fill_between(concentration,concentration*8-188,concentration*20/3-494/3,color="orange",alpha=0.3,edgecolor="none")
    ax.fill_between(concentration,concentration*20/3-494/3,concentration*20/3-494,color="deeppink",alpha=0.3,edgecolor="none")

    plt.grid(True, which='major', linestyle='--', linewidth=0.7, zorder=0)
    plt.grid(True, which='minor', axis='y', linestyle=':', linewidth=0.5, zorder=0)
    plt.xlim(25.5, 32.5)
    plt.ylim(19, 46)
    plt.xlabel("Concentration SSY [wt\%]", fontsize=12)
    plt.ylabel("T [°C]", fontsize=12)
    plt.savefig("/Users/s2862303/PhD/NLBF/paper_draft/figures/experimental_phase_diagram.pdf")
    #plt.legend()
    plt.show()

# phase space --------------------------------------------
Lx = 50
Ly = 50
Lz = 1
directory = "/Users/s2862303/PhD/NLBF/phase_space_k=0_W=0/output"
#plot_phase_space(directory)
plot_phase_space_inv(directory)
#experimental_phase_space()