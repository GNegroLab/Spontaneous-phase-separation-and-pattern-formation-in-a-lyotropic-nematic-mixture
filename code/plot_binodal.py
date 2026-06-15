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
# 

def binodal_curve(phi0,gamma0s,directory_path):

    fig,ax = plt.subplots(1,1,figsize=(10,8))
    violin_data = []
    minmax_data = np.zeros((len(gamma0s),2))
    std_data = np.zeros((len(gamma0s),2))

    for i,gamma0 in enumerate(gamma0s):

        try:
            file = directory_path+"/Qtensor_phi0_{:.3f}_gamma0_{:.3f}_k_0.000_W_0.000.txt".format(phi0,gamma0)
        except Exception as e:
            print("Exception raised: ", e)
            return

        # read
        v, d, phi, d_avg, binder = readQtensor(Lx,Ly,Lz,file)
        phi = phi.flatten()

        # violin plot
        violin_data.append(phi)

        # identify max and min of domains
        tol = 1e-3
        maxi = np.max(phi)
        mini = np.min(phi)

        if abs(maxi-mini)<tol:
            minmax_data[i] = np.array([maxi,maxi])
            std_data[i] = np.array([np.std(phi),np.std(phi)])

        else:
            maxi_domain = phi[phi>(mini+(maxi-mini)/2)]
            mini_domain = phi[phi<(mini+(maxi-mini)/2)]
            
            minmax_data[i] = np.array([np.average(maxi_domain),np.average(mini_domain)])
            std_data[i] = np.array([np.std(maxi_domain),np.std(mini_domain)])

    plt.violinplot(violin_data,positions=gamma0s,orientation="horizontal")
    plt.xlabel(r"local nematic density, $\phi$")
    plt.ylabel(r"bare coupling coeficient, $\gamma_0$")
    plt.title(r"binodal curve for global density $\phi_0={:.3f}$".format(phi0))
    plt.savefig(directory_path+"/binodal.pdf")
    plt.grid()
    #plt.xticks(np.arange(0,2,0.2))
    plt.show()

    plt.figure(figsize=(4,4))
    plt.errorbar(minmax_data[:,0],gamma0s,xerr=std_data[:,0],marker="o",color=my_brown)
    plt.errorbar(minmax_data[:,1],gamma0s,xerr=std_data[:,1],marker="o",color=my_brown)
    plt.xlabel(r"local nematic density, $\phi$")
    plt.ylabel(r"bare coupling coefficient, $\gamma_0$")
    plt.title(r"binodal curve for global density $\phi_0={:.3f}$".format(phi0))
    plt.savefig(directory_path+"/binodal_2.pdf")
    plt.tight_layout()
    plt.show()

    return

# ------------ binodal curves / violin plot ------------ # 
#directory_path = '/Users/s2862303/PhD/NLBF/binodal_curves/phi0_0.75'
#gamma0s = [round(x,3) for x in np.arange(1.0,2.65,0.05)]
#binodal_curve(0.75,gamma0s,directory_path)
