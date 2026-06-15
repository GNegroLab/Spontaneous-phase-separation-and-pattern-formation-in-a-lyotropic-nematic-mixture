import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import patches
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
from reader import readQtensor
from matplotlib.patches import Circle
from tools import calculate_charge_density, calculate_charge_density_2

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
#plt.rcParams['ytick.major.width']='1.25'

#plt.rcParams['xtick.major.width']='1.25'
# plt.rcParams['xtick.major.pad']='20'
# plt.rcParams['ytick.major.pad']='20'
# plt.rcParams['xtick.minor.pad']='25'
# plt.rcParams['ytick.minor.pad']='25'

#plt.rcParams['ytick.major.width']='1.2'
#plt.rcParams['xtick.major.width']='1.2'
plt.rcParams['ytick.major.size']='0'
plt.rcParams['xtick.major.size']='0'
plt.rcParams['xtick.labelbottom']=False
plt.rcParams['ytick.labelleft']=False

# dark mode
# plt.rcParams["figure.facecolor"] = "black"
# plt.rcParams["axes.titlecolor"] = "white"
# plt.rcParams["axes.labelcolor"] = "white"
# #plt.rcParams["axes.facecolor"] = "black"
# plt.rcParams["axes.edgecolor"] = "white"
# plt.rcParams["xtick.color"] = "white"
# plt.rcParams["ytick.color"] = "white"
# #legend.labelcolor : "white"


def sci_tex(x):
    mantissa, exp = "{:e}".format(x).split('e')  # e.g. "1.230000", "+04"
    mantissa = float(mantissa)
    exp = int(exp)                               # convert "+04" → 4
    return rf"${mantissa} \times 10^{{{exp}}}$"

def make_frame_Q(filename,iteration,output,phi0,gamma0,k=0.00,W=0.00,movie=False,defects=False,SHIFTX=0,SHIFTY=0):
     
    width = 5
    height = 4
    gap = 0.05*height
    rel_height = 0.9*height
    fig = plt.figure(figsize=(width,height),constrained_layout=True)
          
    v, d, phi, d_avg, binder, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)

    # shift y
    phi = np.roll(phi, shift=SHIFTY, axis=1)
    v = np.roll(v, shift=SHIFTY, axis=2)
    d = np.roll(d, shift=SHIFTY, axis=1)
    
    
    # shift x
    phi = np.roll(phi, shift=SHIFTX, axis=0)
    v = np.roll(v, shift=SHIFTX, axis=1)
    d = np.roll(d, shift=SHIFTX, axis=0)
    
    # Main plot axes
    ax = fig.add_axes([gap/width, 0.05, rel_height/width, 0.9])  # left, bottom, width, height
    im = ax.imshow(np.transpose(phi[0:lim,0:lim]),interpolation="gaussian",origin="lower",cmap=my_cmap,aspect="equal",vmin=0, vmax=2)
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.text(1.5,1.5,r"$t=$"+sci_tex(iteration),fontsize=15,color="white")

    # Insert axis
    # ax2 = fig.add_axes([0.05, 0.05, rel_height/(2*width), 0.9])  # left, bottom, width, height
    # ax2.set_xticks([])
    # ax2.set_yticks([])
    
    R = lim//2+small_lim//2
    L = lim//2-small_lim//2
    im = ax.imshow(np.transpose(phi)[L:R,L:R],interpolation="gaussian",origin="lower",cmap=my_cmap,aspect="equal",vmin=0, vmax=2)
    qv = ax.quiver(np.transpose(v[0,L:R,L:R]*3*d[L:R,L:R]),
                           np.transpose(v[1,L:R,L:R]*3*d[L:R,L:R]),
                           np.transpose(d[L:R,L:R]),
                           linewidth=1,
                           width=.002,
                           headlength=0,       # No head
                           headaxislength=0,   # No axis part of head
                           headwidth=0,        # No head width
                           minlength=0,        # Draw even short arrows
                           pivot='middle',     # Centered arrows (optional) cmap=my_cmap_inv,
                           scale=100,
                           clim=(0,1))
            
    # cax1 = fig.add_axes([gap/width*2+rel_height/width, 0.05, 0.04, 0.9])  # same height as plot
    # cbar = fig.colorbar(im, cax=cax1, orientation='vertical')
    # cbar.set_label(r"compositional phase, $\phi$", labelpad=7.5, fontsize=20)
    # cbar.ax.xaxis.set_label_position('top')  # move label to top
    # cbar.ax.tick_params(axis='x', bottom=False,labelbottom=True, pad=0.1,labelsize=50)

    if defects:

        charge_density = calculate_charge_density(v)
        one_half_id = np.where(charge_density>0.4)
        minus_one_half_id = np.where(charge_density<-0.4)

        for r, c in zip(one_half_id[0],one_half_id[1]):
            center_x = r+0.5
            center_y = c+0.5

            circle = Circle(
                    (center_x, center_y),fill=False,
                    radius = 2,
                    edgecolor = '#ffff00',
                    linewidth = 2,
                    linestyle = '-',
                    alpha = 1
                    )
            ax.add_patch(circle)

        for r, c in zip(minus_one_half_id[0],minus_one_half_id[1]):
            center_x = r+0.5
            center_y = c+0.5

            circle = Circle(
                    (center_x, center_y),fill=False,
                    radius = 2,
                    edgecolor = '#00ffff', 
                    linewidth = 2,      
                    linestyle = '-',    
                    alpha = 1         
                    )
            ax.add_patch(circle)

    ax.set_xlim([0.5,lim-.5])
    ax.set_ylim([0.5,lim-.5])
    if movie:
        plt.savefig(output+"/phi_0={}_gamma_0={}_k={}_W={}_{}.png".format(phi0,gamma0,k,W,iteration),dpi=300)
    else:
        plt.savefig(output+"/phi_0={}_gamma_0={}_k={}_W={}.pdf".format(phi0,gamma0,k,W),transparent=True)   
        #plt.show()
    return

def make_frame_Qxy(filename,iteration,output,phi0,gamma0,k=0.00,W=0.00,movie=False,defects=False):
     
    width = 5
    height = 4
    gap = 0.05*height
    rel_height = 0.9*height
    fig = plt.figure(figsize=(width,height),constrained_layout=True)
          
    v, d, phi, d_avg, binder, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)
    print(np.max(Qxy))
                        # Main plot axes
    ax = fig.add_axes([gap/width, 0.05, rel_height/width, 0.9])  # left, bottom, width, height
    im = ax.imshow(np.transpose(abs(Qxy[0:lim,0:lim])),interpolation="gaussian",origin="lower",cmap=my_cmap,aspect="equal",vmin=0,vmax=0.35)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(1,1,"Iteration {}".format(iteration),fontsize=15,color="white")
            
    cax1 = fig.add_axes([gap/width*2+rel_height/width, 0.05, 0.04, 0.9])  # same height as plot
    cbar = fig.colorbar(im, cax=cax1, orientation='vertical')
    cbar.set_label(r"$|Q_{xy}|$", labelpad=7.5)
    cbar.ax.xaxis.set_label_position('top')  # move label to top
    cbar.ax.tick_params(axis='x', bottom=False,labelbottom=True, pad=0.1)

    if defects:

        charge_density = calculate_charge_density(v)
        one_half_id = np.where(charge_density>0.4)
        minus_one_half_id = np.where(charge_density<-0.4)

        for r, c in zip(one_half_id[0],one_half_id[1]):
            center_x = r+0.5
            center_y = c+0.5

            circle = Circle(
                    (center_x, center_y),fill=False,
                    radius = 2,
                    edgecolor = '#ffff00',
                    linewidth = 2,
                    linestyle = '-',
                    alpha = 1
                    )
            ax.add_patch(circle)

        for r, c in zip(minus_one_half_id[0],minus_one_half_id[1]):
            center_x = r+0.5
            center_y = c+0.5

            circle = Circle(
                    (center_x, center_y),fill=False,
                    radius = 2,
                    edgecolor = '#00ffff', 
                    linewidth = 2,      
                    linestyle = '-',    
                    alpha = 1         
                    )
            ax.add_patch(circle)

    ax.set_xlim([0.5,lim-.5])
    ax.set_ylim([0.5,lim-.5])
    if movie:
        plt.savefig(output+"/phi_0={}_gamma_0={}_k={}_W={}_{}.pdf".format(phi0,gamma0,k,W,iteration),dpi=300)
    else:
        #plt.savefig(output+"/phi_0={}_gamma_0={}_k={}_W={}_{}.pdf".format(phi0,gamma0,k,W,iteration))   
        plt.show()
    return

def make_frame_charge_density(filename,iteration,output,phi0,gamma0,k=0.00,W=0.00,movie=False,defects=False):
     
    width = 5
    height = 4
    gap = 0.05*height
    rel_height = 0.9*height
    fig = plt.figure(figsize=(width,height),constrained_layout=True)
          
    v, d, phi, d_avg, binder, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)        
    #charge_density = calculate_charge_density(v)
    charge_density = calculate_charge_density_2(Qxx,Qxy)

                        # Main plot axes
    ax = fig.add_axes([gap/width, 0.05, rel_height/width, 0.9])  # left, bottom, width, height
    im = ax.imshow(np.transpose(charge_density),interpolation="gaussian",origin="lower",cmap=my_cmap,aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(1,1,"Iteration {}".format(iteration),fontsize=15,color="white")
            
    cax1 = fig.add_axes([gap/width*2+rel_height/width, 0.05, 0.04, 0.9])  # same height as plot
    cbar = fig.colorbar(im, cax=cax1, orientation='vertical')
    cbar.set_label(r"charge density, $q$", labelpad=7.5)
    cbar.ax.xaxis.set_label_position('top')  # move label to top
    cbar.ax.tick_params(axis='x', bottom=False,labelbottom=True, pad=0.1)

    if defects:

        one_half_id = np.where(charge_density>0.4)
        minus_one_half_id = np.where(charge_density<-0.4)

        for r, c in zip(one_half_id[0],one_half_id[1]):
            center_x = r+0.5
            center_y = c+0.5

            circle = Circle(
                    (center_x, center_y),fill=False,
                    radius = 2,
                    edgecolor = '#ffff00',
                    linewidth = 2,
                    linestyle = '-',
                    alpha = 1
                    )
            ax.add_patch(circle)

        for r, c in zip(minus_one_half_id[0],minus_one_half_id[1]):
            center_x = r+0.5
            center_y = c+0.5

            circle = Circle(
                    (center_x, center_y),fill=False,
                    radius = 2,
                    edgecolor = '#00ffff', 
                    linewidth = 2,      
                    linestyle = '-',    
                    alpha = 1         
                    )
            ax.add_patch(circle)

    ax.set_xlim([0.5,lim-.5])
    ax.set_ylim([0.5,lim-.5])
    if movie:
        plt.savefig(output+"/phi_0={}_gamma_0={}_k={}_W={}_{}.pdf".format(phi0,gamma0,k,W,iteration),dpi=300)
    else:
        #plt.savefig(output+"/phi_0={}_gamma_0={}_k={}_W={}_{}.pdf".format(phi0,gamma0,k,W,iteration))   
        plt.show()
    return

def movie(start,n_frames,input,output,defects,phi0=1.00,gamma0=2.00,k=0.00,W=0.00,step=100):

    input += "/Qtensor_phi0_{:.3f}_gamma0_{:.3f}_k_{:.3f}_W_{:.3f}_".format(phi0,gamma0,k,W)    

    for n in range(start,n_frames,step):
        print(n)
        filename = input+"{}.txt".format(n)

        try:
            make_frame_Q(filename,n,output,phi0,gamma0,k,W,True,defects)
        except Exception as e:
            print(f"An error occurred: {e}")
            continue


# single frame --------------------------------------------
Lz = 1
Lx,Ly = 64, 64
lim = Lx
small_lim = lim

gamma0=2.0
phi0=1.0

# single frame
output_dir="/Users/s2862303/Google Drive/My Drive/PhD/NLBF_main/reply/2D/lamella_phi0_1.0_gamma0_2.0_W_-0.03/kphi=0.01/K_0.01"
input_dir="/Users/s2862303/Google Drive/My Drive/PhD/NLBF_main/reply/2D/lamella_phi0_1.0_gamma0_2.0_W_-0.03/kphi=0.01/K_0.01/Qtensor_phi0_1.000_gamma0_2.000_k_0.010_W_-0.030.txt"
make_frame_Q(input_dir,0,output_dir,1.0,2.0,0.00,0.00,False,False)

#(filename,iteration,output,phi0,gamma0,k=0.00,W=0.00,movie=False,defects=False,SHIFTX=0,SHIFTY=0)

# make_frame_Qxy(input_dir,0,output_dir,1.4,2.0,0.00,0.00,False,False)
# output_dir="/Users/s2862303/PhD/NLBF/paper_draft/fig1"
# input_dir="/Users/s2862303/PhD/NLBF/paper_draft/fig1/256x256/Qtensor_phi0_1.00_gamma0_2.00_k_0.00_W_0.00.txt"
# make_frame_Qxy(input_dir,0,output_dir,1.0,2.0,0.00,0.00,False,False)

# movie ---
# output_dir="/Users/s2862303/PhD/NLBF/dummy.pdf"
# input_dir="/Users/s2862303/PhD/NLBF/domain_growth_fig3/seed1/Qtensor_phi0_1.00_gamma0_2.00_k_0.00_W_0.00_49500.txt"
# frames = [1000, 5000, 10000, 20000, 25000, 30000, 50000]
# for frame in frames:
#     input_dir="/Users/s2862303/PhD/NLBF/domain_growth_fig3/128x128/seed1/Qtensor_phi0_1.00_gamma0_2.00_k_0.00_W_0.00_{}.txt".format(frame)
#     output_dir="/Users/s2862303/PhD/NLBF/domain_growth_fig3/128x128/NST_{}_ND.pdf".format(frame)
#     make_frame_Q(input_dir,frame,output_dir,phi0,gamma0,k=0.01,W=0.00,defects=False,movie=False)#,SHIFTX=10)

# Lx,Ly = 128,128
# lim = 128
# small_lim = lim

# for frame in frames:
#     input_dir="/Users/s2862303/PhD/NLBF/domain_growth_fig3/128x128/seed1/Qtensor_phi0_1.00_gamma0_2.00_k_0.01_W_0.00_{}.txt".format(frame)
#     output_dir="/Users/s2862303/PhD/NLBF/domain_growth_fig3/128x128/ST_{}_ND.pdf".format(frame)
#     make_frame_Q(input_dir,frame,output_dir,phi0,gamma0,k=0.01,W=0.00,defects=False,movie=False)#,SHIFTX=45,SHIFTY=20)

# for x in list:
#     input_dir="/Users/s2862303/PhD/NLBF/paper_draft/spider_defects/seed9585/Qtensor_phi0_1.000_gamma0_2.000_k_0.010_W_-0.030_{}.txt".format(x)
#     output_dir="/Users/s2862303/PhD/NLBF/paper_draft/spider_defects/seed9585"
#     make_frame_Q(input_dir,x,output_dir,phi0,gamma0,k=0.01,W=-0.03,defects=True,movie=False)

#make_frame_charge_density(input_dir,82200,output_dir,phi0,gamma0,k=0.01,W=0.00,defects=False,movie=False)
#make_frame_Q(input_dir,100000,output_dir,phi0,gamma0,k=0.01,W=0.00,defects=False,movie=False)
#output_dir="/Users/s2862303/PhD/NLBF/domain_growth_fig3/128x128/seed1"
#input_dir="/Users/s2862303/PhD/NLBF/paper_draft/fig1/256x256/Qtensor_phi0_1.40_gamma0_2.40_k_0.00_W_0.00.txt"
input_dir="/Users/s2862303/Library/CloudStorage/GoogleDrive-anabpaulino@gmail.com/My Drive/PhD/NLBF_main/fig4_super_smectics/defects_fig4/seed9585"
output_dir=input_dir



#movie(0,10000,input_dir,output_dir,True,phi0=1.0,gamma0=2.0,k=0.01,W=-0.03,step=100)
#movie(10000,100000,input_dir,output_dir,False,phi0=1.0,gamma0=2.0,k=0.01,W=0,step=1000)

# filenames="/Users/s2862303/PhD/NLBF/steady_state_diagnosis/128x128/Qtensor/Qtensor_phi0_1.00_gamma0_2.00_k_0.00_W_0.00_"
# total_charge = []
# total_charge_2 = []

# for i in range(0,50000,1000):
#     print("Iteration: ",i)
#     filename=filenames+str(i)+".txt"
#     v, d, phi, d_avg, binder, Qxx, Qxy = readQtensor(Lx,Ly,Lz,filename)        
#     charge_density = calculate_charge_density(v)
#     charge_density_2 = calculate_charge_density_2(Qxx,Qxy)
#     total_charge.append(np.sum(np.abs(charge_density)))
#     total_charge_2.append(np.sum((charge_density_2)))
    
    
# plt.rcParams['xtick.labelbottom']=True
# plt.rcParams['ytick.labelleft']=True

# fig = plt.figure(figsize=(4,4))
# plt.ylabel(r"q")
# plt.xlabel(r"time")
# plt.plot(np.arange(0,50000,1000),total_charge_2,color=my_brown)
# plt.show()