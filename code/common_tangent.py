import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# matplotlib formatting -----
plt.rcParams.update({
     "text.usetex": True,
     "font.family": "Computer Modern Serif"
})

plt.rcParams["axes.grid"] = "False"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["axes.edgecolor"] = "black"
#plt.rcParams["axes.linewidth"] = "1"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

#plt.rcParams['xtick.major.pad']='20'
#plt.rcParams['ytick.major.pad']='20'
#plt.rcParams['xtick.minor.pad']='25'
#plt.rcParams['ytick.minor.pad']='25'

plt.rcParams["axes.linewidth"]='1.25'
#plt.rcParams["axes.titlepad"]='20'
plt.rcParams["axes.labelpad"]='5'
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['ytick.major.width']='1.25'
plt.rcParams['xtick.major.width']='1.25'
#plt.rcParams['ytick.major.size']='0'
#plt.rcParams['xtick.major.size']='0'
#plt.rcParams['xtick.labelbottom']=True
#plt.rcParams['ytick.labelleft']=True

# dark mode
plt.rcParams["figure.facecolor"] = "black"
plt.rcParams["axes.titlecolor"] = "white"
plt.rcParams["axes.labelcolor"] = "white"
#plt.rcParams["axes.facecolor"] = "black"
plt.rcParams["axes.edgecolor"] = "white"
plt.rcParams["xtick.color"] = "white"
plt.rcParams["ytick.color"] = "white"
#legend.labelcolor : "white"


# ----------- basic functions definitions
def gamma(gamma0,phi):
    return gamma0 + Delta*phi

def q(gamma0,phi):
    Q = np.zeros_like(phi)
    Q[gamma(gamma0,phi)>=2.7] = .25*(1+np.sqrt(9-24/gamma(gamma0,phi)))[gamma(gamma0,phi)>=2.7]
    return Q
    
def f_hom(gamma0,phi):
    gamma_val = gamma(gamma0,phi)
    q_val = q(gamma0,phi)
    return A0/3*(1-gamma_val/3)*(q_val**2) - 2*A0*gamma_val/27*(q_val**3) + A0*gamma_val*q_val**4/9 + a/2*phi**2

# ------ derivatives
def qp(gamma0,phi):
    gamma_val = gamma(gamma0,phi)
    return (3 * Delta / (gamma_val**2) ) / np.sqrt(9-24/gamma_val)

def qpp(gamma0,phi):
    gamma_val = gamma(gamma0,phi)
    return -54 * Delta**2 * (gamma_val-2) / (gamma_val**4 * (9-24/gamma_val)**1.5)
 
def fp(gamma0,phi):

    # Auxiliary variables
    gamma_val = gamma(gamma0,phi)
    q_val = q(gamma0,phi)
    q_prime = qp(gamma0,phi)
    
    # Terms independent of q'
    f = - (A0 * Delta / 9) * (q_val**2)
    f -= (2 * A0 * Delta / 27) * (q_val**3)
    f += (A0 * Delta / 9) * (q_val**4)
    f += a * phi

    # Coefficient of q_val'
    f += ( (2 * A0 / 3) * (1 - gamma_val / 3) * q_val - (2 * A0 * gamma_val / 9) * (q_val**2) + (4 * A0 * gamma_val / 9) * (q_val**3) ) * q_prime    

    if isinstance(phi, (float, np.floating)):
        res = a*phi
        if gamma_val>=2.7:
            res = f
        return res

    res = phi*a
    res[gamma_val>=2.7] = f[gamma_val>=2.7]

    return res       

def fpp(gamma0,phi):

    # Auxiliary variables
    gamma_val = gamma(gamma0,phi)
    q_val = q(gamma0,phi)
    q_prime = qp(gamma0,phi)
    q_double_prime = qpp(gamma0,phi)

    # Constant coefficient
    f = a

    # Coefficients of q'
    f += 2 * ( - (2 * A0 * Delta / 9) * q_val - (2 * A0 * Delta / 9) * (q_val**2) + (4 * A0 * Delta / 9) * (q_val**3)) * q_prime

    # Coefficients of (q')^2
    f += ( (2 * A0 / 3) * (1 - gamma_val / 3) - (4 * A0 * gamma_val / 9) * q_val + (12 * A0 * gamma_val / 9) * (q_val**2)) * (q_prime**2)

    # Coefficient of q''
    f += (  (2 * A0 / 3) * (1 - gamma_val / 3) * q_val - (2 * A0 * gamma_val / 9) * (q_val**2) + (4 * A0 * gamma_val / 9) * (q_val**3) ) * q_double_prime

    if isinstance(phi, (float, np.floating)):
        res = a
        if gamma_val>=2.7:
            res = f
        return res

    res = np.ones_like(phi)*a
    res[gamma_val>=2.7] = f[gamma_val>=2.7]

    return res

# ------- global variables and common tangent computation
Delta = 0.1
A0 = 0.25
a = 0.01

def compute_common_tangent(plot=False,generate=False):

    # local variables
    gamma0s = np.linspace(0.2,2.7,20) # only after 0.7 does the critical point appear (2.7 - delta*maxphi)
    phis = np.arange(-0.5,3.5,0.001)
    spinodal = []
    binodal = []
    binodal_interval = 0
    
    for gamma0 in gamma0s:
    
        # free energy and second derivative
        f_pp = fpp(gamma0,phis)
        f_p = fp(gamma0,phis)
        f = f_hom(gamma0,phis)
        
        if plot:
            fig = plt.figure(figsize=(4,3))
            ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
            ax.plot(phis,f,color="black",label=r"$f_{hom}$")
            ax.plot(phis,f_pp,linestyle="dashed",color="black",label=r"$f_{\phi\phi}$",linewidth=0.9)

        # spinodal
        spinodal_ind = np.argmin(np.abs(f_pp))
    
        # critical point
        critical_ind = np.argmax(np.abs(f_pp))
        f_pp[critical_ind-1] = np.nan
        spinodal.append([phis[spinodal_ind],phis[critical_ind]])
    
        # common tangent
        # greedy solution
        df = fp(gamma0,phis)
        tol=1e-3
        neighb=5
        pairs = []
        for i in range(1,critical_ind-neighb,1):
            for j in range(critical_ind+neighb,len(phis),1):
                phi = phis[i]
                phi2 = phis[j]
                cond2 = abs(f[i] - phi*f_p[i] - ( f[j] - phi2*f_p[j] ))
                cond1 = abs(f_p[i] - f_p[j]) 
                if cond1<tol and cond2<tol:
                    pairs.append((phi,phi2,cond1+cond2))
                
        pairs = sorted(pairs,key=lambda x: x[2])
        if len(pairs):
            m = (f_hom(gamma0,pairs[0][1])-f_hom(gamma0,pairs[0][0]))/(pairs[0][1]-pairs[0][0])
            # equation of a line
            b = -m*pairs[0][0] + f_hom(gamma0,pairs[0][0])


            if plot:
                ax.axvspan(pairs[0][0],pairs[0][1],color="#fadccb",label="spinodal")
                ax.axvspan(phis[spinodal_ind],phis[critical_ind],color="#d49dcb",label="binodal")
                ax.plot([pairs[0][0],pairs[0][1]],[pairs[0][0]*m+b,pairs[0][1]*m+b],color="#d052b6")
                ax.scatter(pairs[0][0],f_hom(gamma0,pairs[0][0]),color="#cf9d79",edgecolor="#cf9d79",zorder=3)
                ax.scatter(pairs[0][1],f_hom(gamma0,pairs[0][1]),color="#cf9d79",edgecolor="#cf9d79",zorder=3)
                ax.scatter(phis[spinodal_ind],f_hom(gamma0,phis[spinodal_ind]),color="#9c512a",edgecolor="#9c512a",zorder=3)
                ax.scatter(phis[critical_ind],f_hom(gamma0,phis[critical_ind]),color="#d052b6",edgecolor="#d052b6",zorder=3)
    
            binodal.append([gamma0,pairs[0][0],pairs[0][1]])
            binodal_interval = pairs[0][1]-pairs[0][0]

    
        # formatting
        if plot:
            ax.axhline(0,color="black",linewidth=0.5)
            ax.set_title(r"$\gamma_0={:.2f}$".format(gamma0))
            ax.set_xlim([0.,2.])
#            ax.set_ylim([min(-.0005,np.min(f)-0.001),a+0.001])
            ax.set_xlabel(r"$\phi$")
            ax.set_ylabel(r"$f_{hom}$")
            ax.legend(edgecolor="none",borderpad=0.5,framealpha=0.4,loc='center left',bbox_to_anchor=(.025, .55))
            plt.savefig("common_tangent_{:.2f}.pdf".format(gamma0))
            plt.show()

    if generate:
        np.savetxt("binodal_X.txt",np.array(binodal))
        np.savetxt("spinodal_X.txt",np.array(spinodal))
    
    return binodal_interval
    
# load data
gamma0s = np.linspace(0.2,2.7,20) # only after 0.7 does the critical point appear (2.7 - delta*maxphi)
phis = np.arange(-0.5,3.5,0.001)

compute_common_tangent(plot=True,generate=True)
spinodal = np.loadtxt("spinodal_X.txt")
binodal = np.loadtxt("binodal_X.txt")

# plot
# fig = plt.figure(figsize=(4,3))
# ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
# ax.plot(gamma0s,spinodal[:,1],color="#bb9d7b",label=r"$\phi_c$",linewidth=2)
# ax.fill_between(binodal[:,0],binodal[:,1],binodal[:,2],color="#f8e0cd",alpha=0.7,edgecolor="none",label='binodal region')
# ax.scatter(2.3,0.07,color="#f8e0cd",edgecolor="#9a92a0")
# ax.scatter(2.3,1.1,color="#f8e0cd",edgecolor="#9a92a0")
# ax.fill_between(gamma0s,spinodal[:,0],spinodal[:,1],color="#f59f76",edgecolor="none",label='spinodal region')
# ax.set_xlim([0.3,2.3])
# ax.set_ylim([0.0,2.])
# ax.set_ylabel(r"$\phi$")
# ax.set_xlabel(r"$\gamma_0$")
# ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5)) 
# ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1)) 
# ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5)) 
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1)) 
# ax.legend(edgecolor="none",borderpad=0.5)
# plt.tight_layout()
# plt.savefig("binodal_tangent.pdf")
# plt.show()

# to invert the graph
fig = plt.figure(figsize=(4,3))
ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
ax.plot(spinodal[:, 1], gamma0s, color="#bb9d7b", label=r"$\phi_c$", linewidth=2) 
ax.fill_betweenx(gamma0s, binodal[:, 1], binodal[:, 2], color="#f8e0cd", alpha=0.7, edgecolor="none", label='binodal region') 
#ax.scatter(0.07, 2.3, color="#f8e0cd", edgecolor="#9a92a0")
#ax.scatter(1.1, 2.3, color="#f8e0cd", edgecolor="#9a92a0") 
ax.fill_betweenx(gamma0s, spinodal[:, 0], spinodal[:, 1], color="#f59f76", edgecolor="none", label='spinodal region') 
ax.set_ylim([0.3, 2.3]) 
ax.set_xlim([0.0, 2.])  
ax.set_xlabel(r"nematic composition, $\phi$") 
ax.set_ylabel(r"bare coupling coefficient, $\gamma_0$") 
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5)) 
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1)) 
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5)) 
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1)) 
ax.invert_yaxis()
ax.legend(edgecolor="none", borderpad=0.5)
plt.savefig("binodal_tangent_swapped_inverted.pdf")
plt.show()

# # ------------ Varying A0/a ------------ # 
# gamma0=2.0
# A0s = np.arange(0.05,0.35,0.025)
# aas = np.arange(0.01,0.02,0.0025)

# calculate = 0
# if calculate:
#     binodal_intervals = np.zeros((len(A0s),len(aas)))
# else:
#     binodal_intervals = np.loadtxt("binodal_intervals.txt")

# if calculate:    
#     for i,A in enumerate(A0s):
#         for j,aa in enumerate(aas):
#             A0=A
#             binodal_intervals[i][j] = compute_common_tangent(plot=False,generate=False)
#     np.savetxt("binodal_intervals.txt",binodal_intervals)
    
# # for i,a in enumerate(aas):
# #     plt.plot(A0s/a, binodal_intervals[:,i],label=r"$a={:.3f}$".format(a),marker="o")

# for i,A in enumerate(A0s):
#     plt.plot(A/aas, binodal_intervals[i,:],label=r"$A_0={:.3f}$".format(A),marker="o")
    
# plt.legend()
# plt.xlabel(r"$A_0/a$")
# plt.ylabel(r"binodal interval")
# plt.savefig("binodal_intervals_2.pdf")
# plt.show()
# # ------------ Varying A0/a ------------ # 