import matplotlib.pyplot as plt
import numpy as np
import reader
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d 
from scipy.signal import savgol_filter


#import shendrukGroupFormat as ed
#plt.style.use("shendrukGroupStyle")

# detect plateau
def detect_plateau(data, bin, tolerance, stepsize=100, plot=False):

    prev = np.average(data[0:bin])
    found = False
    for i in range(1, len(data) // bin):

        if plot:
            plt.axvline(i * bin)

        current = np.average(data[i * bin : (i + 1) * bin])
        var = abs(current - prev) / prev

        if var > tolerance:
            prev = current
        else:
            print(
                "Found plateau at iteration {} with a variation of {}%".format(
                    i * bin * stepsize, round(var * 100, 2)
                )
            )
            found = True
            break

    if plot:
        plt.plot(data)
        plt.show()

    if found:
        return i * bin, data[i * bin]
    else:
        return None
    
# filename = "/Users/s2862303/PhD/NLBF/steady_state_diagnosis/256x256/energy_phi0_1.00_gamma0_2.00_k_0.00_W_0.00.txt"
# data = np.loadtxt(filename)
# detect_plateau(data[:,3],20,0.02,plot=True)

def structure_factor(phi,Lx,Ly):

    phi -= np.average(phi)

    C_k = np.fft.fft2(phi)
    C_k = np.fft.fftshift(C_k)
    C_k = np.real(C_k)**2 + np.imag(C_k)**2

    freqx = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Lx, d=1))  # x frequencies f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
    freqy = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Ly, d=1))  # d is the units of the spacing between gridpoints
                                                              # convert form frequency (f, cycles per unit length) to wave vector k
    kx, ky = np.meshgrid(freqx, freqy)
    
    # low-pass filter to detect the legs patterns
    R = np.sqrt(kx**2 + kx**2)
    mask = R > (2*np.pi/8)
    
    C_k *= mask

    return kx, ky, C_k

def peak_k(phi,Lx,Ly):
    
    kx, ky, C_k = structure_factor(phi,Lx,Ly)

    # find the peak
    magnitude = np.max(C_k)
    max_index = np.unravel_index(np.argmax(C_k), C_k.shape)
    kmod = np.sqrt(np.abs(ky[max_index])**2+np.abs(kx[max_index])**2)
    
    return kmod, magnitude


def correlation_plot(phi,Lx,Ly):

    correlation_map = np.zeros((Lx,Ly),dtype=float)
    mean_phi = np.average(phi)
    
    for i in range(0,Lx,1):
        for j in range(0,Ly,1):
 
            phi_aux = np.roll(phi,shift=i, axis=1)
            phi_shifted = np.roll(phi_aux,shift=j,axis=0)
            correlation_map[i,j] = np.average(phi*phi_shifted)

    correlation_map = np.roll(np.roll(correlation_map,shift=Ly//2, axis=0),shift=Lx//2, axis=1)
    return correlation_map/mean_phi**2

def bilinear_interpolation(grid, x, y):
    
    Lx, Ly = grid.shape 
    
    x1, y1 = int(np.floor(x)) % Lx, int(np.floor(y)) % Ly
    x2, y2 = (x1 + 1) % Lx, (y1 + 1) % Ly
    dx, dy = x - np.floor(x), y - np.floor(y)

    Q11 = grid[x1, y1]
    Q21 = grid[x1, y2]
    Q12 = grid[x2, y1]
    Q22 = grid[x2, y2]

    return Q11 * (1 - dx) * (1 - dy) + Q21 * dx * (1 - dy) + Q12 * (1 - dx) * dy + Q22 * dx * dy


def correlation_oriented(phi,n_director,planar=False,grid_points=600):
    
    Lx, Ly = np.shape(phi)
    
    # dr array
    rs = np.linspace(0,np.sqrt((Lx/2)**2+(Ly/2)**2),grid_points)
    corr = np.zeros_like(rs)
    
    thetas = []
    # compute correlation in perpendicular direction
    for i in range(0,Lx,1):
        for j in range(0,Ly,1):
            
            nx, ny = n_director[0,i,j], n_director[1,i,j]
            if planar:
                nx, ny = -ny, nx # normal vector to the lamellae
            tan_theta = ny/nx
            thetas.append(np.arctan2(ny,nx))
            
            for it, r in enumerate(rs):
                
                rx, ry = r / np.sqrt(1 + tan_theta**2), r * tan_theta / np.sqrt(1 + tan_theta**2)
                corr[it] += phi[i,j]*phi[int(i+rx)%Lx,int(j+ry)%Ly] #bilinear_interpolation(phi, i+rx, j+ry)
    
    #fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    #ax.hist(thetas)
    #plt.show()
    
    return rs, corr/(np.average(phi)**2*Lx*Ly)

def wrap_nematic_angle_difference(dtheta):
    
    dtheta_wrapped = dtheta.copy()
    dtheta_wrapped[dtheta > np.pi / 2] -= np.pi
    dtheta_wrapped[dtheta <= -np.pi / 2] += np.pi
    
    return dtheta_wrapped

def calculate_charge_density(director_field):

    nx, ny = director_field[0,:,:], director_field[1,:,:]
    theta = np.arctan2(ny, nx)

    thetaR = np.roll(theta, shift=-1, axis=0)
    thetaT = np.roll(theta, shift=-1, axis=1)
    thetaTR = np.roll(thetaR, shift=-1, axis=1)

    # counter-clock-wise
    dtheta1 = wrap_nematic_angle_difference(thetaR - theta)
    dtheta2 = wrap_nematic_angle_difference(thetaTR - thetaR)
    dtheta3 = wrap_nematic_angle_difference(thetaT - thetaTR)
    dtheta4 = wrap_nematic_angle_difference(theta - thetaT)
    dtheta = dtheta1 + dtheta2 + dtheta3 + dtheta4

    charge_field = dtheta / (2 * np.pi)

    return charge_field

def calculate_charge_density_2(Qxx,Qxy):

    dxQxx = (np.roll(Qxx, shift=+1, axis=0) - np.roll(Qxx, shift=-1, axis=0))/2
    dyQxy = (np.roll(Qxy, shift=+1, axis=1) - np.roll(Qxy, shift=-1, axis=1))/2
    dyQxx = (np.roll(Qxx, shift=+1, axis=1) - np.roll(Qxx, shift=-1, axis=1))/2
    dxQxy =(np.roll(Qxy, shift=+1, axis=0) - np.roll(Qxy, shift=-1, axis=0))/2
    
    charge_field = (dxQxx*dyQxy - dxQxy*dyQxx) / (2 * np.pi)

    return charge_field


def func(x, a, b, c):
    return a*x[0] + b*x[1] + c

def linear_fit_2(kphis,Ls,Wc,outname,planar=False):

    Wc = np.transpose(Wc)
    X1, X2 = np.meshgrid(kphis,Ls)     
    size = X1.shape
    x1_1d = X1.reshape((1, np.prod(size)))
    x2_1d = X2.reshape((1, np.prod(size)))
   
    xdata = np.vstack((x1_1d, x2_1d))
    
    guess = (2.0, 0.1875, 0)
    if planar:
        guess = (4.0, 0.3750, 0)
        
    ydata = Wc.reshape(np.prod(size))
    popt, pcov = curve_fit(func, xdata, ydata)
    print("original: {}\nfitted: {}".format(guess, popt))
    
    with open(outname+'.txt', 'w') as f:
        f.write("alpha   = {:.3f}\n".format(popt[0]))
        f.write("beta    = {:.3f}\n".format(popt[1]))
        f.write("constant = {:.3f}\n\n".format(popt[2]))
    
    z_fit = func(xdata, *popt)
    Z_fit = z_fit.reshape(Wc.shape) 
    
    fig, ax = plt.subplots(1,2,figsize=(6,3))
    c1 = ax[0].pcolormesh(X1, X2, Wc, shading='auto')
    c2 = ax[1].pcolormesh(X1, X2, Z_fit, shading='auto')

    fig.colorbar(c1, ax=ax[0])
    fig.colorbar(c2, ax=ax[1])
    ax[0].set_title(r'Critical homeotropic anchoring, $W^c_\perp$')
    ax[0].set_xlabel(r'surface tension, $k_\phi$')
    ax[0].set_ylabel(r'elasticity, $L_1$')
    ax[1].set_title(r'Critical planar anchoring, $W^c_\parallel$')
    ax[1].set_xlabel(r'surface tension, $k_\phi$')
    ax[1].set_ylabel(r'elasticity, $L_1$')
    plt.tight_layout()
    plt.show()
    
    f_flat = Wc.ravel()
    f_pred = Z_fit.ravel()
    residuals = f_flat - f_pred

    MSE = np.mean(residuals**2)
    RMSE = np.sqrt(MSE)
    MAE = np.mean(np.abs(residuals))
    MAX_ERR = np.max(np.abs(residuals))

    # R^2 score
    SS_res = np.sum(residuals**2)
    SS_tot = np.sum((f_flat - np.mean(f_flat))**2)
    R2 = 1 - SS_res/SS_tot

    with open(outname+'.txt', 'a') as f:
        f.write("MSE     = {:.3f}\n".format(MSE))
        f.write("RMSE    = {:.3f}\n".format(RMSE))
        f.write("MAE     = {:.3f}\n".format(MAE))
        f.write("Max Err = {:.3f}\n".format(MAX_ERR))
        f.write("R^2     = {:.3f}\n".format(R2))
    
    return 

def linear_fit(Apar, Bpar, data, outname, visualization=False):

    A_grid, B_grid = np.meshgrid(Apar, Bpar, indexing='ij')
    print(np.shape(A_grid),np.shape(data))
    f_grid = data
    A_flat = A_grid.ravel()
    B_flat = B_grid.ravel()
    f_flat = f_grid.ravel()

    # matrix for linear solver
    X = np.column_stack([A_flat, B_flat, np.ones_like(A_flat)])
    # solve with least squares
    coeffs, residuals, rank, s = np.linalg.lstsq(X, f_flat, rcond=None)
    alpha, beta, constant = coeffs

    with open(outname+'.txt', 'w') as f:
        f.write("alpha   = {:.3f}\n".format(alpha))
        f.write("beta    = {:.3f}\n".format(beta))
        f.write("constant = {:.3f}\n\n".format(constant))

    # compute errors
    f_pred = X @ coeffs
    residuals = f_flat - f_pred

    MSE = np.mean(residuals**2)
    RMSE = np.sqrt(MSE)
    MAE = np.mean(np.abs(residuals))
    MAX_ERR = np.max(np.abs(residuals))

    # R^2 score
    SS_res = np.sum(residuals**2)
    SS_tot = np.sum((f_flat - np.mean(f_flat))**2)
    R2 = 1 - SS_res/SS_tot

    with open(outname+'.txt', 'a') as f:
        f.write("MSE     = {:.3f}\n".format(MSE))
        f.write("RMSE    = {:.3f}\n".format(RMSE))
        f.write("MAE     = {:.3f}\n".format(MAE))
        f.write("Max Err = {:.3f}\n".format(MAX_ERR))
        f.write("R^2     = {:.3f}\n".format(R2))

    # visualization
    if visualization:
        f_pred_grid = alpha*A_grid + beta*B_grid + constant
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(A_grid, B_grid, f_grid, alpha=0.7, label='Data')
        ax.plot_surface(A_grid, B_grid, f_pred_grid, alpha=0.7, color='red',label="fit")
        ax.set_xlabel(r'$k_\phi$')
        ax.legend()
        ax.set_ylabel(r'$L_1$')
        ax.set_zlabel(r'$W_c^\parallel$')
        plt.title('Critical Anchoring: Data vs Linear Fit')
        plt.show()
