import numpy as np
from scipy.stats import kstat
import warnings
import os

def binder(x):
    fourth = np.average(x**4)
    second = np.average(x**2)
    return 1 - fourth / (3 * (second**2))

def readQtensor(Lx, Ly, Lz, filename):

    """Reads a Qtensor.txt and outputs each quantities into numpy arrays."""
    
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"{filename} not found")

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("error", category=UserWarning)
            Qtensor = np.loadtxt(filename)
    
    except UserWarning:
        raise ValueError(f"{filename} is empty or invalid (UserWarning).")

    except OSError as e:
        raise ValueError(f"{filename} could not be read: {e}")

    if Qtensor.size == 0:
        raise ValueError(f"{filename} contained no data.")

    # initialize grids
    Qxx = np.zeros((Lx, Ly, Lz), dtype=float)
    Qxy = np.zeros((Lx, Ly, Lz), dtype=float)
    Qxz = np.zeros((Lx, Ly, Lz), dtype=float)
    Qyy = np.zeros((Lx, Ly, Lz), dtype=float)
    Qyz = np.zeros((Lx, Ly, Lz), dtype=float)

    # velocity field
    ux = np.zeros((Lx, Ly, Lz), dtype=float)
    uy = np.zeros((Lx, Ly, Lz), dtype=float)
    uz = np.zeros((Lx, Ly, Lz), dtype=float)

    v = np.zeros((3, Lx, Ly, Lz), dtype=float)  # vector with highest eigenvalues
    d = np.zeros((Lx, Ly, Lz), dtype=float)  # q field (scalar orientational field)

    mu = np.zeros((Lx, Ly, Lz), dtype=float)
    density = np.zeros((Lx, Ly, Lz), dtype=float)
    phi = np.zeros((Lx, Ly, Lz), dtype=float)
    flowparameter = np.zeros((Lx, Ly, Lz), dtype=float)
    compressibilitty = np.zeros((Lx, Ly, Lz), dtype=float)

    for row in Qtensor:
        x, y, z = int(row[0]), int(row[1]), int(row[2])

        Qxx[x, y, z] = row[3]
        Qxy[x, y, z] = row[4]
        Qxz[x, y, z] = row[5]
        Qyy[x, y, z] = row[6]
        Qyz[x, y, z] = row[7]

        ux[x, y, z] = row[8]
        uy[x, y, z] = row[9]
        uz[x, y, z] = row[10]
        
        v[0, x, y, z] = row[11]
        v[1, x, y, z] = row[12]
        v[2, x, y, z] = row[13]
        d[x, y, z] = 3 * row[14] / 2
        
        mu[x, y, z] = row[15]
        density[x, y, z] = row[16]
        phi[x, y, z] = row[17]
        flowparameter[x, y, z] = row[18]
        compressibilitty[x, y, z] = row[19]

    return v[:, :, :, 0], d[:, :, 0], phi[:, :, 0], np.average(d), binder(phi), Qxx[:,:,0], Qxy[:,:,0]

# v, d, phi, d_avg, U, Qxx, Qxy = readQtensor(50,50,1,"../Wc/kphi_0.010_L1_0.010/Qtensor_phi0_1.000_gamma0_2.000_k_0.010_W_-0.030.txt")

# import matplotlib.pyplot as plt
# print(np.max(phi),np.min(phi))
# print(np.max(d),np.min(d))

# fig, axs = plt.subplots(1,2)
# axs[0].imshow(phi,interpolation="gaussian")
# axs[1].imshow(d,interpolation="gaussian")
# plt.show()
# print(np.max(d[abs(phi-1.0)<0.05]))