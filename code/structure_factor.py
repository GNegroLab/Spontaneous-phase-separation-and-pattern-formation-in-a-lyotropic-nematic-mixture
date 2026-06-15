import numpy as np
import sys
from reader import readQtensor

def structure_factor(phi):

    Lx, Ly = np.shape(phi)

    phi -= np.average(phi)

    C_k = np.fft.fft2(phi)
    C_k = np.fft.fftshift(C_k)
    C_k = np.real(C_k)**2 + np.imag(C_k)**2
    #C_k /= np.sum(phi)

    kx = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Lx, d=1))  # x frequencies f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
    ky = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(Ly, d=1))  # d is the units of the spacing between gridpoints
                                                              # convert form frequency (f, cycles per unit length) to wave vector k
    return kx, ky, C_k

filepath = sys.argv[1]

try:
    Lx = int(sys.argv[2])
    Ly = int(sys.argv[3])
except IndexError:
    print("Error: Not enough arguments provided!", file=sys.stderr)
    sys.exit(1)
except ValueError:
    print("Error: The inputs provided were not valid integers!", file=sys.stderr)
    sys.exit(1)
    
v, d, phi, d_avg, U = readQtensor(Lx,Ly,1,filepath)
kx, ky, C_k = structure_factor(phi)

output_matrix = np.zeros((C_k.shape[0] + 1, C_k.shape[1] + 1))

output_matrix[1:, 1:] = C_k
output_matrix[1:, 0] = ky
output_matrix[0, 1:] = kx
output_matrix[0, 0] = len(kx)

output_filepath = "/Users/s2862303/PhD/NLBF/code/structure_factor_tmp.dat"
np.savetxt(output_filepath, output_matrix, fmt='%.4e', delimiter=' ')