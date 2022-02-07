import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt

# Define system size
size = 7
deltax = 1./6.

# Define constants
d = np.array([0.17] * (size + 2))
siga = np.array([0.56] * size)
nu = 2.4
sigf = np.array([.42] * size)

# Build the matrix
matrice = []
for i in range(size):
    row = [0] * size
    if i == 0:
        row[i] = - (1/(2*deltax**2)) * (d[i] + 2 * d[i+1] + d[i+2]) - siga[i]
        row[i+1] = (1/(2* deltax**2)) * (d[i+1] + d[i+2])
    elif i == size - 1:
        row[i-1] = (1/(2*deltax**2)) * (d[i] + d[i+1])
        row[i] = -(1/(2*deltax**2)) * (d[i] + 2 * d[i+1] + d[i+2]) - siga[i]
    else:
        row[i - 1] = (1 / (2 * deltax ** 2)) * (d[i] + d[i + 1])
        row[i] = -(1 / (2 * deltax ** 2)) * (d[i] + 2 * d[i + 1] + d[i + 2]) - siga[i]
        row[i + 1] = (1 / (2 * deltax ** 2)) * (d[i + 1] + d[i + 2])
    matrice.append(row)

matrice = np.array(matrice)
# Define our solver
sigf_m = []
for i in range(size):
    row = [0]*size
    row[i] = sigf[i]*nu
    sigf_m.append(row)
sigf_m = np.array(sigf_m)

#idsigf = np.identity(size) * sigf * nu

w, v = linalg.eig(matrice, b=sigf)
#print(sum(v[:,0])*deltax)
plt.plot(v[:,0])
plt.show()

plt.plot([i for i in v[:,0]])
plt.plot([i for i in v[:,1]])
plt.show()

plt.plot(list(range(size)), [i for i in v[:,0]])
plt.plot(list(range(size)), [i for i in v[:,1]])
plt.plot(list(range(size)), [i for i in v[:,2]])
plt.show()

plt.plot(list(range(size)), [i for i in v[:,0]])
plt.plot(list(range(size)), [i for i in v[:,1]])
plt.plot(list(range(size)), [i for i in v[:,2]])
plt.plot(list(range(size)), [i for i in v[:,3]])
plt.plot(list(range(size)), [i for i in v[:,4]])
plt.plot(list(range(size)), [i for i in v[:,5]])
plt.show()

#plt.plot(v[0])
#plt.show()