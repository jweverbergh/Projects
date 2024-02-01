#epsilon = 10^(-20)

from matplotlib import pyplot as plt
from numpy import *

U = []
X = []
x = []

l = []
with open('conv_rate_residual_prec.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_residual.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_residual_prec_2.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if i==27:
            break
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_residual_prec_2_IC_2.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if i==27:
            break
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b_prec.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b_prec_2.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if i==27:
            break
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b_prec_2_IC_2.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if i==27:
            break
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)


plt.figure()

plt.plot(X[1],U[1],'oc', label='sans préconditionnement')
plt.plot(X[0],U[0],'ob', label='ILU0', markersize=7)
plt.plot(X[2],U[2],'o', color='orange', label='IC', markersize=4)
#plt.plot(X[3],U[3],'oy', label='IC-algo2', markersize=4)

plt.plot(X[5],U[5],'c')
plt.plot(X[4],U[4],'b')
plt.plot(X[6],U[6],'--',color='orange')
#plt.plot(X[7],U[7],'y--', label='IC-algo2', markersize=4)

plt.annotate("Points : ||r_n||/||r_0||\nLignes continues : ||Ax_n-b||/||b||", (9.5, 1e-7))

plt.title('Vitesse de convergence des résidus successifs \n à la première itération temporelle, pour les variations jour/nuit')
plt.legend(loc='upper right')
plt.xlabel('Itération CG [-]')
plt.ylabel('Vitesse de convergence [-]')
plt.yscale('log')

plt.show()