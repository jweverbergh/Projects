#epsilon = 10^(-20)

from matplotlib import pyplot as plt
from numpy import *

U = []
X = []
x = []

l = []
with open('conv_rate_residual_prec_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_residual_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_residual_prec_2_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_residual_prec_2_IC_2_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b_prec_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b_prec_2_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
with open('conv_rate_Ax-b_prec_2_IC_2_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
x.append(linspace(1, len(l), 100))
X.append([i for i in range(1,len(l)+1)])
U.append(l)



plt.figure()

plt.plot(X[1],U[1],'c', label='sans préconditionnement')
plt.plot(X[0],U[0],'b', label='ILU0', markersize=7)
plt.plot(X[2],U[2], color='orange', label='IC', markersize=4)
#plt.plot(X[3],U[3],'oy', label='IC-algo2', markersize=4)

plt.plot(X[5],U[5],'c')
plt.plot(X[4],U[4],'b')
plt.plot(X[6],U[6],'--',color='orange')
#plt.plot(X[7],U[7],'y--', label='IC-algo2', markersize=4)

plt.annotate("1. ||r_n||/||r_0||\n2. ||Ax_n-b||/||b||", (120, 1e-5))
plt.annotate("2", (50,1e-13))
plt.annotate("2", (160, 1e-13))
plt.annotate("1", (30,1e-16))
plt.annotate("1", (150,1e-16))

plt.title('Vitesse de convergence des résidus successifs \n à la première itération temporelle, pour les variations été/hiver')
plt.legend(loc='upper right')
plt.xlabel('Itération CG [-]')
plt.ylabel('Vitesse de convergence [-]')
plt.yscale('log')

plt.show()