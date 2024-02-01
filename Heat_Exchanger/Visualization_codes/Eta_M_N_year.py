from matplotlib import pyplot as plt
from numpy import *

U = []
X = []
x = []
eff = []
ind = []
l=[]
with open("eta_M_year.out", 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        ind.append(20*(i+1))


plt.figure()

plt.plot(ind, l, label="N=20")
plt.title('Efficacité en fonction de M\npour les variations été/hiver')
plt.xlabel('Valeur de M [-]')
plt.ylabel('Efficacité [-]')
plt.legend()
plt.show()


l=[]
ind = []
with open("eta_N_year.out", 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        ind.append(10*(i+1))
with open("eta_N_year_suite.out", 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        ind.append(140+40*(i+1))


plt.figure()

plt.plot(ind, l, label="M=20")
plt.plot([ind[0], ind[-1]],[0.14005, 0.14005], 'purple')
plt.title('Efficacité en fonction de N\npour les variations été/hiver')
plt.xlabel('Valeur de N [-]')
plt.ylabel('Efficacité [-]')
plt.legend()
plt.show()
