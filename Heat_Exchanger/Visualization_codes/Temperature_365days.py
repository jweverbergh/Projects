#epsilon = 10^(-20)

from matplotlib import pyplot as plt
from numpy import *

U = []
X = []
x = []

l = []
with open('temp_hour_365days_day.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if i==24*200: #stop après 200 jours
            break
x.append(linspace(1, len(l), 100))
X.append([i/24 for i in range(1,len(l)+1)])
U.append(l)


plt.figure()

plt.plot(X[0],U[0],'b', label='au centre de la conduite')

plt.title('Evolution de la température au centre de la conduite\n pour les variations jour/nuit pendant 200 jours')
#plt.legend(loc='upper right')
plt.xlabel('Nombre de jours [-]')
plt.ylabel('Température [°C]')

plt.show()