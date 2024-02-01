from matplotlib import pyplot as plt
from numpy import *

U = []
X = []
seuil = 2919 #pour prendre 8 ans uniquement

P = 365
surface = []
for i in range(8*P):
    t = i
    surface.append((10) + (15)*sin(2*pi*t/P))


l = []
maximum = 0
with open('temp_day_8years_CI_15.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
            ind = i
        if i==seuil:
            break
X.append([i for i in range(1,len(l)+1)])
U.append(l)

l = []
maximum = 0
with open('temp_day_8years_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
            ind = i
        if i==seuil:
            break
U.append(l)

l = []
maximum = 0
with open('temp_day_8years_CI_8_007.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
            ind = i
        if i==seuil:
            break
U.append(l)

l = []
maximum = 0
with open('temp_day_8years_CI_5.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
            ind = i
        if i==seuil:
            break
U.append(l)

l = []
maximum = 0
with open('temp_day_8years_CI_0.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
            ind = i
        if i==seuil:
            break
U.append(l)

l = []
maximum = 0
with open('temp_day_8years_CI_-_5.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
            ind = i
        if i==seuil:
            break
U.append(l)


T_max = 12.07
T_min = 7.85
P = 365
x = linspace(1, X[0][-1], 1000)
u = []
for i in range(len(x)):
    u.append(10 + 2.07*sin((2*pi*x[i])/(P)+ 4.3767))


plt.figure()

plt.plot(X[0], surface, label="en surface (f(t))")
#plt.plot(X[0],U[0],'c', label='CI : 15°C')
plt.plot(X[0],U[1], label='au centre de la conduite')
#plt.plot(X[0],U[2],'black', label='CI : 8°C')
#plt.plot(X[0],U[3],'r', label='CI : 5°C')
#plt.plot(X[0],U[4],'y', label='CI : 0°C')
#plt.plot(X[0],U[5],'orange', label='CI : -5°C')
plt.plot([X[0][0], X[0][-1]],[7.85, 7.85], 'purple')
plt.plot([X[0][0], X[0][-1]],[12.07, 12.07], 'purple')
#plt.plot(x,u,'--', color="greenyellow", label = 'pattern')


plt.title('Evolution de la température \n pour les variations été/hiver pendant 8 ans')
plt.legend(loc='lower right')
plt.xlabel('Nombre de jours [-]')
plt.ylabel('Température [°C]')

plt.show()