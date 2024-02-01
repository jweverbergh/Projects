#epsilon = 10^(-20)

from matplotlib import pyplot as plt
from numpy import *

U = []
X = []
x = []
max_l = []
max_l_1 = []
max_l_2 = []
max_l_3 = []
max_l_4 = []
max_l_5 = []
P = 60*60*24*365

l = []
maximum = 0
j=1
liste_index = []
with open('temp_day_8years_CI_15.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
            index = i
        if (i>(j*P)/(60*60*24)):
            max_l.append(maximum)
            liste_index.append(index)
            maximum = 0
            j += 1
        if i==53100: #ne sert à rien
            break
X.append([i for i in range(1,len(l)+1)])
U.append(l)
max_l.append(maximum)
liste_index.append(index)

l = []
maximum = 0
j=1
with open('temp_day_8years_year.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
        if (i>(j*P)/(60*60*24)):
            max_l_1.append(maximum)
            maximum = 0
            j += 1
        if i==53100:
            break
X.append([i for i in range(1,len(l)+1)])
U.append(l)
max_l_1.append(maximum)

l = []
maximum = 0
j = 1
with open('temp_day_8years_CI_8_007.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
        if (i>(j*P)/(60*60*24)):
            max_l_2.append(maximum)
            maximum = 0
            j += 1
        if i==53100:
            break
X.append([i for i in range(1,len(l)+1)])
U.append(l)
max_l_2.append(maximum)

l = []
maximum = 0
j=1
with open('temp_day_8years_CI_0.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
        if (i>(j*P)/(60*60*24)):
            max_l_3.append(maximum)
            maximum = 0
            j += 1
        if i==53100:
            break
X.append([i for i in range(1,len(l)+1)])
U.append(l)
max_l_3.append(maximum)

l = []
maximum = 0
j=1
with open('temp_day_8years_CI_-_5.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
        if (i>(j*P)/(60*60*24)):
            max_l_4.append(maximum)
            maximum = 0
            j += 1
        if i==53100:
            break
X.append([i for i in range(1,len(l)+1)])
U.append(l)
max_l_4.append(maximum)

l = []
maximum = 0
j=1
with open('temp_day_8years_CI_5.out', 'r') as fd:
    b = fd.readlines()
    for i in range(len(b)):
        output = b[i]
        l.append(float(output))
        if (maximum<float(output)):
            maximum = float(output)
        if (i>(j*P)/(60*60*24)):
            max_l_5.append(maximum)
            maximum = 0
            j += 1
        if i==53100:
            break
X.append([i for i in range(1,len(l)+1)])
U.append(l)
max_l_5.append(maximum)

print(max_l_1)
plt.figure()

plt.plot(liste_index, max_l, 'c', label='15')
plt.plot(liste_index, max_l_1, 'y', label='10')
#plt.plot([liste_index[0], liste_index[-1]], [12.07, 12.07], '--', label='max du pattern')
plt.plot(liste_index, max_l_2, 'black', label='8')
plt.plot(liste_index, max_l_5, "r", label='5')
#plt.plot(liste_index, max_l_3, "y", label='0')
plt.plot(liste_index, max_l_4, 'orange', label='-5')

plt.title('Evolution des maximums locaux de température \n au centre de la conduite pour les variations été/hiver, pendant 8 ans')
plt.legend(loc='upper right')
plt.xlabel('Nombre de jours [-]')
plt.ylabel('Température [°C]')
plt.xscale('log')

plt.show()