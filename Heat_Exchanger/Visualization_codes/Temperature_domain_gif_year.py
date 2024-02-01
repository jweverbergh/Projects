import matplotlib.cm as cm
import matplotlib.cbook as cbook
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import numpy as np
import matplotlib.pyplot as plt

#enlever vmin et vmax dans ax.imshow pour ne pas fixer le range des colors
def fct(doc, days, ind):
    l = [[] for i in range(41)]
    j=0
    with open(doc, 'r') as fd:
        b = fd.readlines()
        for i in range(len(b)):
            if (i%81 == 0 and i!=0):
                j += 1
            output = b[i]
            l[j].append(float(output))

    P = 365
    for i in range(81):
        l[-1].append((10) + (15)*np.sin(2*np.pi*days/P))

    Z= np.array(l)
    fig, ax = plt.subplots()
    im = ax.imshow(Z, interpolation='antialiased',cmap=plt.cm.jet, vmin=-5, vmax=25,
    origin='lower', extent=[0, 10, 0, 5])
    plt.title("Température après {} jours [°C]".format(days))
    plt.xlabel("L [m]")
    plt.ylabel("H [m]")
    plt.colorbar(im)
    plt.savefig("Domain_{}".format(ind))
    plt.show()
    plt.close()


liste_docs = ['temp_domain_2550days.out', 'temp_domain_2575days.out', 'temp_domain_2600days.out','temp_domain_2625days.out','temp_domain_2650days.out', 'temp_domain_2675days.out', 'temp_domain_2700days.out','temp_domain_2725days.out','temp_domain_2750days.out','temp_domain_2775days.out','temp_domain_2800days.out','temp_domain_2825days.out','temp_domain_2850days.out','temp_domain_2875days.out','temp_domain_2900days.out','temp_domain_2925days.out']
days = [2550, 2575, 2600, 2625, 2650, 2675, 2700, 2725, 2750, 2775, 2800, 2825, 2850, 2875, 2900, 2925]
ind = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

for i in range(len(liste_docs)):
    fct(liste_docs[i], days[i], ind[i])