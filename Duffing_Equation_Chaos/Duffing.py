import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

##Duffing system
def sys(var, t, a, b):
    x, x_dot = var
    return (x_dot, -a[0]*x_dot + x - x**3 + a[1]*np.cos(a[2]*t))


##parameters of the Duffing equation
delta = 0.3
F = 0.305
omega = 1.2
vec = [delta, F, omega]

##parameters of the simulation
T = 2*np.pi/omega               #period of the sinusoidal excitation
n = 150                         #number of intervals dt in one period T
dt = T/n                        #time step
K = 1500                        #number of periods T
frac = 0.3                      #fraction of the iterates we do not
                                #plot (to skip the transient response)

t = np.linspace(0,K*T,K*n)      #time abscissae
x0 = [1,0]                      #initial point


##solve the system
sol = integrate.odeint(sys, x0, t, args=(vec,0))

"""
##solution (t, x)
plt.plot(t/T, sol[:,0])
plt.xlabel('t/T')
plt.ylabel('x')
plt.savefig('x_t.pdf', bbox_inches='tight')
plt.show()


##trajectory (x, x_dot)
plt.plot(sol[int(frac*K)*n:,0], sol[int(frac*K)*n:,1])
plt.scatter(list([sol[i,0]] for i in range(int(frac*K)*n, K*n, n)), \
            list([sol[i,1]] for i in range(int(frac*K)*n, K*n, n)), \
            color="red",label = "(x_k, y_k) from Poincare map")
plt.scatter(1,0, color='green',label="Initial point (1,0)")
plt.xlabel("x")
plt.ylabel("x_dot")
plt.legend()
plt.savefig('traj.pdf', bbox_inches='tight')
plt.show()


##Poincare map (x_k, x_{k+1})
x_0 = np.linspace(0,1,100)                  #initial points between 0 and 1
x_0 = x_0.reshape((len(x_0),1))
y_0 = np.zeros(len(x_0)).reshape((len(x_0),1))
x0 = np.concatenate([x_0,y_0], axis=1)      #initial iterates

for i in range(len(x0)):                    #solve system for each initial point
    period = 1                              #number of times the Poincare map is
                                            #applied (period=2 for P^2)
    sol = integrate.odeint(sys, x0[i], t, args=(vec,0))
    plt.scatter(list([sol[i,0]] for i in range(int(frac*K)*n,\
                                               K*n-period*n, period*n)),\
                np.array(list([sol[i,0]] for i in range(int(frac*K)*n+period*n, \
                                                        K*n, period*n))), s = 10)

plt.scatter(list([sol[i,0]] for i in range(int(frac*K)*n, K*n-period*n, \
                                           period*n))[-1], \
            np.array(list([sol[i,0]] for i in range(int(frac*K)*n+period*n,\
                                                    K*n, period*n)))[-1], \
            s = 60, color='black', label='fixed points')
plt.plot(1.2*x_0,1.2*x_0, label='bissector')    
plt.xlabel("$x_k$")
plt.ylabel("$x_{k+1}$")
plt.legend()
plt.savefig('poincare_map_w12_F029.pdf', bbox_inches='tight')
plt.show()


##strange attractor at phi
l=1/6                           #fraction such that phi=l*T
plt.scatter(list([sol[i,0]] for i in range(int(frac*K)*n+ int(l*n), \
                                           K*n, n)), \
            list([sol[i,1]] for i in range(int(frac*K)*n+ int(l*n), \
                                           K*n, n)), s = 10)
plt.xlabel('x_k')
plt.ylabel('y_k')
plt.savefig('strange_attractor_6_6.pdf', bbox_inches='tight')
plt.show()
"""

##strange attractor animation   (source : https://towardsdatascience.
#com/intro-to-dynamic-visualization-with-python-animations-and-
#interactive-plots-f72a7fb69245)
fig = plt.figure()
ax = fig.add_subplot(111)
colors = plt.get_cmap('winter',150)
x, y, c = np.random.random((3, 10))
f_d = ax.scatter([], [], s=10)
phi = np.arange(n)
plt.xlim(-2, 2)
plt.ylim(-1, 1)

def animate(i):
    phi_i = phi[i]
    x = list(sol[i,0] for i in range(int(frac*K)*n+phi_i, K*n-n+phi_i, n))
    y = list(sol[i,1] for i in range(int(frac*K)*n+phi_i, K*n-n+phi_i, n))
    x = np.array(x).reshape((len(x),1))
    y = np.array(y).reshape((len(x),1))
    array = np.concatenate([x,y], axis=1)    
    f_d.set_offsets(array)
    f_d.set_color(colors(i))
    return f_d,

#decomment the following to have the animation by running directly, and comment the next block of the code
"""
animation.FuncAnimation(fig, animate, blit=True, frames=len(phi),\
                        interval=10, repeat=True)

plt.xlabel('x_k')
plt.ylabel('y_k')

plt.show()
"""

plt.xlabel('x_k')
plt.ylabel('y_k')
for j in range(0,150):
    animate(j)
    #plt.savefig("strange_attractor_" + str(j))
plt.show()

"""
##bifurcation diagram
F = np.linspace(0.25, 0.55, 3000)
x0 = [1,0] 
for i in range(len(F)):
    vec[1] = F[i]
    sol = integrate.odeint(sys, x0, t, args=(vec,0))
    strobes_per_period = list([sol[i,0]] for i in range(int(frac*K)*n, \
                                                        K*n, n))
    
    plt.scatter(F[i]*np.ones(len(strobes_per_period)), strobes_per_period,\
                s = 0.7, color = "black")

plt.xlabel('F')
plt.ylabel('x_k')
plt.savefig("bifurcation_diagram.pdf", bbox_inches='tight')
plt.show()

"""
