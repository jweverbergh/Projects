import numpy as np
import matplotlib.pyplot as plt
import math


## Functions


#return the value at t of the i-th B-spline of degree d with knots T
def B_spline_recursive(deg, T, i, t):
    if deg==0:
        if t>=T[i] and t<T[i+1]:
            return 1
        else:
            return 0
    else:
        first = 0
        second = 0
        den_first = T[i+deg]-T[i]
        den_second = T[i+1+deg]-T[i+1]
        if den_first!=0:
            first = (t-T[i])/den_first * B_spline_recursive(deg-1, T, i, t)
        if den_second!=0:
            second = (T[i+1+deg]-t)/den_second * B_spline_recursive(deg-1, T, i+1, t)
        return first + second

#return the coordinates of the point on the curve at t;
#the curve is composed with B-splines of degree deg at
#knots T, with control points P
def curve_with_B_splines(deg, T, t, P):
    coord = np.zeros(len(P[0]))
    for i in range(len(P)):
        coord += B_spline_recursive(deg, T, i, t) * P[i]
    return coord

#compute mu
#if param = 1 : compute the approximate weight 2/C described in section 3 in the article
#if param = 2 : compute the mu found by exact line search
def compute_mu(A, param):
    n = len(A[0])-1
    m = len(A.T[0])-1
    
    if param==1:
        c = np.zeros(n+1)
        for i in range(n+1):
            for j in range(m+1):
                c[i] += A[j][i]
        C = max(c)
        return 2/C
    return 0.1

#compute the error between Q and b-spline curves
def compute_error(Q, P, A):
    err = 0
    m = len(Q)-1
    n = len(P)-1
    
    for i in range(m+1):
        err += np.linalg.norm(Q[i] - A[i]@P)**2
    return err

#LSPIA method on data Q with m+1 data points,
#taking n+1 control points and B-splines of degree deg
#return P the control points from the k-th iteration, and T the knots
#param gives the choice of mu
def algo(m, n, deg, Q, k, param, normal):
    dim = 2 #only 2 is possible, otherwise we must change Q defined later on
    
    #initialization of control points
    P = np.zeros((n+1, dim))
    P[0] = Q[0]
    P[-1] = Q[-1]
    for i in range(1,n):
        j = math.floor((m+1)*i/n)
        P[i] = Q[j]

    #initialization of parameters for Q
    t_Q = np.zeros(len(Q))
    t_Q[0] = 0
    t_Q[-1] = 1
    D = 0
    for i in range(1,m+1):
        D += np.linalg.norm(Q[i]-Q[i-1])
    for i in range(1,m):
        t_Q[i] = t_Q[i-1] + np.linalg.norm(Q[i]-Q[i-1])/D

    #initialization of knots for the B-splines
    T = np.zeros(len(P)+deg+1)
    d = (m+1)/(n-deg+1)
    for i in range(deg+1):
        T[i] = 0
        T[len(T)-1-i] = 1
    for j in range(1,n-deg+1):
        i = math.floor(j*d)
        alpha = j*d - i
        T[j+deg] = (1-alpha)*t_Q[i-1] + alpha*t_Q[i]

    #initialization of A
    A = np.zeros((m+1, n+1))
    for i in range(m+1):
        for j in range(n+1):
            A[i][j] = B_spline_recursive(deg, T, j, t_Q[i])
    A[-1][-1] = 1 #since the last line of A is full of 0 (by definition of the bsplines at 1) and the sum of all b-splines must be 1
    
    #verification of Schoenberg-Whitney condition
    for i in range(n+1):
        boolean = False
        for j in range(m+1):
            if A[j][i] !=0 :
                boolean = True
        if boolean == False:
            print('error : Schoenberg-Whitney condition not verified')
    
    #iterations
    errors = [10]
    it = 0
    while True:
        #compute error
        err = compute_error(Q, P, A)
        errors.append(err)
        
        #stop criterion : number of iterations
        if k!=None:
            if k==it:
                break
            it+=1
        #stop criterion : decrease of error
        else:
            if np.abs(errors[-1]-errors[-2])<1e-7:
                break
        
        #algorithm
        mu = compute_mu(A, param)
        P = P + mu*A.T@(Q-A@P)

    #solve the normal equations
    P_sol=0
    err_sol=0
    if normal ==1:
        P_sol = np.linalg.solve(A.T@A, A.T@Q)
        err_sol = compute_error(Q, P_sol, A)

    return P, T, errors, P_sol, err_sol

#plot
def plot(m, n, deg, Q, k, param, normal, plot_err):
    #apply LSPIA method and solve normal equations
    P, T, errors, P_sol, err_sol = algo(m, n, deg, Q, k, param, normal)
    dim = 2
    
    #compute the b-spline curve for the sol of LSPIA method
    t = np.linspace(T[deg],T[n+1]-0.0001,500) #the curve is not defined on T[n+1]
    u = np.zeros((len(t),dim))
    for j in range(len(t)):
        u[j] = curve_with_B_splines(deg, T, t[j], P)
    
    #compute the b-spline curve for the sol of normal equations
    if normal==1:
        u_sol = np.zeros((len(t),dim))
        for j in range(len(t)):
            u_sol[j] = curve_with_B_splines(deg, T, t[j], P_sol)
    
    
    #plot the b-spline curve
    plt.plot(u[:,0], u[:,1], color='red')
    plt.scatter(Q.T[0], Q.T[1], s=15, color='blue')
    plt.scatter(P.T[0], P.T[1], s=25, color='orange')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Step {}'.format(len(errors)-2))
    plt.show()
    
    #plot the solution of normal equations
    if normal == 1:
        plt.plot(u_sol[:,0], u_sol[:,1], color='red')
        plt.scatter(Q.T[0], Q.T[1], s=15, color='blue')
        plt.scatter(P_sol.T[0], P_sol.T[1], s=25, color='orange')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Least squares solution')
        plt.show()
    
    #plot the error at each iteration
    if plot_err == 1:
        plt.scatter(np.arange(len(errors)-1), errors[1:], s=15, label='LSPIA method')
        plt.axhline(y = err_sol, color = 'r', linestyle = '-', label='least squares solution')
        plt.xlabel('k-th iteration')
        plt.ylabel('$E_k$')
        plt.title('Error at each iteration')
        plt.legend()
        plt.show()
    
    #plot the error at each iteration
    if plot_err == 1:
        plt.scatter(np.arange(len(errors)-1), errors[1:], s=15, label='LSPIA method')
        plt.axhline(y = err_sol, color = 'r', linestyle = '-', label='least squares solution')
        plt.plot(np.linspace(1,len(errors)-2), 1.5/np.linspace(1,len(errors)-2), label = '1.5/k')
        plt.plot(np.linspace(1,len(errors)-2), 5/(np.linspace(1,len(errors)-2))**2,color='green', label='5/$k^2$')
        plt.xlabel('k-th iteration')
        plt.ylabel('$E_k$')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Error at each iteration')
        plt.legend()
        plt.show()
    
#create the curve to be fitted : letter J (divided into 28 segments)
def Q_create(m):
    Q = np.zeros((m+1,2))
    s = (m+1)/28  #number of points in each segment
    for i in range(m+1):
        if i<4*s: #4 segments for the upper long bar of J
            Q[i][0] = i/(4*s)
            Q[i][1] = 1
        elif i<8*s: #4 segments for the right long bar of J
            Q[i][0] = 1
            Q[i][1] = 1 - (i-4*s)/(4*s)
        elif i<16*s: #8 segments for the lower demi-circle of J
            Q[i][0] = 0.5 + 0.5*np.cos((i-8*s)*np.pi/(8*s))
            Q[i][1] = -0.5*np.sin((i-8*s)*np.pi/(8*s))
        elif i<17*s: #1 segment for the little horizontal bar of J
            Q[i][0] = (i-16*s)*0.25/s
            Q[i][1] = 0
        elif i<21*s: #4 segments for the upper demi-circle of J
            Q[i][0] = 0.5 + 0.25*np.cos((i-17*s)*np.pi/(4*s)+np.pi)
            Q[i][1] = 0.25*np.sin((i-17*s)*np.pi/(4*s)+np.pi)
        elif i<24*s: #3 segments for the left long bar of J
            Q[i][0] = 0.75
            Q[i][1] = (i-21*s)*0.75/(3*s)
        elif i<27*s: #3 segments for the lower long bar of J
            Q[i][0] = 0.75 - (i-24*s)*0.75/(3*s)
            Q[i][1] = 0.75
        elif i<28*s: #1 segment for the little vertical bar of J
            Q[i][0] = 0
            Q[i][1] = 0.75 + (i-27*s)*0.25/s
    Q[-1][1] = 1
    return Q


## Simulations

#parameters
m = 200        #number of data points-1
n = 50         #number of control points-1
deg = 3        #degree of B-splines
k = 0          #number of iterations of LSPIA; if k=None, then iterate until the errors difference is less than 1e-7
param = 1      #if param=1, then compute mu with the approximation 2/C from the paper
normal = 0     #if =0, do nothing; if =1, then compute the solution of normal equations
plot_err = 0   #if =0, do nothing; if =1, then plot the error at each iteration

#curve to be fitted : letter J (divided into 28 segments)
Q = Q_create(m)


#produce graphs of section 4.1 in the report
k = [0,1,3,5,7]
for i in range(len(k)):
    plot(m, n, deg, Q, k[i], param, normal, plot_err)

k = None
normal = 1
plot(m, n, deg, Q, k, param, normal, plot_err)

k = 20
normal = 0
plot_err = 1
plot(m, n, deg, Q, k, param, normal, plot_err)