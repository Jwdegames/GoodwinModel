import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.special import lambertw
#This program graphs employment and wage share, which are functions of time, from the Goodwin Model of Class Strucggle
#mu is the employment rate
#u is the wage share
#v is the capital output ratio
#theta is the natural labor productivity growth rate
#n is the natural labor supply growth rate
#alpha and beta are the parameters from a linear phillips curve
#or g_w = - alpha + beta * mu, where g_w is the growth in wages or dw/dt

#For the exponential phillips curve
# g_w = -alpha + beta*e^(k*mu)

#Linear phillips curve
def dZ_dt(Z, t, v=5, theta=0.009, n=0.075, alpha=0.6, beta = 1): # a,b,c,d optional arguments.
    u, mu = Z[0], Z[1]
    dudt, dmudt = (-(alpha+theta)+beta*mu)*u, (1/v-(theta+n)-u/v)*mu
    return [dudt, dmudt]

#exponential phillips curve
def dZ2_dt(Z, t, v=5, theta=0.009, n=0.075, alpha=0.6, beta = 1, k = 0.002): # a,b,c,d optional arguments.
    u, mu = Z[0], Z[1]
    dudt, dmudt = (-(alpha+theta)+beta*np.exp(k*mu))*u, (1/v-(theta+n)-u/v)*mu
    return [dudt, dmudt]

#A value for linear phillips curve  
def calc_A(v, theta, n, alpha):
    return (1/v - (theta+n))/(alpha+theta)

#Find first max in data set
#Data set should have a max and there should be at least 3 data points
def find_first_max(data_set):
    for i in range(1,len(data_set)-1):
        if (data_set[i-1] < data_set[i] and data_set[i+1] < data_set[i]):
            return data_set[i]

#Same as function above, but this finds the first min in the data set
def find_first_min(data_set):
    for i in range(1,len(data_set)-1):
        if (data_set[i-1] > data_set[i] and data_set[i+1] > data_set[i]):
            return data_set[i]

#Factorial
def factorial(n):
    if n <= 1:
        return 1
    else:
        return n * factorial(n-1)

#Calculate the approximation of the integral of e^(k*mu)/mu
def approximate_integral_mu(k, mu, terms):
    if (terms == 1):
        return np.log(mu)
    elif (terms == 2):
        return np.log(mu)+k*mu
    else:
        n = terms - 1
        return k**n*mu**n/(factorial(n)*(n))+approximate_integral_mu(k, mu, n)

#Calculates the value of the mu side - u side
def mu_function(coefficients, mu, u, terms_to_use):
    alpha = coefficients[0]
    theta = coefficients[1]
    beta = coefficients[2]
    k = coefficients[3]
    v = coefficients[4]
    n = coefficients[5]
    c = coefficients[6]
    mu_side = -(alpha+theta)*np.log(mu)+beta*approximate_integral_mu(k, mu, terms_to_use)
    u_side = (1/v-(theta+n))*np.log(u)-u/v+c
    return mu_side - u_side

#Finds the root of the mu function via bisection given the value of u and the upper and lower bounds of mu
def bisect_mu_root(coefficients, u, lowerBound, upperBound, TOL, terms_to_use):
    '''Given the coefficients of a cubic, the lower and upper bounds of an interval, and a tolerance, returns an 
    interval with the root such that the interval's bounds are within the tolerance'''
    iterations = 0
    #Determine the left and right values to be used for the while loop conditional
    leftEval = mu_function(coefficients, lowerBound, u, terms_to_use)
    rightEval = mu_function(coefficients, upperBound, u, terms_to_use)
    #Continue looping until we are in tolerance and a root is in the interval
    while (upperBound-lowerBound) > TOL or (not ((leftEval <= 0 and rightEval >= 0) or leftEval >= 0 and rightEval <= 0)):
        #Bisect the interval boundaries at the start of the loop, so code doesn't need to be rewritten for each conditional pathway
        p = (lowerBound+upperBound)/2
        yp = mu_function(coefficients,p, u, terms_to_use)
        #f(x) is increasing
        if (mu_function(coefficients, lowerBound, u, terms_to_use) < mu_function(coefficients, upperBound, u, terms_to_use)):
            #Reassign the boundary as needed
            if yp < 0:
                lowerBound = p
            elif yp > 0:
                upperBound = p
            else:
                return (p-TOL/2,p+TOL/2,iterations + 1)
        #f(x) is decreasing
        else:
            #Reassign the boundary as needed
            if yp < 0:
                upperBound = p
            elif yp > 0:
                lowerBound = p
            else:
                return (p-TOL/2,p+TOL/2,iterations + 1)
        iterations += 1
        #Find the left and right y values to be used for the while loop conditional
        leftEval = mu_function(coefficients, lowerBound, u, terms_to_use)
        rightEval = mu_function(coefficients, upperBound, u, terms_to_use)
    return (lowerBound, upperBound, iterations)

v=5
theta=0.009
n=0.075
alpha=0.6
beta = 1
k = 0.2
A = calc_A(v, theta, n ,alpha)
ts = np.linspace(0, 100, 175)
Z0 = [0.4, 0.55] # initial conditions for u and mu
Zs = odeint(dZ_dt, Z0, ts, args = (v,theta,n,alpha,beta))
# use optional argument ’args’ to pass parameters to dZ_dt
#args=(2,0.2,0.2,0.3,0.2)
wage_share = Zs[:,0] # first column
employment = Zs[:,1] # second column

M0 = Z0[1]*beta/(alpha+theta)
U0 = Z0[0]/(v*(alpha+theta))
c2 = -np.log(M0*U0**A)+M0+U0

# Let’s plot ’wage share’ and ’employment rate’ as a function of time
plt.plot(ts, wage_share, "+", label="wage share")
plt.plot(ts, employment, "x", label="employment rate")
#Draw the boundary lines
plt.axhline(y = 1 - v * (theta + n), c = 'g', linestyle = '--')
plt.axhline(y = (alpha+theta)*v*-A*np.real(lambertw(-np.exp((1-c2)/A)/A)), c = 'c', linestyle = '--')
plt.axhline(y = (alpha+theta)*v*-A*np.real(lambertw(-np.exp((1-c2)/A)/A,-1)), c = 'c', linestyle = '--')
plt.axhline(y = -(alpha+theta)/beta*np.real(lambertw(-np.exp((A-c2))/(A**A))), c = 'r', linestyle = '--')
plt.axhline(y = -(alpha+theta)/beta*np.real(lambertw(-np.exp((A-c2))/(A**A),-1)), c = 'r', linestyle = '--')
plt.axhline(y = (alpha+theta)/beta, c = 'y', linestyle = '--')
plt.xlabel("Time", fontsize=14)
plt.title("Goodwin's Model of Class Struggle (Linear Phillips Curve) - Parametric Time Plots")
plt.ylabel("Ratio", fontsize=14)
plt.legend();

#Graph wage share and employment together in a single cycle
fig = plt.figure(2)
plt.plot(wage_share,employment)

#Plot the equilibrium lines
plt.axvline(x = 1 - v * (theta + n), c = 'y', linestyle = '--')
plt.axhline(y = (alpha+theta)/beta, c = 'g', linestyle = '--')

#Plot the boundary lines
plt.axvline(x = (alpha+theta)*v*-A*np.real(lambertw(-np.exp((1-c2)/A)/A)), c = 'c', linestyle = '--')
plt.axvline(x = (alpha+theta)*v*-A*np.real(lambertw(-np.exp((1-c2)/A)/A,-1)), c = 'c', linestyle = '--')
plt.axhline(y = -(alpha+theta)/beta*np.real(lambertw(-np.exp((A-c2))/(A**A))), c = 'r', linestyle = '--')
plt.axhline(y = -(alpha+theta)/beta*np.real(lambertw(-np.exp((A-c2))/(A**A),-1)), c = 'r', linestyle = '--')



plt.xlabel("Wage Share", fontsize=14)
plt.ylabel("Employment Rate", fontsize=14)
plt.title("A Cycle of Goodwin's Model of Class Struggle (Linear Phillips Curve)")

#Graph the respective Phillips Curve 
fig = plt.figure(3)
unemployment = np.linspace(0,1,100)
inflation = -alpha+beta*(1-unemployment)
plt.plot(unemployment,inflation)
plt.xlabel("Unemployment")
plt.ylabel("Inflation")
plt.title("Linear Phillips Curve")
plt.axhline(y = 0, c = 'k')

#Now let's do the same for the expoential phillips curve
alpha= 1
beta = 0.1
k = 4.5
Z2s = odeint(dZ2_dt, Z0, ts, args = (v,theta,n,alpha,beta, k))
# use optional argument ’args’ to pass parameters to dZ2_dt

terms_to_use = 20
A_mu0 = approximate_integral_mu(k,Z0[1],terms_to_use)
#print("A_mu0",A_mu0)
c = -(alpha+theta)*np.log(Z0[1])+beta*A_mu0-(1/v-(theta+n))*np.log(Z0[0])+Z0[0]/v


'''
#The Code below proves the relationship between u and mu for the exponential phillips curve
#Uncomment this to test it out
print(c)
print("Equation 1:",-(alpha+theta)*np.log(Z0[1])+beta*A_mu0,"=",(1/v-(theta+n))*np.log(Z0[0])-Z0[0]/v+c)
print("Equation 2:",1/(1/v-(theta+n))*(-(alpha+theta)*np.log(Z0[1])+beta*A_mu0-c),"=",np.log(Z0[0])+Z0[0]/(v*(n+theta)-1))
print("Equation 3:",np.exp(np.log(Z0[1]**(-(alpha+theta)/(1/v-(theta+n)))))*np.exp((beta*A_mu0-c)/(1/v-(theta+n))),"=",np.exp(np.log(Z0[0]))*np.exp(Z0[0]/(v*(theta+n)-1)))

f_mu = 1/v-(theta+n)
f_mu2 = v*(theta+n)-1
print("Equation 4:",Z0[0]/f_mu2*np.exp(Z0[0]/f_mu2),"=",Z0[1]**(-(alpha+theta)/f_mu)*np.exp((beta*A_mu0-c)/f_mu)/f_mu2)
print("Equation 5:",Z0[0]/f_mu2,"=",np.real(lambertw(Z0[1]**(-(alpha+theta)/f_mu)*np.exp((beta*A_mu0-c)/f_mu)/f_mu2)))
u_lamb = Z0[1]**(-(alpha+theta)/f_mu)*np.exp((beta*A_mu0-c)/f_mu)/f_mu2
print("Final Equation:",Z0[0],"=",f_mu2*np.real(lambertw(u_lamb)))
'''

#args=(2,0.2,0.2,0.3,0.2)
wage_share2 = Z2s[:,0] # first column
employment2 = Z2s[:,1] # second column
fig = plt.figure(4)
plt.plot(ts, wage_share2, "+", label="wage share")
plt.plot(ts, employment2, "x", label="employment rate")

#Plot equilibrium lines
plt.axhline(y = 1 - v * (theta + n), c = 'y', linestyle = '--')
plt.axhline(y = np.log((alpha+theta)/beta)/k, c = 'g', linestyle = '--')


equi_mu = np.log((alpha+theta)/beta)/k
A_mu = approximate_integral_mu(k,equi_mu,terms_to_use)
#print('equi_mu',equi_mu)
#print(A_mu)

#Common factors used in equation
f_mu = (1/v-(theta+n))
f_mu2 = v*(theta+n)-1

#Store the bounds as variables
u_lamb = equi_mu**(-(alpha+theta)/f_mu)*np.exp((beta*A_mu-c)/f_mu)/f_mu2
#print("Inner Lambert:",u_lamb)
u_lower_bound = f_mu2*np.real(lambertw(u_lamb))
u_upper_bound = f_mu2*np.real(lambertw(u_lamb,-1))
coefficients = [alpha, theta, beta, k, v, n, c]

TOL = 1e-5
mu_lower_bound_tuple = bisect_mu_root(coefficients, 1-v*(theta+n), 1e-10, equi_mu, TOL, terms_to_use)
mu_lower_bound = (mu_lower_bound_tuple[0]+mu_lower_bound_tuple[1])/2
mu_upper_bound_tuple = bisect_mu_root(coefficients, 1-v*(theta+n), equi_mu, 1, TOL, terms_to_use)
mu_upper_bound = (mu_upper_bound_tuple[0]+mu_upper_bound_tuple[1])/2

#print(mu_lower_bound)

#print(u_lower_bound)
#print(u_upper_bound)

#print(find_first_min(wage_share2))
#print(find_first_max(wage_share2))
#print(find_first_min(employment2))
#print(find_first_max(employment2))

#Actually plot
plt.axhline(y = u_lower_bound, c = 'c', linestyle = '--')
plt.axhline(y = u_upper_bound, c = 'c', linestyle = '--')
plt.axhline(y = mu_lower_bound, c = 'r', linestyle = '--')
plt.axhline(y = mu_upper_bound, c = 'r', linestyle = '--')

plt.xlabel("Time", fontsize=14)
plt.title("Goodwin's Model of Class Struggle (Exponential Phillips Curve) - Parametric Time Plots")
plt.ylabel("Ratio", fontsize=14)
plt.legend();

#Graph wage share and employment together in a single cycle
fig = plt.figure(5)
plt.plot(wage_share2,employment2)

#Plot the equilibrium lines
plt.axvline(x = 1 - v * (theta + n), c = 'y', linestyle = '--')
plt.axhline(y = np.log((alpha+theta)/beta)/k, c = 'g', linestyle = '--')

#Plot the boundary lines
plt.axvline(x = u_lower_bound, c = 'c', linestyle = '--')
plt.axvline(x = u_upper_bound, c = 'c', linestyle = '--')
plt.axhline(y = mu_lower_bound, c = 'r', linestyle = '--')
plt.axhline(y = mu_upper_bound, c = 'r', linestyle = '--')
plt.xlabel("Wage Share", fontsize=14)
plt.ylabel("Employment Rate", fontsize=14)
plt.title("A Cycle of Goodwin's Model of Class Struggle (Exponential Phillips Curve)")

#Graph the respective Phillips Curve 
fig = plt.figure(6)
#unemployment = np.linspace(0,1,100)
inflation2 = -alpha+beta*np.exp(k*(1-unemployment))
plt.plot(unemployment,inflation2)
plt.xlabel("Unemployment")
plt.ylabel("Inflation")
plt.title("Exponential Phillips Curve")
plt.axhline(y = 0, c = 'k')
