from numpy import linspace, exp, log, real
from scipy.integrate import odeint
from scipy.special import lambertw
from Plot import Plot

class GoodwinModel:
    '''Contains data about the Goodwin Model'''

    def __init__(self):
        self.Z0 = [0.4, 0.55]
        self.v = 5
        self.theta = 0.009
        self.n = 0.075
        self.alpha = 0.6
        self.beta = 1
        self.k = 0.2


    def updateVals(self, u0, mu0, v, theta, n, alpha, beta, k):
        '''Updates the values'''
        self.Z0 = [u0, mu0]
        self.v = v
        self.theta = theta
        self.n = n
        self.alpha = alpha
        self.beta = beta
        self.k = k


    def updatePlots(self):
        '''Updates the plots'''
        self.clearPhillipsPlot()
        self.populatePhillipsPlot()
        self.phillipsPlot.plot()
        self.phillipsPlot.draw()
        # print ("Updated Phillips Plot")

        self.clearParametricPlot()
        self.populateParametricPlot()
        self.parametricPlot.show()
        self.parametricPlot.showLegend()
        self.parametricPlot.draw()

        self.clearGoodwinPlot()
        self.populateGoodwinPlot()
        self.goodwinPlot.show()
        self.goodwinPlot.showLegend(bbox_to_anchor=(1, 0.5))
        self.goodwinPlot.draw()


    def dZ_dt(self, Z, t, v=5, theta=0.009, n=0.075, alpha=0.6, beta=1):
        '''Used to get the necessary coordinates for wage share and employment (linear version)'''
        u, mu = Z[0], Z[1]
        dudt, dmudt = (-(alpha + theta) + beta * mu) * u, (1 / v - (theta + n) - u / v) * mu
        return [dudt, dmudt]


    # exponential phillips curve
    def dZ2_dt(self, Z, t, v=5, theta=0.009, n=0.075, alpha=0.6, beta=1, k=0.002):
        '''Used to get the necessary coordinates for wage share and employment (exponential version)'''
        u, mu = Z[0], Z[1]
        dudt, dmudt = (-(alpha + theta) + beta * exp(k * mu)) * u, (1 / v - (theta + n) - u / v) * mu
        return [dudt, dmudt]


    def calc_A(self, v, theta, n, alpha):
        '''Calculates the A value for linear phillips curve'''
        return (1 / v - (theta + n)) / (alpha + theta)


    def find_first_max(self, data_set):
        '''Find first max in data set (used for debugging)
        Data set should have a max and there should be at least 3 data points'''
        for i in range(1, len(data_set) - 1):
            if (data_set[i - 1] < data_set[i] and data_set[i + 1] < data_set[i]):
                return data_set[i]


    def find_first_min(self, data_set):
        '''Find first min in data set (used for debugging)
        Data set should have a min and there should be at least 3 data points'''
        for i in range(1, len(data_set) - 1):
            if (data_set[i - 1] > data_set[i] and data_set[i + 1] > data_set[i]):
                return data_set[i]


    def factorial(self, n):
        '''Finds the factorial of n'''
        if n <= 1:
            return 1
        else:
            return n * self.factorial(n - 1)


    def approximate_integral_mu(self, k, mu, terms):
        '''Calculate the approximation of the integral of e^(k*mu)/mu'''
        if (terms == 1):
            return log(mu)
        elif (terms == 2):
            return log(mu) + k * mu
        else:
            # General formula
            n = terms - 1
            return k ** n * mu ** n / (self.factorial(n) * (n)) + self.approximate_integral_mu(k, mu, n)


    def mu_function(self, coefficients, mu, u, terms_to_use):
        '''Calculates the value of the mu side - u side'''
        alpha = coefficients[0]
        theta = coefficients[1]
        beta = coefficients[2]
        k = coefficients[3]
        v = coefficients[4]
        n = coefficients[5]
        c = coefficients[6]
        mu_side = -(alpha + theta) * log(mu) + beta * self.approximate_integral_mu(k, mu, terms_to_use)
        u_side = (1 / v - (theta + n)) * log(u) - u / v + c
        return mu_side - u_side


    # Finds the root of the mu function via bisection given the value of u and the upper and lower bounds of mu
    def bisect_mu_root(self, coefficients, u, lowerBound, upperBound, TOL, terms_to_use):
        '''Given the coefficients of a cubic, the lower and upper bounds of an interval, and a tolerance, returns an
        interval with the root such that the interval's bounds are within the tolerance'''
        iterations = 0
        # Determine the left and right values to be used for the while loop conditional
        leftEval = self.mu_function(coefficients, lowerBound, u, terms_to_use)
        rightEval = self.mu_function(coefficients, upperBound, u, terms_to_use)
        # Continue looping until we are in tolerance and a root is in the interval
        while (upperBound - lowerBound) > TOL or (
                not ((leftEval <= 0 and rightEval >= 0) or leftEval >= 0 and rightEval <= 0)):
            # Bisect the interval boundaries at the start of the loop, so code doesn't need to be rewritten for each conditional pathway
            p = (lowerBound + upperBound) / 2
            yp = self.mu_function(coefficients, p, u, terms_to_use)
            # f(x) is increasing
            if (self.mu_function(coefficients, lowerBound, u, terms_to_use) < self.mu_function(coefficients, upperBound,
                                                                                               u, terms_to_use)):
                # Reassign the boundary as needed
                if yp < 0:
                    lowerBound = p
                elif yp > 0:
                    upperBound = p
                else:
                    return (p - TOL / 2, p + TOL / 2, iterations + 1)
            # f(x) is decreasing
            else:
                # Reassign the boundary as needed
                if yp < 0:
                    upperBound = p
                elif yp > 0:
                    lowerBound = p
                else:
                    return (p - TOL / 2, p + TOL / 2, iterations + 1)
            iterations += 1
            # Find the left and right y values to be used for the while loop conditional
            leftEval = self.mu_function(coefficients, lowerBound, u, terms_to_use)
            rightEval = self.mu_function(coefficients, upperBound, u, terms_to_use)
        return (lowerBound, upperBound, iterations)


    def makePhillipsPlot(self, parent, exponential):
        '''Makes the Phillips Curve Plot - should only be called once'''
        self.parent = parent
        self.exponential = exponential
        self.unemployment = linspace(0, 1, 100)
        self.phillipsPlot = Plot(parent, 5, 4, 100)
        self.populatePhillipsPlot()
        self.phillipsPlot.figure.tight_layout()
        return self.phillipsPlot


    def clearPhillipsPlot(self):
        '''Clears Phillips Plot'''
        self.phillipsPlot.clear()


    def populatePhillipsPlot(self):
        '''Adds all necessary data to Phillips Plot'''
        # 1 - unemployment rate = employment rate
        self.inflation = -self.alpha + self.beta * (1 - self.unemployment)
        self.inflation2 = -self.alpha + self.beta * exp(self.k * (1 - self.unemployment))
        if self.exponential:
            self.phillipsPlot.setY(self.inflation2)
            print("Making exponential plot")
        else:
            self.phillipsPlot.setY(self.inflation)
        self.phillipsPlot.setX(self.unemployment)
        # Add axes and title
        self.phillipsPlot.setTitle("Phillips Curve")
        self.phillipsPlot.setAxesTitles("Unemployment", "Inflation")
        # Add reference line
        self.phillipsPlot.makeHLine(0, 'k')


    def makeParametricPlot(self, u=0.4, mu=0.55):
        '''Makes parametric time plots - Phillips Curve plot should have been made already'''
        self.Z0 = [u, mu]  # initial conditions for u and mu
        self.ts = linspace(0, 100, 175)
        self.terms_to_use = 20
        self.parametricPlot = Plot(self.parent, 5, 4, 100)
        self.populateParametricPlot()
        self.parametricPlot.figure.tight_layout()
        return self.parametricPlot


    def clearParametricPlot(self):
        '''Clears Parametric Plot'''
        self.parametricPlot.clear()


    def populateParametricPlot(self):
        '''Adds all necessary data to parametric plot'''
        self.Zs = odeint(self.dZ_dt, self.Z0, self.ts, args=(self.v, self.theta, self.n, self.alpha, self.beta))
        self.Z2s = odeint(self.dZ2_dt, self.Z0, self.ts, args=(self.v, self.theta, self.n, self.alpha, self.beta, self.k))
        # use optional argument ’args’ to pass parameters to dZ_dt
        # args=(2,0.2,0.2,0.3,0.2)
        self.A = self.calc_A(self.v, self.theta, self.n, self.alpha)

        if self.exponential:
            self.wage_share = self.Z2s[:, 0]  # first column
            self.employment = self.Z2s[:, 1]  # second column
        else:
            self.wage_share = self.Zs[:, 0]  # first column
            self.employment = self.Zs[:, 1]  # second column

        self.parametricPlot.setX(self.ts)
        self.parametricPlot.setX2(self.ts)
        self.parametricPlot.setY(self.wage_share)
        self.parametricPlot.setY2(self.employment)

        self.M0 = self.Z0[1] * self.beta / (self.alpha + self.theta)
        self.U0 = self.Z0[0] / (self.v * (self.alpha + self.theta))
        self.c2 = -log(self.M0 * self.U0 ** self.A) + self.M0 + self.U0

        # Add axes and title
        self.parametricPlot.setTitle("Parametric Plots of Goodwin Model")
        self.parametricPlot.setAxesTitles("Time", "Percentage As Decimal")
        # Add labels
        self.parametricPlot.setLabel("Wage Share (W.S.)")
        self.parametricPlot.setLabel2("Employment Rate (E.R.)")
        # Draw the equilibrium lines where intersection on the plot means inflection point
        self.findEquiAndExtrema()
        self.parametricPlot.makeHLine(y=self.erEqui, c='g', linestyle='--',
                                      label="E.R. Equilibrium")

        self.parametricPlot.makeHLine(y=self.wsEqui, c='y', linestyle='--',
                                      label="W.S. Equilibrium")

        # Draw the boundary lines where intersection on the plot means derivative is 0

        self.parametricPlot.makeHLine(y=self.wsUpperB, c='c', linestyle='--', label="W.S. Boundary")

        self.parametricPlot.makeHLine(y=self.wsLowerB, c='c', linestyle='--')

        self.parametricPlot.makeHLine(y=self.erUpperB, c='r', linestyle='--',
                                      label="E.R. Boundary")

        self.parametricPlot.makeHLine(y=self.erLowerB, c='r', linestyle='--')
        # print("WS Boundaries:", self.wsUpperB,"-",self.wsLowerB)
        # print("ER Boundaries:", self.erUpperB,"-",self.erLowerB)
        self.parametricPlot.plot()


    def findEquiAndExtrema(self):
        '''Finds the equilibrium and extrema values'''
        if self.exponential == True:
            # Find equilibrium points
            self.erEqui = 1 - self.v * (self.theta + self.n)
            self.wsEqui = log((self.alpha + self.theta) / self.beta) / self.k

            # Approximated integral values
            self.A_mu0 = self.approximate_integral_mu(self.k, self.Z0[1], self.terms_to_use)
            self.A_mu = self.approximate_integral_mu(self.k, self.wsEqui, self.terms_to_use)

            # Common factors used in equation
            self.f_mu = (1 / self.v - (self.theta + self.n))
            self.f_mu2 = self.v * (self.theta + self.n) - 1

            self.c = -(self.alpha + self.theta) * log(self.Z0[1]) + self.beta * self.A_mu0 - self.f_mu * \
                     log(self.Z0[0]) + self.Z0[0] / self.v

            coefficients = [self.alpha, self.theta, self.beta, self.k, self.v, self.n, self.c]

            # Find the extrema
            # Bisection method required for employment rate (mu) because integral can't be integrated
            TOL = 1e-5
            mu_lower_bound_tuple = self.bisect_mu_root(coefficients, self.erEqui, 1e-10, self.wsEqui, TOL,
                                                       self.terms_to_use)
            mu_upper_bound_tuple = self.bisect_mu_root(coefficients, self.erEqui, self.wsEqui, 1, TOL, self.terms_to_use)
            self.u_lamb = self.wsEqui ** (-(self.alpha + self.theta) / self.f_mu) * \
                          exp((self.beta * self.A_mu - self.c) / self.f_mu) / self.f_mu2
            # print("Inner Lambert:",u_lamb)
            self.wsUpperB = self.f_mu2 * real(lambertw(self.u_lamb, -1))
            self.wsLowerB = self.f_mu2 * real(lambertw(self.u_lamb))
            self.erUpperB = (mu_upper_bound_tuple[0] + mu_upper_bound_tuple[1]) / 2
            self.erLowerB = (mu_lower_bound_tuple[0] + mu_lower_bound_tuple[1]) / 2
        else:
            # Find equilibrium points
            self.erEqui = 1 - self.v * (self.theta + self.n)
            self.wsEqui = (self.alpha + self.theta) / self.beta

            # Find extrema
            self.wsUpperB = (self.alpha + self.theta) * self.v * -self.A * real(
                lambertw(-exp((1 - self.c2) / self.A) / self.A))
            self.wsLowerB = (self.alpha + self.theta) * self.v * -self.A * real(
                lambertw(-exp((1 - self.c2) / self.A) / self.A, -1))
            self.erUpperB = -(self.alpha + self.theta) / self.beta * real(
                lambertw(-exp((self.A - self.c2)) / (self.A ** self.A)))
            self.erLowerB = -(self.alpha + self.theta) / self.beta * real(
                lambertw(-exp((self.A - self.c2)) / (self.A ** self.A), -1))


    def makeGoodwinPlot(self):
        '''Makes a plot of the Goodwin Model - should be called after parametric plots are made'''
        self.goodwinPlot = Plot(self.parent, 5, 4, 100)
        self.populateGoodwinPlot()
        self.goodwinPlot.figure.tight_layout()
        return self.goodwinPlot


    def clearGoodwinPlot(self):
        '''Clears Goodwin Plot'''
        self.goodwinPlot.clear()


    def populateGoodwinPlot(self):
        '''Adds all necessary data to Goodwin plot'''
        self.goodwinPlot.setX(self.wage_share)
        self.goodwinPlot.setY(self.employment)
        self.goodwinPlot.plot()

        # Add axes and title
        self.goodwinPlot.setTitle("A Cycle of the Goodwin Model")
        self.goodwinPlot.setAxesTitles("Wage Share", "Employment Rate")

        # Draw the equilibrium lines where intersection on the plot means inflection point
        self.goodwinPlot.makeVLine(x=self.erEqui, c='g', linestyle='--',
                                   label="E.R. Equilibrium")

        self.goodwinPlot.makeHLine(y=self.wsEqui, c='y', linestyle='--',
                                   label="W.S. Equilibrium")

        # Draw the boundary lines where intersection on the plot means derivative is 0
        self.goodwinPlot.makeVLine(x=self.wsUpperB, c='c', linestyle='--', label="W.S. Boundary")
        self.goodwinPlot.makeVLine(x=self.wsLowerB, c='c', linestyle='--')
        self.goodwinPlot.makeHLine(y=self.erUpperB, c='r', linestyle='--', label="E.R. Boundary")
        self.goodwinPlot.makeHLine(y=self.erLowerB, c='r', linestyle='--')
