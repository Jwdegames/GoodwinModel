# GoodwinModel
This Program simulates the relationship between wage share and employment in an economy. Learn more about the simulator on the [Guide](https://github.com/Jwdegames/GoodwinModel/blob/main/Goodwin%20Model%20Simulator%20Guide.pdf).

To use the simulator, you can either download the exe at [the Releases Page](https://github.com/Jwdegames/GoodwinModel/releases), or you can clone this repo and run the Python files using the Python command. Make sure to have PyQt6, PyQt6-WebEngine, numpy, matplotlib, and scipy installed.

This was written in Python and utilizes advanced Calculus and Algebra concepts such as McLaurin Series, Differential Equations, and Lambert Functions. The simulator has a model based off of a linear and exponential Phillips Curve. The solutions to the linear model were based off of here: [Analytic Physics - Interactive Lokta-Volterra Equations](https://analyticphysics.com/Differential%20Equations/Interactive%20Lotka-Volterra%20Equations.htm). McLaurin Series were needed for the exponential version because it is not possible to integrate. 

This simulator shows the equilibrium and extrema, along with 2 variations - a linear Phillips curve, and an exponential Phillips curve. The Goodwin model is based off the Keynesian Harrod-Domar and Keynesian Phillips Curve models, and demonstrates how Employment and Wages creates a predator-prey cycle. This program solves for the equilibrium and extrema using the Lokta-Volterra Model equations. See more info on it's wikipedia page: https://en.wikipedia.org/wiki/Goodwin_model_(economics)

