datafit.py
=========

This python package is a framework for fitting experimental data to 
mathematical models and analyzing the quality of fit. The fitting
itself is performed by a Levenberg-Marquardt algorithm via the
scipy.optimize.leastsq function, and the main role of this code is
to facilitate viewing and analysis of the fitting results, including

    Standard deviations of best-fit parameter values
    Sum of squares of residuals
    Parameter covariance matrix
    Parameter correlation matrix
    Parameter partial correlation matrix
    Plotting of model against data
    Plotting of fit residuals

This code was developed over the course of a graduate data analysis
class in order to analyze a series of specific case studies, which
are included as examples to explain the various outputs of the fits
and illustrate how they can be interpreted for real data. 

The original code was written in the sage mathematics environment 
(http://www.sagemath.org), and was later adapted to a pure python 
implementation to make it more easily accessible. As such, it retains 
some characteristics of the original sage implementation, including 
the use of symbolic expressions to define model functions (via the 
sympy module).


Dependencies
-------------

numpy, scipy, matplotlib, and sympy.


Getting Started
-------------

The included case study scripts illustrate the basic use of the
code, while the corresponding PDFs discuss the background theory and 
interpretation of fitting results. 

In general, csv data is imported and the x and y columns selected using 
the get_data and get_xy functions. The model expression is defined 
symbolically, using the sympy.symbols method to initialize the variables 
and fitting parameters, and passed to a Model object. The x and y 
data are then passed to the find_fit method of the model along with a 
dictionary of initial guess values. If the fit is successful, an instance
of the Fit class will be returned. 

The Fit class includes summary, plot, and plot_residuals methods. The
summary prints the best-fit parameter values as well as several outputs 
useful for analyzing the quality of fit. By default the partial correlation
matrix is omitted, as it can become relatively slow to calculate as the
number of parameters increases, but it can be requested by setting the
partial_correlations keyword of the Model's find_fit call to True. 
