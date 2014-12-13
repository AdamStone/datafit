datafit.py
=========

This is a python implementation for fitting experimental data to 
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

The original code was written in the Sage mathematics environment (http://www.sagemath.org) during a graduate data analysis class in order to analyze a series of specific case studies, which are included in the repository along with some of the code adapted to a pure python implementation. Because it began as Sage worksheets, it retains some aspects of the Sage implementation such as the use of symbolic expressions to define model functions (via the sympy module).

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




Dependencies
-------------

numpy, scipy, matplotlib, and sympy.
