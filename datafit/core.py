from __future__ import division
import matplotlib.pyplot as plt
import numpy
import sympy
import csv
from sympy import Symbol

from optimize import find_fit


# TODO setup clean trial-and-error and plotting of guess values
# TODO account for non-symbolic (python function) model expressions


def get_data(csv_filename, delimiter=','):
    """ Shortcut for importing csv data from file.
    
    Parameters
    ----------
    
    csv_filename : string
        Location/filename of the data file (including extension).
        
    delimiter : string
        Specifies how columns are separated (usually ',' or '\t')
        
    Returns
    ----------
    
    raw_data : list
        CSV data as a list of all rows in the data file.
    """
    with open(csv_filename, 'rb') as csvfile:
        raw_data = [row for row in csv.reader(csvfile, delimiter=delimiter)]
        return raw_data

def get_xy(csv_data, x_col, y_col, plot=False):
    """ Given CSV data, returns numpy arrays of the specified x and y columns
    (excluding the first row of labels).
    
    Parameters
    ----------
    
    csv_data : list
        CSV data as a list of all rows in the data file.
        
    x_col : int
        Integer index of column of x data.
        
    y_col : int
        Integer index of column of y data.
        
    plot : bool
        If True, execution will be paused while a plot of the data is shown.
        
    Returns
    ----------
    
    x : `numpy.ndarray`
        Numpy 1D array containing x values.
        
    y : `numpy.ndarray`
        Numpy 1D array containing y values.
    """
    x = numpy.array([float(row[x_col]) for row in csv_data[1:] if row[y_col] != ''])
    y = numpy.array([float(row[y_col]) for row in csv_data[1:] if row[y_col] != ''])

    if plot:
        plt.clf()
        plt.scatter(x, y)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Refractive Index')
        print 'Close plot to continue...'
        plt.show()
    return x, y






class Fit(object):
    """ Class for storing and displaying fit results. """
    def __init__(self, model, data, parameters, results, stdev, covariance, 
                 correlations, partial_correlations):
        """ 
        Parameters
        ----------
        
        model : `Model`
            The Model instance with which the fit was obtained.
            
        data : list
            List of (x, y) tuples representing experimental data, e.g. 
            obtained by zip(x, y).
            
        parameters : list
            List of model parameters to be optimized by the fit. 
            
        results : dict
            Dictionary mapping parameters to best-fit values (the first
            item returned by optimize.find_fit)
        
        stdev : dict
            Dictionary mapping parameters to best-fit standard deviations
            (the second item returned by optimize.find_fit)
            
        covariance : `numpy.ndarray`
            Array representing the covariance matrix corresponding to the
            best-fit results (the third item returned by optimize.find_fit).
            
        correlations : `numpy.ndarray`
            Array representing the correlation matrix corresponding to the
            best-fit results (the fourth item returned by optimize.find_fit).
            
        partial_correlations: `numpy.ndarray` or None
            Array representing the partial correlation matrix corresponding
            to the best-fit results (the fifth item returned by optimize.find_fit).
        """
        self.model = model
        self.data = data
        self.parameters = parameters
        self.results = results
        self.stdev = stdev
        self.covariance = covariance
        self.correlations = correlations
        self.partial_correlations = partial_correlations

        self.x, self.y_data = numpy.array(zip(*self.data))
        self.y_model = self.model.function(self.x, self.results)        
        self.residuals = self.y_data - self.y_model
        self.ssq = sum(self.residuals**2)
        
    def summary(self):
        """ Prints a summary of results of the fit. Matrices are printed with
        the order of parameters corresponding to each column indicated above, 
        and the transpose should be taken to correspond to the order of rows.
        """
        print 'Fit results:\n', self.results, '\n'
        print 'Parameter standard deviation:\n', self.stdev, '\n'
        print 'Sum of squares of residuals:\n', self.ssq, '\n'
        print 'Covariance matrix:\n', self.parameters, '\n', self.covariance, '\n'
        print 'Correlation matrix:\n', self.parameters, '\n', self.correlations, '\n'
        if self.partial_correlations != None:
            print 'Partial correlation matrix:\n', self.parameters, '\n', self.partial_correlations, '\n'
    
    def plot(self):
        """ Plots the model fit as a curve against the original scatter data. 
        """
        plt.clf()
        x, y_data, y_model = self.x, self.y_data, self.y_model
        
        xpad = (max(x) - min(x)) * 0.05
        xmin = min(x) - xpad
        xmax = max(x) + xpad
        
        ymin = min(min(y_data), min(y_model))
        ymax = max(max(y_data), max(y_model))
        ypad = (ymax - ymin) * 0.05
        ymin -= ypad
        ymax += ypad
                
        plt.scatter(x, y_data, facecolor='blue', edgecolor='black')
        plt.plot(x, y_model)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Refractive Index')
        plt.title('Best fit results')
        plt.axis([xmin, xmax, ymin, ymax])
        print 'Close plot to continue...'
        plt.show()

    def plot_residuals(self):
        """ Plots the residuals (difference between data and model at each
        point), assuming errors in y.
        """
        plt.clf()
        x, residuals = self.x, self.residuals
        
        xpad = (max(x) - min(x)) * 0.05
        xmin = min(x) - xpad
        xmax = max(x) + xpad
        
        ypad = (max(residuals) - min(residuals)) * 0.05
        ymin = min(residuals) - ypad
        ymax = max(residuals) + ypad

        plt.hlines(0, xmin, xmax)
        plt.scatter(x, residuals, facecolor='blue', edgecolor='black')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Residuals')
        plt.title('Best fit residuals')
        plt.axis([xmin, xmax, ymin, ymax])
        print 'Close plot to continue...'
        plt.show()
        
        
        
        

class Model(object):
    """ Class for storing model expression and performing fit against
    provided data. """
    def __init__(self, expression, name=None):
        """ 
        Parameters
        ----------
        
        expression : `sympy.Expr`
            Symbolic expression representing the mathematical model to which
            the data should be fit.
            (see http://docs.sympy.org/dev/gotchas.html#variables)
            
        name : string
            Name of the model.
        
        """
        self.expression = expression
        self.name = name
        self.variables = expression.atoms(Symbol)

    def function(self, x, parameter_values):
        """ Generates a python function that accepts a single argument (x) and
        returns an evaluation of the model expression at the value of x 
        for the given set of parameter values. 
        
        Parameters
        ----------
        
        x : float
            The value (or array of values) of x at which the model should
            be evaluated.
            
        parameter_values : dict
            Dictionary mapping the model parameters (the sympy.Symbol objects
            used in creating the symbolic expression) to the values to 
            which they should be set (e.g. initial guess values, or 
            values returned by the fit).
        
        Returns
        ----------
        
        ret : function
            The model expression as a (python) function of x with the provided
            parameter values assigned.
        """
        
        parameters = []
        variables = []
        args = [x]

        for var in self.variables:
            if var in parameter_values:
                parameters.append(var)
                args.append(parameter_values[var])
            else:
                variables.append(var)
                
        return sympy.lambdify(variables + parameters, self.expression)(*args)
        
    def find_fit(self, x, y, initial_guess, partial_correlations=False):
        """ Method which arranges the variables, parameters, and guesses 
        in a consistent order and passes them to the find_fit function.
        
        Parameters
        ----------

        x : float
            Array of data x coordinates.

        y : float
            Array of data y coordinates.
        
        initial_guess : dict
            Dictionary mapping parameter Symbols to guess values. Guesses 
            should be physically reasonable, and some trial-and-error
            is often necessary to find one that produces a successful fit. 

        partial_correlations : bool
            Set True for fit to return a partial correlation matrix. Default
            is false since this can be slow for large models, and since
            this is redundant with the standard correlation matrix for models
            with less than three parameters.
        
        Returns
        ----------

        ret : `Fit` or None
            If successful, returs a Fit instance containing details of 
            the fit results.
                    
        """
        output = 'Attempting fit with model {}'.format(self.name)
        print '#'*len(output)
        print output
        print '#'*len(output)
        
        parameters = []
        variables = []
        guess = []        
        for var in self.variables:
            if var in initial_guess:
                parameters.append(var)
                guess.append(initial_guess[var])
            else:
                variables.append(var)
        
        data = zip(x, y)
        try:
            fit = find_fit(data, self.expression, guess, parameters, 
                   variables, partial_correlations=partial_correlations)
        except:
            print '\nFit failed; guess parameters may be unphysical'
            return
        return Fit(self, data, parameters, *fit)