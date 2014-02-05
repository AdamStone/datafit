from __future__ import division
import matplotlib.pyplot as plt
import numpy
import sympy
from sympy import Symbol

from optimize import find_fit


# TODO setup clean trial-and-error and plotting of guess values
# TODO account for non-symbolic (python function) model expressions

def get_xy(csv_data, x_col, y_col, plot=False):
    """ Given csv data, returns (x, y) tuples from the specified x and y columns
    (excluding the first row).
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
    def __init__(self, model, data, parameters, results, stdev, covariance, 
                 correlations, partial_correlations):
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
        print 'Fit results:\n', self.results, '\n'
        print 'Sum of squares of residuals:\n', self.ssq, '\n'
        print 'Parameter standard deviation:\n', self.stdev, '\n'
        print 'Covariance matrix:\n', self.parameters, '\n', self.covariance, '\n'
        print 'Correlation matrix:\n', self.parameters, '\n', self.correlations, '\n'
        if self.partial_correlations != None:
            print 'Partial correlation matrix:\n', self.parameters, '\n', self.partial_correlations, '\n'
    
    def plot(self):
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
        
        
        
        

class Model2D(object):
    def __init__(self, expression, name=None):
        self.expression = expression
        self.name = name
        self.variables = expression.atoms(Symbol)

    def function(self, x, parameter_values):
        """ parameter_values is a dict of the form {symbol:value} for all 
        parameters (excluding the x variable). """

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