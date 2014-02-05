from __future__ import division
import numpy
import sympy
from sympy.core.expr import Expr

def partial_correlation(x,y,const_list,par_list,correlation_matrix, level=1):
    if len(const_list) == 1:
        z = const_list[0]
        ix = par_list.index(x)
        iy = par_list.index(y)
        iz = par_list.index(z)
        r_xy = correlation_matrix[ix][iy]
        r_xz = correlation_matrix[ix][iz]
        r_yz = correlation_matrix[iy][iz]
    else:
        z = const_list.pop()
        r_xy = partial_correlation(x,y,const_list[:],par_list,correlation_matrix,level=level+1)
        r_xz = partial_correlation(x,z,const_list[:],par_list,correlation_matrix,level=level+1)
        r_yz = partial_correlation(y,z,const_list[:],par_list,correlation_matrix,level=level+1)
    correlation = (r_xy - r_xz * r_yz)/numpy.sqrt((1-r_xz**2)*(1-r_yz**2))
    return correlation
    
def partial_correlation_matrix(par_list, correlation_matrix):   
    import numpy
    size = len(correlation_matrix)
    if size == 2: return correlation_matrix
    pcor = numpy.zeros((size,size))
    for i, row in enumerate(correlation_matrix):
        for j, cor in enumerate(row):
            if i == j:
                pcor[i][j] = 1
            else:
                pvars = [par_list[i],par_list[j]]
                pconst = par_list[:]
                pconst.remove(pvars[0])
                pconst.remove(pvars[1])
                
                pcor[i][j] = partial_correlation(par_list[i],par_list[j], pconst, par_list, correlation_matrix)
    return pcor

### Modify sage.numerical.optimize.find_fit to return parameter errors, covariance matrix, correlation matrix, and partial correlation matrix ###
def find_fit(data, model, initial_guess = None, parameters = None, variables = None, solution_dict = True, partial_correlations=False, weights=None):

    if not isinstance(data, numpy.ndarray):
        try:
            data = numpy.array(data, dtype = float)
        except (ValueError, TypeError):
            raise TypeError, "data has to be a list of lists, a matrix, or a numpy array"
    elif data.dtype == object:
        raise ValueError, "the entries of data have to be of type float"

    if data.ndim != 2:
        raise ValueError, "data has to be a two dimensional table of floating point numbers"

    if isinstance(model, Expr):
        if variables is None:
            variables = list(model.arguments())
        if parameters is None:
            parameters = list(model.variables())
            for v in variables:
                parameters.remove(v)

    if data.shape[1] != len(variables) + 1:
        raise ValueError, "each row of data needs %d entries, only %d entries given" % (len(variables) + 1, data.shape[1])

    if parameters is None or len(parameters) == 0 or \
       variables is None or len(variables) == 0:
        raise ValueError, "no variables given"

    if initial_guess == None:
        initial_guess = len(parameters) * [1]

    if not isinstance(initial_guess, numpy.ndarray):
        try:
            initial_guess = numpy.array(initial_guess, dtype = float)
        except (ValueError, TypeError):
            raise TypeError, "initial_guess has to be a list, tuple, or numpy array"
    elif initial_guess.dtype == object:
        raise ValueError, "the entries of initial_guess have to be of type float"

    if len(initial_guess) != len(parameters):
        raise ValueError, "length of initial_guess does not coincide with the number of parameters"

    if isinstance(model, Expr):
        var_list = variables + parameters
        func = sympy.lambdify(var_list, model)
    else:
        func = model

    def function(x_data, params):
        result = numpy.zeros(len(x_data))
        for row in xrange(len(x_data)):
            fparams = numpy.hstack((x_data[row], params)).tolist()
            result[row] = func(*fparams)
        return result

    def error_function(params, x_data, y_data, weights=weights):
        result = numpy.zeros(len(x_data))
        for row in xrange(len(x_data)):
            fparams = x_data[row].tolist() + params.tolist()
            result[row] = func(*fparams)
        if weights:
            return (result - y_data)/weights
        else: 
            return result - y_data

    x_data = data[:, 0:len(variables)]
    y_data = data[:, -1]

    from scipy.optimize import leastsq
    args = (x_data, y_data)
    output = leastsq(error_function, initial_guess, args = args, full_output=True)
    estimated_params = output[0]
    
    s_sq = (error_function(estimated_params, *args)**2).sum()/(len(y_data)-len(initial_guess))
    if output[1] != None: 
        pcov = output[1]*s_sq

    if isinstance(estimated_params, float):
        estimated_params = [estimated_params]
    else:
        estimated_params = estimated_params.tolist()

    if solution_dict:
        solution_dict = {}
        for item in zip(parameters, estimated_params):
            solution_dict[item[0]] = item[1]
        stdev_dict = {}
        if output[1]!=None:
            for item in zip(parameters, pcov.diagonal()):
                stdev_dict[item[0]] = numpy.sqrt(item[1])
           
            size=len(pcov)
            pcor = numpy.zeros((size,size))
            for i, row in enumerate(pcov):
                std_i = numpy.sqrt(row[i])
                for j, cov in enumerate(row):
                    std_j = numpy.sqrt(pcov[j][j])
                    pcor[i][j] = cov/(std_i*std_j)
                    
            if partial_correlations:
                partial_pcor = partial_correlation_matrix(parameters,pcor)
            else: partial_pcor = None        
 
        else:
            stdev_dict = 'Unavailable (singular matrix encountered).'
            pcov = 'Unavailable (singular matrix encountered).'
            pcor = 'Unavailable (singular matrix encountered).'
            partial_pcor = 'Unavailable (singular matrix encountered).'
        return solution_dict, stdev_dict, pcov, pcor, partial_pcor

    else:
        return [item[0] == item[1] for item in zip(parameters, estimated_params)]