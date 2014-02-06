from __future__ import division

from datafit import get_data, get_xy, Model
from sympy import symbols
    
    
# Initialize data
raw_data = get_data('data/case1/Refractive Index vs Temperature.csv', delimiter='\t')
x, y = get_xy(raw_data, 0, 1, plot=False)



#####################################
""" TRY FITTING WITH 1 OSCILLATOR """
#####################################

# Set up model
S1, L1 = symbols('S1, L1', domain='positive')
X = symbols('X')

osc1 = S1/(1 - L1**2/X**2)   # expression for 1 oscillator
n = (osc1 + 1)**(1/2)        # Sellmeier eq. with 1 oscillator

model1 = Model(n, name='Sellmeier (1 oscillator)')


# Perform fit
guess = {S1: 1, L1: 200}

fit = model1.find_fit(x, y, guess)

if fit:
    # if a fit was successful...
    fit.summary()
    fit.plot()
    fit.plot_residuals()

#raw_input( 'Press return to continue...' )
   







######################################
""" TRY FITTING WITH 2 OSCILLATORS """
######################################

# Set up model
S2, L2 = symbols('S2, L2', domain='positive') # new parameters
    
osc2 = S2/(1 - L2**2/X**2)     # new oscillator
n = (osc1 + osc2 + 1)**(1/2)   # Sellmeier eq. with 2 oscillators

model2 = Model(n, name='Sellmeier (2 oscillators)')


# Perform fit
guess2 = fit.results.copy()      # use previous result in new guess
guess2.update({S2:1, L2:5000})   # guess low frequency oscillator

fit2 = model2.find_fit(x, y, guess2, partial_correlations=True)

if fit2:
    # if a fit was successful...
    fit2.summary()
    fit2.plot()
    fit2.plot_residuals()

#raw_input( 'Press return to continue...' )
    
    
    
    
    
    
    
    
######################################
""" TRY FITTING WITH 3 OSCILLATORS """
######################################

# set up model
S3, L3 = symbols('S3, L3', domain='positive') # new parameters

osc3 = S3/(1 - L3**2/X**2)                    # new oscillator
n = (osc1 + osc2 + osc3 + 1)**(1/2)           # Sellmeier eq. with 3 oscillators

model3 = Model(n, name='Sellmeier (3 oscillators)')


# Perform fit
guess3 = fit2.results.copy()
guess3.update({S3:1, L3:250})   # guess higher frequency oscillator

fit3 = model3.find_fit(x, y, guess3, partial_correlations=True)

if fit3:
    # if a fit was successful...
    fit3.summary()
    fit3.plot()
    fit3.plot_residuals()
    
#raw_input( 'Press return to continue...' )
      
    
    
    
    
    
    
    
#####################################
""" APPROXIMATE FAR-IR OSCILLATOR """
#####################################

# set up model
D = symbols('D', domain='positive')      # new parameter

oscIR = -D*X**2                          # IR oscillator approximation
n = (osc1 + oscIR + osc3 + 1)**(1/2)     # Sellmeier eq. with 2.5 oscillators

model4 = Model(n, name='Sellmeier (2 oscillators + IR approximation)')


# Perform fit
guess4 = fit3.results.copy()
guess4.update({D:1e-8})

fit4 = model4.find_fit(x, y, guess4, partial_correlations=True)

if fit4:
    # if a fit was successful...
    fit4.summary()
    fit4.plot()
    fit4.plot_residuals()