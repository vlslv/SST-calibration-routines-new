#!/usr/bin/python
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array, pi
import numpy as np
import numpy.linalg as la

from scipy import optimize
from matplotlib import pyplot as plt, cm, colors

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

if __name__=="__main__":
    from matplotlib.pyplot import plot, scatter, show
    series = [0,0,0,2,0,0,0,-2,0,0,0,2,0,0,0,-2,0]
    maxtab, mintab = peakdet(series,.3)
    plot(series)
    scatter(array(maxtab)[:,0], array(maxtab)[:,1], color='blue')
    scatter(array(mintab)[:,0], array(mintab)[:,1], color='red')
    show()


from scipy.odr import odrpack as odr
from scipy.odr import models

def poly_lsq(x,y,n,verbose=False,itmax=200):
    ''' Performs a polynomial least squares fit to the data,
    with errors! Uses scipy odrpack, but for least squares.

    IN:
       x,y (arrays) - data to fit
       n (int)      - polinomial order
       verbose      - can be 0,1,2 for different levels of output
                      (False or True are the same as 0 or 1)
       itmax (int)  - optional maximum number of iterations

    OUT:
       coeff -  polynomial coefficients, lowest order first
       err   - standard error (1-sigma) on the coefficients

    --Tiago, 20071114i
    https://github.com/tiagopereira/python_tips/wiki/Scipy:-curve-fitting
    '''

    # http://www.scipy.org/doc/api_docs/SciPy.odr.odrpack.html
    # see models.py and use ready made models!!!!

    func   = models.polynomial(n)
    mydata = odr.Data(x, y)
    myodr  = odr.ODR(mydata, func,maxit=itmax)

    # Set type of fit to least-squares:
    myodr.set_job(fit_type=2)
    if verbose == 2: myodr.set_iprint(final=2)

    fit = myodr.run()

    # Display results:
    if verbose: fit.pprint()

    if fit.stopreason[0] == 'Iteration limit reached':
        print '(WWW) poly_lsq: Iteration limit reached, result not reliable!'

    # Results and errors
    coeff = fit.beta[::-1]
    err   = fit.sd_beta[::-1]

    return coeff,err

##################################################################################################
### Nice work of Luciano Riano
### https://gist.github.com/lorenzoriano/6799568
##################################################################################################

def calc_R(x,y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f(c, x, y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()

def leastsq_circle(x,y):
    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
    xc, yc = center
    Ri       = calc_R(x, y, *center)
    R        = Ri.mean()
    residu   = np.sum((Ri - R)**2)
    return xc, yc, R, residu

def plot_data_circle(x,y, xc, yc, R,leyenda=True):
    '''
    Version modificada, entrada opcional de leyenda
    '''
    f = plt.figure(facecolor='white',figsize=(7, 5.4), dpi=72)
    plt.axis('equal')

    theta_fit = np.linspace(-pi, pi, 180)

    x_fit = xc + R*np.cos(theta_fit)
    y_fit = yc + R*np.sin(theta_fit)
    plt.plot(x_fit, y_fit, 'b-' , label="fitted circle", lw=2)
    plt.plot([xc], [yc], 'bD', mec='y', mew=1)
    plt.xlabel('x')
    plt.ylabel('y')   
    # plot data
    plt.plot(x, y, 'ro', label='data', mew=1)
    if leyenda==True:
	 plt.legend(loc='best',labelspacing=0.1,frameon=False)
    plt.grid()
    plt.title('Least Squares Circle')

##############################################################################
# The angle between two vectors, Python version
#https://newtonexcelbach.wordpress.com/2014/03/01/the-angle-between-two-vectors-python-version/
##############################################################################
 
#@xl_func("numpy_row v1, numpy_row v2: float")
def py_ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)
