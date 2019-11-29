#!/usr/bin/env python

"""
Reading in and manipulating fits files

"""

import numpy as np
import pandas as pd
from scipy import optimize
from scipy.interpolate import interpn

def read_exopars(dataframe,par1='a',par2='Mstar',meth='RV'):
    M_earth = 5.97
    M_jupiter = 1898
    method = dataframe[:,1]
    filter = 'Radial Velocity'

    # Return only data acquired with input method
    if meth == 'RV':
        filter = 'Radial Velocity'
    elif meth == 'Transit':
        filter = 'Transit'
    elif meth == 'PT':
        filter = 'Pulsar Timing'
    elif meth == 'ML':
        filter = 'Microlensing'
    elif meth == 'Im':
        filter = 'Imaging'
    
    # Read in data from input dataframe, using filter on method
    a = pd.to_numeric(dataframe[:,7][method==filter])
    a_perr = pd.to_numeric(dataframe[:,8][method==filter])
    a_nerr = pd.to_numeric(dataframe[:,9][method==filter])
    a_err = [-1*a_nerr,a_perr]
    Mstar = pd.to_numeric(dataframe[:,19][method==filter]) * (M_jupiter/M_earth)
    Mstar_perr = pd.to_numeric(dataframe[:,20][method==filter]) * (M_jupiter/M_earth)
    Mstar_nerr = pd.to_numeric(dataframe[:,21][method==filter]) * (M_jupiter/M_earth)
    Mstar_err = [-1*Mstar_nerr,Mstar_perr]

    return a,a_err,Mstar,Mstar_err


def linfit(xdata, ydata, yerr=None, pinit=[1.0,-1.0]):

    #logx = np.log10(xdata)
    #logy = np.log10(ydata)
    #logyerr = yerr/ydata
    pinit=[5.8,0.3]
    yerr = np.full(len(xdata),0.1)
    
    # Define function for calculating a power law
    linfit = lambda p, x: p[0] + p[1] * x
    linerr = lambda p, x, y, err: (y-linfit(p,x))/err

    # Fit data with function defined above
    out = optimize.leastsq(linerr,pinit,args=(xdata,ydata,yerr[0]), full_output=1)

    # Determine best-fit parameters and associated errors
    pfinal = out[0]
    covar = out[1]

    index = pfinal[1]
    intercep = pfinal[0]
    indexErr = np.sqrt( covar[1][1] )
    intercepErr = np.sqrt( covar[0][0] )
    
    return index,indexErr,intercep,intercepErr

def Zrecal(R23_met):
    a = 664.8453
    b=-225.7533
    c=25.76888
    d=-0.9761368
    
    O3N2_met = a + (b*x) + (c * x**2) + (d * x**3)
    
    return O3N2_met

def conflevels(x,y,nbins,confints=[0.99,0.95,0.68]):
    # Make a 2d normed histogram
    H,xedges,yedges=np.histogram2d(x,y,bins=nbins,normed=True)

    norm=H.sum() # Find the norm of the sum
    # Set contour levels
    contour1=0.99
    contour2=0.95
    contour3=0.68

    # Set target levels as percentage of norm
    target1 = norm*contour1
    target2 = norm*contour2
    target3 = norm*contour3

    # Take histogram bin membership as proportional to Likelihood
    # This is true when data comes from a Markovian process
    def objective(limit, target):
        w = np.where(H>limit)
        count = H[w]
        return count.sum() - target

    # Find levels by summing histogram to objective
    level1= optimize.bisect(objective, H.min(), H.max(), args=(target1,))
    level2= optimize.bisect(objective, H.min(), H.max(), args=(target2,))
    level3= optimize.bisect(objective, H.min(), H.max(), args=(target3,))

    # For nice contour shading with seaborn, define top level
    level4=H.max()
    levels=[level1,level2,level3,level4]

    return levels

def density_scatter(x , y, ax=None, sort=True, bins=20,):
    """
    Scatter plot colored by 2d histogram
    """
    #if ax is None :
    #    fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    #ax.scatter( x, y, c=z, **kwargs )
    return z
