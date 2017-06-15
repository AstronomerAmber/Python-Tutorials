#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:07:08 2017

@author: cburns

Fit the sky background at LCO using CSP I and II data
"""

from sqlalchemy import create_engine
import pandas as pd
import ephem
from numpy import power,absolute,cos,sin,pi,exp,log10,array
from scipy.optimize import leastsq
import sys,os,string

# I like putting all my functions at the top of the script and save the main
# script for below
def X(Z):
    # return the airmass as a function of zenith angle in radians
    return power(1-0.96*power(sin(Z),2),-0.5)

def modelsky(alpha, rho, kx, Z, Zm, mdark):
   # Given the variables, compute the sky brightness in mag/square-arc-sec
   # due to the moon
   Istar = power(10, -0.4*(3.84+0.026*absolute(alpha)+4e-9*power(alpha,4)))
   frho = power(10, 5.36)*(1.06 + power(cos(rho),2))+power(10, 6.15-rho*180./pi/40)
   Bmoon = frho*Istar*power(10,-0.4*kx*X(Zm))*(1-power(10,-0.4*kx*X(Z)))
   Bdark = 34.08*exp(20.723 - 0.92104*mdark)*power(10,-0.4*kx*(X(Z)-1))*X(Z)
   return mdark - 2.5*log10((Bmoon+Bdark)/Bdark)

def fitfunc(p, alpha, rho, Z, Zm, magsky):
    # function needed by scipiy.optimize.leastsq
    mdark,kx = p
    return magsky - modelsky(alpha, rho, kx, Z, Zm, mdark)

# __MAIN__

# Which filter do we want to fit?
filt = sys.argv[1]
engine = create_engine('mysql+pymysql://bootcamp:pmactoob@kepler.obs.carnegiescience.edu/Phot')

query = '''
select MAGINS.jd,MAGINS.alpha as RA, MAGINS.delta as Decl, MAGINS.airm,
       -2.5*log10(MAGINS.sky / 0.435/0.435) + 2.5*log10(MAGINS.expt) + MAGFIT1.zp + 25 as magsky
from (MAGINS,MAGFIT1)
where MAGINS.night=MAGFIT1.night and MAGINS.filt=MAGFIT1.filt
      and MAGINS.filt='{}' and MAGINS.field like 'SA%%' and MAGINS.sky > 0'''
      
# do the query and get the data as a Pandas dataFrame
data = pd.read_sql_query(query.format(filt), engine)

# Now we need a loop to compute all the angles using ephem
OBS = ephem.Observer()
OBS.long = "-70.6926"      # longitude (negative --> West)
OBS.lat = "-20.0146"       # latitude (negative --> South)
OBS.elev = 2380            # elevation in meters
JDoff = 2415020

sky = ephem.FixedBody()
moon = ephem.Moon()
sun = ephem.Sun()
alpha = []
rho = []
Z = []
Zm = []
sun = ephem.Sun()
for JD,RA,DEC in zip(data['jd'],data['RA'],data['Decl']):
    sky._ra = RA*pi/180    # in radians
    sky._dec = DEC*pi/180
    sky._epoch = ephem.J2000
    OBS.date = JD-JDoff
    sky.compute(OBS)
    moon.compute(OBS)
    sun.compute(OBS)
    alpha.append((pi - ephem.separation(moon,sun))*180./pi)   # in degrees
    rho.append(ephem.separation(moon,sky))                    # in radians
    Z.append(pi/2 - sky.alt)                              # in radians
    Zm.append(pi/2 - moon.alt)                            # in radians
    
data['alpha'] = pd.Series(array(alpha), index=data.index)
data['rho'] = pd.Series(array(rho), index=data.index)
data['Z'] = pd.Series(array(Z), index=data.index)     # radians
data['Zm'] = pd.Series(array(Zm), index=data.index)

# Do the least-squares fit.
pars,stat = leastsq(fitfunc, [22,0.2], args=(data['alpha'],data['rho'],data['Z'],
                    data['Zm'],data['magsky']))
print pars
data['modelsky']=pd.Series(modelsky(data['alpha'],data['rho'],pars[1],
                           data['Z'],data['Zm'],pars[0]), index=data.index)

# Visualize!
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.io import show,output_file
from bokeh.models import ColumnDataSource

output_file('fitsky_%s.html' % filt)
source = ColumnDataSource(data)
TOOLS = ['box_select','lasso_select','reset','box_zoom','help']
vars = [('modelsky','magsky'),('alpha','rho'),('alpha','Zm'),
        ('jd','alpha'),('Z','Zm'),('RA','Decl')]
plots = []
for var in vars:
   s = figure(tools=TOOLS, plot_width=300, plot_height=300)
   s.circle(*var, source=source, selection_color='red', alpha=0.1)
   s.xaxis.axis_label = var[0]
   s.yaxis.axis_label = var[1]
   plots.append(s)
plots[0].line([17.8,22.3],[17.8,22.3], line_color='orangered')

p = gridplot([plots[0:3],plots[3:]])
show(p)
