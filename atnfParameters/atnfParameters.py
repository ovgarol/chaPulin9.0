#!/usr/bin/env python
""" Aqui voy a calcular los parametros para un modelo basado
    en las estadisticas de ATNF.
"""

###########################################################################################
# http://www.scipy.org/Cookbook/OptimizationAndFitDemo1?highlight=%28fit%29
# NO ES LO MISMO EL TIPO LIST Y EL TIPO ARRAY-> ARRAY ES NUMERICO!!!!
# https://www.cfa.harvard.edu/~jbattat/computer/python/science/
# http://matplotlib.org/api/pyplot_api.html#module-matplotlib.pyplot
#
# LOS DATOS SON TOMADOS DE
# http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?P0=P0&W50=W50&R_lum=R_lum&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=p0%3E0.05&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Short+without+errors&no_value=nan&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=40&table_bottom.y=39
###########################################################################################

import numpy as np
import scipy as sp
#import time as tm
import matplotlib.pyplot as plt
from scipy import optimize, array, stats
from matplotlib import rc, rcParams

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})

###########################################################################################
#                                DEFINICIONES PRELIMINARES
###########################################################################################

def dist_lorimer(x,A,B,C):
    return A*np.exp(-(x-B)**2/(2*C**2))     #BASE NATURAL
    #return A*10**(-(x-B)**2/(2*C**2))      #BASE DECIMAL

def dist_norm(x,media,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-media)**2/(2*sigma**2))

def log_Brillo(x,D,E):
    return D*x + E

def log_W(x,F,G):
    #return np.log(F*np.exp(x)**G)
    return G*x + F    

###########################################################################################
#                                       MAIN
###########################################################################################
if __name__ == "__main__":
    np.random.seed()
    mask = 1e6

###########################################################################################
#                              Getting data from file
###########################################################################################

    f = open('atnf-data-clear-pwlq-2.dat', 'r') # open file to read
    logP = []
    logL = []
    logW = []
    logQ = [] #ESTO PARA INCLUIR P_PUNTO EN EL ANalisis p_punto es Q !!!! 
    oP = []
    oQ = []

# ENCABEZADO DE atnf-data-pwlq.dat
# P0 P1 W50 R_LUM  
# (s) (adim) (ms) (mJy kpc^2) 

    for line in f:      # iterate over the lines in the file

        columns = line.split(' ')       # split the line into a list of column values
        columns = [col.strip() for col in columns]  # clean any whitespace off the items

        if columns[1]!= 'NAN':logP.append(np.log(float(columns[1])))      # [P]-> s
        else: logP.append(mask)        

        if columns[3]!= 'NAN' and float(columns[3])!= 0.0:logW.append(np.log(1e-3*1.06*float(columns[3]))) # [W]-> s
        else: logW.append(mask)

        if columns[4]!= 'NAN' and float(columns[4])!= 0.0:logL.append(np.log(1e-3*float(columns[4]))) # [L]-> Jy kpc**2
        else: logL.append(mask)

        if columns[2]!= 'NAN' and float(columns[2])> 0.0:logQ.append(np.log(float(columns[2])))       # [Q]-> s/s (adim)
        else: logQ.append(mask)

        if columns[1]!= 'NAN' and columns[2]!= 'NAN' and float(columns[2])> 0.0 and float(columns[1]) > 0.05: 
            oP.append(np.log(float(columns[1])))     
            oQ.append(np.log(float(columns[2]))) 
    
    f.close()   
    
    logP = array(logP)
    logL = array(logL)
    logW = array(logW)
    logQ = array(logQ)
    
    print 'Pulsares entre ', np.exp(min(logP)), ' y ', np.exp(max(logP)), '[s]'
 
    f.close()

###########################################################################################
#                              Fitting Lorimer's Model
###########################################################################################    

    logP1 = [] # Escoger valores no nulos de P
    for i in range(len(logP)):    
        if logP[i] != mask: 
            logP1.append(logP[i])
    logP1 = array(logP1)
             
    k = 50  # k-> numero de bines del histograma 
    xi = [] # marca de clase del histograma x_i
    a = np.histogram(logP1, bins=k, range=None, normed=False, weights=None, density=None)

    for i in range(len(a[0])): xi.append((a[1][i]+a[1][i+1])/2.)      
    xi = array(xi)
  
    fitfunc = lambda p, x: dist_lorimer(x,p[0],p[1],p[2]) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y     # Distance to the target function
    p0 = [100.0,2.10,0.3]  # (A,B,C) Initial guess for the parameters

    p1, success = optimize.leastsq(errfunc, p0[:], args=(xi, a[0]))

    print 'M. Lorimer -> [A,B,C] ', success, p1, len(logP1), 'datos'

###########################################################################################
#                               Fitting Scaled Gaussian Model
###########################################################################################

    fitfunc = lambda p, x: p[0]*dist_norm(x,p[1],p[2]) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y     # Distance to the target function
    p0 = [100.0,2.10,0.3]  # (A,B,C) Initial guess for the parameters

    p2, success = optimize.leastsq(errfunc, p0[:], args=(xi, a[0]))

    print 'Gaussiana -> [A*,B,C] ', success, p2,  len(logP1), 'datos'

###########################################################################################
#                            Fitting Scaled Gaussian Model for P1
###########################################################################################

    logQ1 = [] # Escoger valores no nulos de Q
    for i in range(len(logQ)):    
        if logQ[i] != mask: 
            logQ1.append(logQ[i])
    logQ1 = array(logQ1)
             
    k = 50  # k-> numero de bines del histograma 
    xi = [] # marca de clase del histograma x_i
    a = np.histogram(logQ1, bins=k, range=None, normed=False, weights=None, density=None)

    for i in range(len(a[0])): xi.append((a[1][i]+a[1][i+1])/2.)      
    xi = array(xi)
  
    fitfunc = lambda p, x: p[0]*dist_norm(x,p[1],p[2]) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y     # Distance to the target function
    p0 = [500.0,-31.0,2.0]  # (A,B,C) Initial guess for the parameters

    p1_1, success = optimize.leastsq(errfunc, p0[:], args=(xi, a[0]))

    print 'M. G. en p1 -> [H*,I,J] ', success, p1_1, len(logQ1), 'datos'


###########################################################################################
#                          Fitting Rubio-Herrera L(logP) Model
###########################################################################################
    
    logP2 = [] # Escoger valores no nulos de P simultaneamente 
    logL2 = [] # con valores no nulos de L, para poder trabajar
    for i in range(len(logP)): 
        if (logP[i] != mask and logL[i] != mask):          
            logP2.append(logP[i])
            logL2.append(logL[i])               
    logP2 = array(logP2)
    logL2 = array(logL2)
    if len(logP2) != len(logL2): print 'AQUI HAY UN GRAVE ERROR'

    fitfunc = lambda p, x: log_Brillo(x,p[0],p[1])  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y     # Distance to the target function
    p0 = [-0.310,-2.10]  # (D,E) Initial guess for the parameters

    p3, success = optimize.leastsq(errfunc, p0[:], args=(logP2, logL2))

    # Fitting the random parameters of the model
    dev = np.sqrt((logL2-log_Brillo(logP2,p3[0],p3[1])).var())  # El mejor estimador de sigma es sqrt(varianza)
    p3 = list(p3)    
    p3.append(dev)
    p3 = array(p3)

    print 'M. R-H L(logP) -> [D,E,sigma] ', success, p3, len(logL2), 'datos'

    # Calcula L para todos los pulsares con Luninosidad conocida    
    randomL =[]
    for i in range(len(logP2)): 
        randomL.append(log_Brillo(logP2[i],p3[0],p3[1])+np.random.normal(0,dev))
    randomL = array(randomL)
    
###########################################################################################
#                          Fitting Parkes et al W(logP) Model
###########################################################################################
    
    logP3 = [] # Escoger valores no nulos de P simultaneamente 
    logW3 = [] # con valores no nulos de W, para poder trabajar
    for i in range(len(logP)): 
        if (logP[i] != mask and logW[i] != mask):          
            logP3.append(logP[i])
            logW3.append(logW[i])               
    logP3 = array(logP3)
    logW3 = array(logW3)
    if len(logP3) != len(logW3): print 'AQUI HAY UN GRAVE ERROR'

    fitfunc = lambda p, x: log_W(x,p[0],p[1])  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y     # Distance to the target function
    p0 = [-3.6,0.7]  # (F,G) Initial guess for the parameters

    p4, success = optimize.leastsq(errfunc, p0[:], args=(logP3, logW3))
    
    # Fitting the random parameters of the model
    devW = np.sqrt((logW3-log_W(logP3,p4[0],p4[1])).var())  
    # El mejor estimador de sigma es sqrt(varianza)
    p4 = list(p4)
    p4.append(devW)
    p4 = array(p4)

    print 'M. Parkes W(logP) -> [F,G,sigma] ', success, p4, len(logW3), 'datos'

    # Calcula L para todos los pulsares con Luninosidad conocida    
    randomW =[]
    for i in range(len(logP3)): 
        randomW.append(log_W(logP3[i],p4[0],p4[1])+np.random.normal(0,devW))
    randomW = array(randomW)
    
###########################################################################################
#                            Parte donde grafico
#########################################################################################

    #matplotlib.pyplot.hist(x, bins=10, range=None, normed=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, hold=None, **kwargs)
    """
    k = 1./np.log(10)
    xi = k*xi    
    logP = k*logP
    logP2 = k*logP2
    logL2 = k*logL2
    logP1 = k*logP1
    randomL = k*randomL
    randomL1 = k*randomL1
    logP3 = k*logP3
    logW3 = k*logW3
    randomW = k*randomW
    randomW1 = k*randomW1
    """
    """
    plt.plot(xi,a[0],color='b',linewidth=1,label='Observado',ls='steps')
    plt.plot(xi,dist_lorimer(xi,p1[0],p1[1],p1[2]),color='g',lw=2,label='Modelado',ls='steps')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel(r'$\ln P$ [s]')
    plt.ylabel(r'$N/N_{Tot}$')
    #plt.show()
    plt.savefig('N_logP.pdf')
    plt.clf()
    """
    dumP = sp.linspace(min(logP), max(logP), 10)

    cP = np.random.normal(p1[1],p1[2],len(oP))
    cQ = np.random.normal(p1_1[1],p1_1[2],len(oP))

    deathX = [-4,3]
    deathline = [0]*2
    for i in range(len(deathX)): deathline[i] = 2.0*deathX[i]-38.0      
 
    ax1 = plt.subplot(1,2,1)
    plt.plot(oP,oQ,',',color='k',label='Observado',alpha=.950)
    plt.plot(deathX, deathline, ls='--', color='k', alpha=.50)
    #plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel(r'$\ln P$ [s]')
    plt.ylabel(r'$\ln \dot{P} $ [adim]')  
    ax2 = plt.subplot(1,2,2, sharex=ax1, sharey=ax1)  
    plt.plot(cP,cQ,',',color='k',label='Modelado',alpha=.950)
    plt.plot(deathX, deathline, ls='--', color='k', alpha=.50)
    plt.legend(loc='upper left')
    #plt.grid(True)
    plt.xlabel(r'$\ln P$ [s]')
    plt.setp(ax2.get_yticklabels(), visible=False)
    #plt.show()
    plt.subplots_adjust(wspace=.05)
    plt.savefig('logP_logQ.pdf')
    plt.clf()
    
    ax1 = plt.subplot(1,2,1)
    plt.plot(dumP,log_Brillo(dumP,p3[0],p3[1]),color='black',linewidth=2,ls='--', alpha=.5)
    plt.plot(logP2,logL2,'.',color='b',label='Observado',alpha=.50)
    plt.legend(loc='upper left')
    #plt.grid(True)
    plt.xlabel(r'$\ln P$ [s]')
    plt.ylabel(r'$\ln L$ [Jy kpc$^2$]')  
    ax2 = plt.subplot(1,2,2, sharex=ax1, sharey=ax1)  
    plt.plot(dumP,log_Brillo(dumP,p3[0],p3[1]),color='black',linewidth=2,ls='--', alpha=.5)
    plt.plot(logP2,randomL,'.',color='g',label='Modelado',alpha=.50)
    #plt.grid(True)
    plt.xlabel(r'$\ln P$ [s]')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.legend(loc='upper left')
    #plt.show()
    plt.subplots_adjust(wspace=.05)
    plt.savefig('logP_logL.pdf')
    plt.clf()

    plt.hist(logL2,10,color='b',alpha=1.,label='Observado',normed=True,histtype='step')
    plt.hist(randomL,10,color='g',alpha=1.,label='Modelado',lw=2,normed=True,histtype='step')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.xlabel(r'$\ln L$ [Jy kpc$^2$]')
    plt.ylabel(r'$N/N_{Tot}$')
    #plt.show()
    plt.savefig('N_logL.pdf')
    plt.clf()
    
    ax1 = plt.subplot(1,2,1)
    plt.plot(dumP,log_W(dumP,p4[0],p4[1]),color='black',linewidth=2,ls='--', alpha=.5)
    plt.plot(logP3,logW3,'.',color='b',label='Observado',alpha=.50)
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.xlabel('$\ln P$ [Jy kpc$^2$]')
    plt.ylabel('$\ln W$ [s]')
    ax2 = plt.subplot(1,2,2, sharex=ax1, sharey=ax1)
    plt.plot(dumP,log_W(dumP,p4[0],p4[1]),color='k',linewidth=2,ls='--', alpha=.5)
    plt.plot(logP3,randomW,'.',color='g',label='Modelado',alpha=.50)
    plt.grid(True)
    plt.xlabel('$\ln(P)$')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.legend(loc='upper left')
    #plt.show()
    plt.subplots_adjust(wspace=.05)
    plt.savefig('logP_logW.pdf')
    plt.clf()

    plt.hist(logW3,10,color='b',alpha=1.,label='Observado',normed=True,histtype='step')
    plt.hist(randomW,10,color='g',alpha=1.,label='Modelado',normed=True,lw=2,histtype='step')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel('$\ln W$ [s]')
    plt.ylabel('$N/N_{Tot}$')
    #plt.show()
    plt.savefig('N_logW.pdf')
    plt.clf()
    #"""

    plt.plot(oP,oQ,'.',color='b',label='Observado',alpha=.250)
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel(r'$\ln P$ [s]')
    plt.ylabel(r'$\ln \dot{P} $ [adim]')  
    plt.plot(cP,cQ,'.',color='g',label='Modelado',alpha=.250)
    plt.savefig('logP_logQ_mismo.pdf')
    plt.clf()

###########################################################################################
#                             DEBUGGING (y salida de datos)
###########################################################################################

    #print >>s,locals()
    #print >>s,globals()
"""
"""
s = file("atnfParameters.dat", "w") # archivo de los datos

print >>s, 'VALORES DE LAS CONSTANTES DEL MODELO' 
print >>s, 'Usados ', len(logP), 'pulsares en total' 
print >>s, 'con periodos entre ', np.exp(min(logP)), ' y ', np.exp(max(logP)), '[s]'
print >>s, 'M. Lorimer  [A,B,C]     ', success, p1, len(logP1), 'datos'
print >>s, 'Gaussiana   [A*,B,C]    ', success, p2, len(logP1), 'datos'
print >>s, 'G. en p1    [H*,I,J]    ', success, p1_1, len(logQ1), 'datos'
print >>s, 'M. L(logP)  [D,E,sigma] ', success, p3, len(logL2), 'datos'
print >>s, 'M. W(logP)  [F,G,sigma] ', success, p4, len(logW3), 'datos'
print >>s, 'Datos de: http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?P0=P0&W50=W50&R_lum=R_lum&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=p0%3E0.05&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Short+without+errors&no_value=nan&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=40&table_bottom.y=39' 

s.close()

print 'YA ESTUVO BUENO... LLEGUE AL FINAL'
