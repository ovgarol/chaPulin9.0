#!/usr/bin/env python
""" Calcula y grafica los parametros para un modelo basado
    en las estadisticas de ATNF.
"""

import numpy as np
#import scipy as sp
#import time as tm
import matplotlib.pyplot as plt
from scipy import optimize, array, stats
from matplotlib import rc, rcParams

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})

###########################################################################################
#                                       MAIN
###########################################################################################
if __name__ == "__main__":

###########################################################################################
#                              Getting data from file
###########################################################################################

    metal = [[-1.,-0.5,0.0,0.5,1.0],\
    [7128,7179,7194,6778,6660],\
    [689,588,639,726,611],\
    [33.18,33.86,33.82,31.32,30.36],\
    [6.33,6.41,6.79,6.38,6.26]]

    f = open('sim_parametros.gal', 'r') # open file to read
    index = []
    Np = []
    #C0 = []
    #C1 = []
    Nc = []
    M = []
    D = []

    # index C0 C1 Np  

    for line in f:      # iterate over the lines in the file

        columns = line.split(' ')       # split the line into a list of column values
        columns = [col.strip() for col in columns]  # clean any whitespace off the items

        index.append(float(columns[0])) 
        #C0.append(np.log(float(columns[1])))      
        #C1.append(np.log(float(columns[2])))     
        Np.append((float(columns[3])))
        Nc.append((float(columns[4])))
        M.append((float(columns[5])))
        D.append((float(columns[6])))

    f.close()   

    index = array(index)
    Np = array(Np)
    #C0 = array(C0)
    #C1 = array(C0)
    Nc = array(Nc)
    M = array(M)
    D = array(D)
    
    logNp = np.log10(Np)
    logNc = np.log10(Nc)


    """
    plt.hist(C0,10,color='k',linewidth=2,histtype='step', label='C0')#, range=rango)
    plt.legend(loc='upper left')
    plt.xlabel('$\log C0$')
    plt.ylabel('$N$')
    plt.grid(True)
    #plt.show()
    plt.savefig('N_logC0_.pdf')
    plt.clf()

    plt.hist(C1,10,color='k',linewidth=2,histtype='step', label='C1')#, range=rango)
    plt.legend(loc='upper left')
    plt.xlabel('$\log C1$')
    plt.ylabel('$N$')
    plt.grid(True)
    #plt.show()
    plt.savefig('N_logC1_.pdf')
    plt.clf()
    """

    plt.hist(logNp,10,color='k',linewidth=2,histtype='step')#, range=rango)
    #plt.legend(loc='upper left')
    plt.xlabel('$\log N_T$ [adim]')
    plt.ylabel('$N$')
    plt.grid(True)
    #plt.show()
    plt.savefig('N_logNT_.pdf')
    plt.clf()
   
    plt.plot(index,Np,',',color='k',alpha=1.)
    #plt.legend(loc='best')
    plt.xlabel(r'$\alpha_3$ [adim]')
    plt.ylabel('$N_T$ [adim]')
    #plt.xlim([-2,2])
    #plt.ylim([-18,-12])
    plt.grid(True)
    #plt.show()
    plt.savefig('index_NT.pdf')
    plt.clf()
    

    plt.hist(logNc,10,color='k',linewidth=2,histtype='step')#, range=rango)
    #plt.legend(loc='upper left')
    plt.xlabel('$\log N_p$ [adim]')
    plt.ylabel('$N$')
    plt.grid(True)
    #plt.show()
    plt.savefig('N_logNp_.pdf')
    plt.clf()
    ""
    plt.plot(index,Nc,',',color='k',alpha=1.)
    #plt.legend(loc='best')
    plt.xlabel(r'$\alpha_3$ [adim]')
    plt.ylabel('$N_p$ [adim]')
    #plt.xlim([-2,2])
    #plt.ylim([-18,-12])
    plt.grid(True)
    #plt.show()
    plt.savefig('index_Np.pdf')
    plt.clf()

    plt.plot(M,Np,',',color='k',alpha=1.)
    plt.errorbar(metal[0],metal[1], yerr = metal[2], fmt='o', ls='--',color='k')
    #plt.legend(loc='best')
    plt.xlabel(r'$\Delta$ [Z] [dex]')
    plt.ylabel('$N_T$ [adim]')
    plt.xlim([-1.5,1.5])
    #plt.ylim([-18,-12])
    plt.grid(True)
    #plt.show()
    plt.savefig('M_NT.pdf')
    plt.clf()

    plt.plot(M,Nc,',',color='k',alpha=1.)
    plt.errorbar(metal[0],metal[3], yerr=metal[4],fmt='o',color='k',ls='--')
    #plt.legend(loc='best')
    plt.xlabel(r'$\Delta$ [Z] [dex]')
    plt.ylabel('$N_p$ [adim]')
    plt.xlim([-1.5,1.5])
    #plt.ylim([-18,-12])
    plt.grid(True)
    #plt.show()
    plt.savefig('M_Np.pdf')
    plt.clf()


    plt.plot(D,Np,',',color='k',alpha=1.)
    #plt.legend(loc='best')
    plt.xlabel(r'$\Delta t$ [yr]')
    plt.ylabel('$N_T$ [adim]')
    #plt.xlim([-2,2])
    #plt.ylim([-18,-12])
    plt.grid(True)
    #plt.show()
    plt.savefig('D_NT.pdf')
    plt.clf()

    plt.plot(D,Nc,',',color='k',alpha=1.)
    #plt.legend(loc='best')
    plt.xlabel(r'$\Delta t$ [yr]')
    plt.ylabel('$N_p$ [adim]')
    #plt.xlim([-2,2])
    #plt.ylim([-18,-12])
    plt.grid(True)
    #plt.show()
    plt.savefig('D_Np.pdf')
    plt.clf()

    
    DD = open('sim.gal', 'a') # open file to read
    print>>DD, 'N_T - ', Np.mean(), '-', Np.std()
    print>>DD, 'N_p - ', Nc.mean(), '-', Nc.std()
    DD.close()
    


#plt.plot(age[0],age[1],'o',color='k',alpha=.25)

plt.legend(loc='best')
plt.xlabel(r'$\Delta  Z$ [dex]')
plt.ylabel(r'$N_{T}$ [adim]')
plt.xlim([-1.5,1.5])
plt.grid(True)
#plt.show()
plt.savefig('metal_Nt.jpg')
plt.clf()


#plt.errorbar(index[0],index[1], yerr = index[2], fmt='o', ls='--')
plt.legend(loc='best')
plt.xlabel(r'$\Delta Z$ [dex]')
plt.ylabel(r'$N_{p}$ [adim]')
plt.xlim([-1.5,1.5])
plt.grid(True)
#plt.show()
plt.savefig('metal_Np.jpg')
plt.clf()

dumi = []
dumd = []
for i in range(len(metal[0])):
    dumi.append(metal[3][i]/metal[1][i])
    dumd.append(metal[4][i]/metal[1][i])

#plt.plot(age[0],age[1],'o',color='k',alpha=.25)
plt.errorbar(metal[0],dumi, yerr = dumd, fmt='o', ls='--',color='k')
plt.legend(loc='best')
plt.xlabel(r'$\Delta  Z$ [dex]')
plt.ylabel(r'$N_{p}/N_{T}$ [adim]')
plt.xlim([-1.5,1.5])
plt.grid(True)
#plt.show()
plt.savefig('metal_Np-Nt.jpg')
plt.clf()

del dumi, dumd

