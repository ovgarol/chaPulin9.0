#!/usr/bin/env python
""" Tercera implementacion del codigo
    Segunda modulariozacion
    Datos obtenidos de atnfParameters
    Segunda construccion completa
    Debug astronomico completo: L max a 27 Jy kpc2
    Indice espectral a -1.82
    maximo brillo de burst 10e5 L mean
    Incluye envejecimiento 
    Incluye deathline
    Incluye la SFH y metalicidad como parte del modelo
    Reescala la metalicidad de [Fe/H] a [Z]
"""
import classesPulsar as K
import numpy as np
#import scipy as sp
import time as tm
import matplotlib.pyplot as plt
from scipy import optimize, array, stats
from matplotlib import rc, rcParams

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})

###########################################################################################
#                                       MAIN
###########################################################################################

def main():
    print '____Simulacion de deteccion de pulsares (8.2)____'
    identificador = tm.strftime("%H%M%S") 
    np.random.seed()
    graf = 0    # 1 es para graficar
    sfh = 0     # 1 es para el modelo completo

    #######################################################################################
    #             PARAMETROS DEL MODELO, AQUI VA LO DE atnfParameters.py
    #######################################################################################
    """
    M. Lorimer  [A,B,C]      2 [ 95.95958925  -0.5033699    0.83728586] 1886 datos
    Gaussiana   [A*,B,C]     2 [ 201.39647687   -0.50337007    0.83728508] 1886 datos
    G. en p1    [H*,I,J]     2 [ 737.5521382   -33.65963024    2.10570097] 1727 datos
    M. L(logP)  [D,E,sigma]  2 [-0.41439392 -2.50427137  1.86696577] 612 datos
    M. W(logP)  [F,G,sigma]  2 [-3.68794224  0.59786374  0.77952714] 1577 datos

    M. Lorimer  [A,B,C]      2 [ 85.26059596  -0.50643039   0.78643468] 1642 datos
    Gaussiana   [A*,B,C]     2 [ 168.07411787   -0.50643038    0.78643427] 1642 datos
    G. en p1    [H*,I,J]     2 [ 590.57158678  -33.75617791    2.03112054] 1642 datos
    M. L(logP)  [D,E,sigma]  2 [-0.22611374 -4.68845047  1.35203616] 563 datos
    M. W(logP)  [F,G,sigma]  2 [-3.01493877  0.61616299  0.6309611 ] 853 datos
    """

    p1 = [ 168.07411787,   -0.50643038,    0.78643427]    #Dist normal N_logP
    p2 = [-0.22611374, -4.68845047,  1.35203616]  #Dist L-logP
    p3 = [-3.01493877,  0.61616299,  0.6309611 ]   #Dist W-logP 
    p4 = [ 590.57158678,  -33.75617791,    2.03112054]   #Dist normal N-logQ
   
        ##########################
        # MODELO SENCILLO
        ##########################
    if sfh == 0:    # toy model
        distance = 90     # [kpc] distancia a la dSphG (90 para Sculptor, 778 para M 31)  
        totalPop = int(np.random.uniform(4752.,447.))   # numero de pulsares esperados PARAMETRO POBLACIONAL
        age = 22.01e6       # edad de la poblacion simulada en anos
        galaxy = K.genericModel('galaxia SdSph (segun atnf)',totalPop,distance,p1,p2,p3,p4,age)
    
        ##########################
        # Parametros para jugar
        ##########################

    """
      a       b       c       d       e       f       g       h       i       j       k       l       m
    $2.3$   & $2.5$ & $2.7$ & $3.0$ & $3.5$ & $2.5$ & $2.5$ & $2.5$ & $2.5$ & $2.5$ & $2.5$ & $2.5$ & $2.5$ \\ 
    $0$     & $0$   & $0$   & $0$   & $0$   & $-1$  & $-0.5$& $0.5$ & $1.0$ & $0$   & $0$   & $0$   & $0$   \\
    $-12$   & $-12$ & $-12$ & $-12$ & $-12$ & $-12$ & $-12$ & $-12$ & $-12$ & $-10$ & $-8$  & $-6$  & $-4$  \\ 

    ""
    A = 1.0e-4  # Escala en la SFR (original 1.0e-4) (adim)
    index = 2.5 #np.random.uniform(2.3,3.5)#np.random.uniform(2.3,3.5)    # indice espectral \alpha_3
    M = .0 # np.random.uniform(-1.,1.)         # Corrimiento en la metalicidad
    Y = -12.5e9#np.random.uniform(-13.0e9,-4e9)     # Corrimiento de la SFH (edad a la que observo) (yrs)

        ##########################
        # MODELO CON sfh
        ##########################

    #distance = 90     # [kpc] distancia a la dSphG (90 para Sculptor, 778 para M 31)  
    #SFH = [[13.5e9+Y,12.5e9+Y,11.5e9+Y,10.5e9+Y,9.5e9+Y,8.5e9+Y,7.5e9+Y,6.5e9+Y,5.5e9+Y],\
    #        [965*A,702*A,351*A,279*A,170*A,102*A,68*A,50*A,33*A]] # [Edad], [Npulsares]
    #SFH = [[80e6,90e6,100e6],[100,100,100]]    
    #totalPop = array(SFH[1]).sum()
    #age = array(SFH[0])
    #galaxy = K.SFHModel('galaxia SdSph (segun deBoer)',distance,p1,p2,p3,p4,SFH)

        ##########################
        # MODELO CON completo
        ##########################
    if sfh == 1: # los datod de sculptor
        distance = 90     # [kpc] distancia a la dSphG (90 para Sculptor, 778 para M 31)  
        datos_g = [[13.5e9+Y,12.5e9+Y,11.5e9+Y,10.5e9+Y,9.5e9+Y,8.5e9+Y,7.5e9+Y,6.5e9+Y,5.5e9+Y],\
            [27.5*A,19.5*A,11.*A,8.*A,5.5*A,2.5*A,2.0*A,1.75*A,1.0*A],\
            [-2.40+M,-2.02+M,-1.71+M,-1.61+M,-1.49+M,-1.26+M,-1.19+M,-1.19+M,-1.12+M],\
            [1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9]] # [Edad], [sfr], [Z], [ancho de burst] 
        age = array(datos_g[0])
        galaxy = K.CompletModel('galaxia SdSph (segun deBoer)',distance,p1,p2,p3,p4,datos_g,index)
        totalPop = galaxy.TotalPulsares()

    if sfh == 2: # un solo burst para pruebas de evolucion
        distance = 90     # [kpc] distancia a la dSphG (90 para Sculptor, 778 para M 31)  
        datos_g = [[1e7+Y],\
            [30*A],\
            [0.0+M],\
            [1e6]] # [Edad], [sfr], [Z], [ancho de burst] 
        age = array(datos_g[0])
        galaxy = K.CompletModel('galaxia SdSph (segun deBoer)',distance,p1,p2,p3,p4,datos_g,index)
        totalPop = galaxy.TotalPulsares()
    """
    #######################################################################################
    #                      AQUI CORRE LAS RUTINAS DE clasesPulsar.py
    #######################################################################################

    startTime = tm.time()

    """
    Radiotelescopio G [K/Jy] Tsys [K] np nu [MHz] delta nu [MHz] beta
    GBT     2.0     350     2   350     50   1.5
    WT      1.0     140     2   328     10   1.3
    ATA 42  0.28    44      2   1420    100  1.16
    ATCA    0.86    30      26  1420    300  1.16
    LOFAR   5.8     1000    2   140     25.6 1.0
    SKA1    7.91    28      2   400     700  1.16

    """


    #Creacion de los Radio telescopios    

    name, G, beta, npol, deltaF, Tsys, f_obs = 0, 2, 1.5, 2, 50.*1e6, 70, 350.  #parametros de GBT
    #name, G, beta, npol, deltaF, Tsys, f_obs = 1, 1, 1.3, 2, 1.*1e7, 140, 328.  #parametros de WT
    #name, G, beta, npol, deltaF, Tsys, f_obs = 2, 0.28, 1.16, 2, 100.*1e6, 44, 1420.  #parametros de ATA 42
    #name, G, beta, npol, deltaF, Tsys, f_obs = 3, 0.86, 1.16, 2, 300.*1e6, 30, 1420.  #parametros de ATCA
    #name, G, beta, npol, deltaF, Tsys, f_obs = 4, 5.8, 1.0, 2, 25.6*1e6, 1000, 140  #parametros de LOFAR
    #name, G, beta, npol, deltaF, Tsys, f_obs = 5, 7.91, 1.16, 2, 700*1e6, 28, 400  #parametros de SKA1

    SN = 5.     # Signal/Noise para la detectabilidad
    totalTime = 6./60.       # tiempo de integracion [hrs] -> 720 s -> 12 min 

    RT = K.genericTelescope(name,G,beta,npol,deltaF,Tsys,f_obs) 
           
    p, q, l, w, b, D, sigma, period, pulsos, totals = RT.pulsarDetect(SN, galaxy, totalTime)
               
    endTime = tm.time()
    tTime = (endTime - startTime)/(totalPop*totalTime)
    
    #AQUI LOS RESULTADOS LOS MUESTRA EN PANTALLA
    """
    print "Tiempo de ejecucion: ",(endTime - startTime)," s - ", tTime, "[s/pulsar*hr]"    
    print "Se generaran ",totalPop ," pulsares." 
    print "Se detectan ", array(period).sum(), " por periodicidad."
    print "Se detectan ", array(pulsos).sum(), " por pulsos."
    print "Se detectan ", array(totals).sum(), " en total"
    """
    #######################################################################################
    #                      AQUI VOY A GUARDAR LOS DATOS
    #######################################################################################
    #"""
    nombre = 'out/data/sim_' + identificador + '_' + str(startTime) +'.dat'
    s = file(nombre, "w") # archivo de los datos

    print >>s,'# ESTOS SON LOS RESULTADOS DE LA SIMULACION'    
    print >>s,"# Tiempo de ejecucion: ",(endTime - startTime),"s-",tTime, "[s/pulsar*hr]"  
    print >>s,'# Galaxia ', galaxy.name, " con ", galaxy.nTot, " pulsares con un burst hace", age, "yrs" 
    print >>s,'# A una distancia de ', galaxy.distance, "[kpc]" 
    print >>s,"# Durante una observacion de ", totalTime, " [hrs]" 
    print >>s,"Se detectan ", array(period).sum(), " por periodicidad."
    print >>s,"Se detectan ", array(pulsos).sum(), " por pulsos."
    print >>s,"Se detectan ", array(totals).sum(), " en total"
    print >>s,'#i,  P0,   P1,   L,   W,  B,  Tipo,   Sigma,   Periodic, Pulse, Detectado'
    print >>s,'##, [s],  [adim],  [Jy kpc^2], [G],   [s], [012D],  [adim],   [01],  [01],   [01]'  

    for i in range(totalPop): 
        if D != 'D': print >>s, i, p[i], q[i], l[i], w[i], b[i], D[i], sigma[i], period[i], pulsos[i], totals[i]
    
    print >>s, '# #### ### ## # ------NADA MAS------ # ## ### #### #'

    s.close()
    #"""
    #######################################################################################
    #                      AQUI VOY A GRAFICAR LOS DATOS
    #######################################################################################

    pPeriodic=[]
    pPulse=[]
    pN=[]
    qPulse=[]
    qPeriodic=[]
    qN=[]
    lPeriodic=[]
    lPulse=[]
    lN=[]    
    wPeriodic=[]
    wPulse=[]
    wN=[]
    bPeriodic=[]
    bPulse=[]
    bN=[]
    DPeriodic=[]
    DPulse=[]
    DN=[]
    sPeriodic=[]
    sPulse=[]
    sN=[]

    if D != 'D':
        for i in range(totalPop):
            if not np.isinf(q[i]): 
                pN.append(p[i])
                qN.append(q[i])            
                lN.append(l[i])
                wN.append(w[i])
                bN.append(b[i])
                DN.append(D[i])
                sN.append(sigma[i])
            if period[i] == 1: 
                pPeriodic.append(p[i])
                qPeriodic.append(q[i])            
                lPeriodic.append(l[i])
                wPeriodic.append(w[i])
                bPeriodic.append(b[i])
                DPeriodic.append(D[i])
                sPeriodic.append(sigma[i])
            if pulsos[i] == 1: 
                pPulse.append(p[i]) 
                qPulse.append(q[i])        
                lPulse.append(l[i]) 
                wPulse.append(w[i]) 
                bPulse.append(b[i])
                DPulse.append(D[i])
                sPulse.append(sigma[i])

    tipo0, tipo1, tipo2, tipo2_a = 0., 0., 0., 0.

    for i in range(len(DPeriodic)):
        if DPeriodic[i] == 0: tipo0 = tipo0 + 1
        elif DPeriodic[i] == 1: tipo1 = tipo1 + 1
        elif DPeriodic[i] == 2 and sPeriodic[i]<=1.: tipo2 = tipo2 + 1
        elif DPeriodic[i] == 2 and sPeriodic[i]>1.: tipo2_a = tipo2_a + 1

    zipo0, zipo1, zipo2, zipo2_a = 0., 0., 0., 0.

    for i in range(len(DPulse)):
        if DPulse[i] == 0: zipo0 = zipo0 + 1
        elif DPulse[i] == 1: zipo1 = zipo1 + 1
        elif DPulse[i] == 2 and sPulse[i]<=1.: zipo2 = zipo2 + 1
        elif DPulse[i] == 2 and sPulse[i]>1.: zipo2_a = zipo2_a + 1

    N_tipo_periodicos   = [tipo0, tipo1, tipo2, tipo2_a]
    N_tipo_pulsos   = [zipo0, zipo1, zipo2, zipo2_a]

    if int(identificador)%100 == 0: graf = 0
    if graf == 1:
        for i in range(5): 
            if galaxy.pulsar[i].dist != 'D': K.pulseHist(galaxy.pulsar[i],totalTime,identificador)

        deathX = [-3,2]
        deathline = [0]*2
        for i in range(len(deathX)): deathline[i] = 2.0*deathX[i]-16.52+1    
        plt.plot(deathX, deathline, ls='--', color='k', alpha=.50)
        plt.plot(pN,qN,',',color='k',alpha=.850)
        plt.plot(pPeriodic,qPeriodic,'^',color='orange',label='Periodicidad',lw='1',alpha=1.)
        plt.plot(pPulse,qPulse,'v',color='red',label='Pulsos',lw='1',alpha=1.)
        plt.legend(loc='best')
        plt.xlabel(r'$\log P$ [s]')
        plt.ylabel(r'$\log \dot{P}$ [s/s]')
        plt.xlim([-2,2])
        plt.ylim([-18,-12])
        plt.grid(True)
        #plt.show()
        plt.savefig('out/img/logP_logPdot_' + identificador + '.pdf')
        plt.clf()
    
        #plt.plot(deathX, deathline, ls='-', color='k')
        plt.plot(pN,bN,',',color='k',alpha=.850)
        plt.plot(pPeriodic,bPeriodic,'^',color='orange',label='Periodicidad',lw='1',alpha=1.)
        plt.plot(pPulse,bPulse,'v',color='red',label='Pulsos',lw='1',alpha=1.)
        plt.legend(loc='best')
        plt.xlabel(r'$\log P$ [s]')
        plt.ylabel(r'$\log B$ [G]')
        #plt.xlim([-2,2])
        #plt.ylim([-18,-12])
        plt.grid(True)
        #plt.show()
        plt.savefig('out/img/logP_logB_' + identificador + '.pdf')
        plt.clf()
        
        #plt.plot(deathX, deathline, ls='-', color='k')
        plt.plot(pN,lN,',',color='k',alpha=.850)
        plt.plot(pPeriodic,lPeriodic,'^',color='orange',label='Periodicidad',lw='1',alpha=1.)
        plt.plot(pPulse,lPulse,'v',color='red',label='Pulsos',lw='1',alpha=1.)
        plt.legend(loc='best')
        plt.xlabel(r'$\log P$ [s]')
        plt.ylabel(r'$\log L$ [Jy kpc$^2$]')
        #plt.xlim([-2,2])
        #plt.ylim([-18,-12])
        plt.grid(True)
        #plt.show()
        plt.savefig('out/img/logP_logL_detected_' + identificador + '.pdf')
        plt.clf()
        
        #plt.plot(deathX, deathline, ls='-', color='k')
        plt.plot(pN,wN,',',color='k',alpha=.850)
        plt.plot(pPeriodic,wPeriodic,'^',color='orange',label='Periodicidad',lw='1',alpha=1.)
        plt.plot(pPulse,wPulse,'v',color='red',label='Pulsos',lw='1',alpha=1.)
        plt.legend(loc='best')
        plt.xlabel(r'$\log P$ [s]')
        plt.ylabel(r'$\log W$ [s]')
        #plt.xlim([-2,2])
        #plt.ylim([-18,-12])
        plt.grid(True)
        #plt.show()
        plt.savefig('out/img/logP_logW_detected_' + identificador + '.pdf')
        plt.clf()
    
        if len(pN) !=0:
            #rango = (array(pN).min(),array(pN).min())
            plt.hist(pN,10,color='k',linewidth=1,histtype='step', label='Total')#, range=rango)
            if len(pPeriodic) !=0: plt.hist(pPeriodic,10,color='orange',linewidth=2,histtype='step', label='Periodico')#, range = rango)
            if len(pPulse) !=0: plt.hist(pPulse,10,color='red',linewidth=2,histtype='step', label='Pulsos')#,range = rango)
            plt.legend(loc='upper left')
            plt.xlabel('$\log P [s]$')
            plt.ylabel('$N$')
            plt.grid(True)
            #plt.show()
            plt.savefig('out/img/N_logP_' + identificador + '.pdf')
            plt.clf()

            #rango = (array(qN).min(),array(qN).min())
            plt.hist(qN,10,color='k',linewidth=1,histtype='step', label='Total')#,range=rango)
            if len(qPeriodic) !=0: plt.hist(qPeriodic,10,color='orange',linewidth=2,histtype='step', label='Periodico')#, range = rango)
            if len(qPeriodic) !=0: plt.hist(qPulse,10,color='red',linewidth=2,histtype='step', label='Pulsos')#, range = rango)
            plt.legend(loc='upper left')
            plt.xlabel('$\log \dot{P} [s/s]$')
            plt.ylabel('$N$')
            plt.grid(True)
            #plt.show()
            plt.savefig('out/img/N_logPdot_' + identificador + '.pdf')
            plt.clf()
    
    #######################################################################################
    #                      AQUI VOY A GUARDAR LOS DATOS
    #######################################################################################

    nombre = 'out/data/sim_' + identificador + '_' + str(startTime)+'.det'
    s = file(nombre, "w") # archivo de los datos

    print >>s,'# ESTOS SON LOS RESULTADOS DE LA SIMULACION'    
    print >>s,'# Galaxia ', galaxy.name, " con ", galaxy.nTot, " pulsares con un burst hace", age, "yrs" 
    print >>s,'# A una distancia de ', galaxy.distance, "[kpc]" 
   
    print >>s,'#i,  P0,   P1,   L,   W,  B,  Tipo,   Sigma,   Periodic, Pulse'
    print >>s,'##, [s],  [adim],  [Jy kpc^2], [G],   [s], [012D],  [adim],   [01],  [01]'  
   
    n = 0
    for i in range(totalPop):
        if totals[i] == 1 or not np.isinf(p[i]): 
            print >>s, i, p[i], q[i], l[i], w[i], b[i], D[i], sigma[i], period[i], pulsos[i]
            n = n+1        
        
    print >>s, '# #### ### ## # ----------- # ## ### #### #'
    print >>s, 'En total son', n, ' pulsares generados.'
      
    pN = array(pN)    
    qN = array(qN)    
    lN = array(lN)    
    wN = array(wN)
        
    print >>s, 'P', pN.mean(), '-', pN.std() 
    print >>s, 'Q', qN.mean(), '-', qN.std() 
    print >>s, 'L', lN.mean(), '-', lN.std() 
    print >>s, 'W', wN.mean(), '-', wN.std() 
    
    s.close()
    if sfh == 1:
        nombregal = 'sim_parametros.gal'
        gal = file(nombregal, "a") # archivo de los datos
        print >>gal, galaxy.showParametros()[0], galaxy.showParametros()[1], galaxy.showParametros()[2], galaxy.showParametros()[3], n, M, -Y
        gal.close()
       
    #######################################################################################
    # AQUI VOY A GUARDAR LOS DATOS DE dumC, periodicos, dumD, pulsos Y dumT EN OTRO ARCHIVO
    #######################################################################################
    
    resumen = 'out/out_simPulsar.dat'
    ss = file(resumen, "a") # archivo de resumen
    print >>ss, nombre, "\t", array(period).sum(), "\t", array(pulsos).sum(), "\t", array(totals).sum()
    ss.close()    

    #######################################################################################
    #                             return 'ok' si termina bien
    #######################################################################################

    return 'ok', array(N_tipo_periodicos), array(N_tipo_pulsos)

###########################################################################################
#                          AQUI SE CORRE MAIN()
###########################################################################################

if __name__ == "__main__":
    main()


