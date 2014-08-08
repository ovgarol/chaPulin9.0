#!/usr/bin/env python
""" Tercera implementacion del codigo
    Segunda modulariozacion
    Datos obtenidos de atnfParameters
    Segunda construccion completa
    Debug astronomico completo: L max a 27 Jy kpc2
    Indice espectral a -1.82
    maximo brillo de burst 10e5 L mean
    Ahora modela p1 para sacar el nivel de nulling en B, 
    basado en Bsurf=3.2e19*(p0*p1)^(1/2)<2 Gauss
    solo para envejecer la poblacion y decir cuales no se verian
    chaPulin_8.2
    Primera revision con modelos completos
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
#         CLASE genericPulsar ALMACENA LA INFORMACION DE LOS PULARES DE UNA GALAXIA
###########################################################################################

class genericPulsar(object):
    def __init__(self, name, logP, logQ, age, p2, p3): # BASE NATURAL
        """Input: Caracteristicas galaxyales P, Q, L, W, age (name es el no. de iteracion!!) 
        """
        self.spectral_index = np.random.normal(-1.8,0.2)  # correccion por banda de observacion  SpectralIndex_Maron
        self.name = int(name)   # entero solo para id
        self.sigma = 0          # este numero va a ser sigmaN, sigmaLN o alpha segun toque       
        self.pulses = []        # vector vacio donde van los pulsos 
        self.maxburst = np.random.normal(100,10000)
       
        ##########################################################
        #  Un primer modelo de envejecimiento, TOMADO DE GUSINOV #
        T0 = float(np.exp(logP))        # periodo inicial al tiempo tc (s)                              
        T1 = float(np.exp(logQ))        # primera derivada inicial (s/s)
        tc = T0/(2*T1)/(3.6*2.4*3.6e6)  # edad caracteristica (yr) 

        age = age - 1.e6*np.random.uniform(7.,22.)  # Consideracion del tiempo de secuencia principal

        tm = np.random.uniform(5e6,5e7) # tiempo de evolucion (yr)
  
        if age > 0.0 and age > tc: 
            B = np.sqrt(T0*T1)*np.exp(-2*(age-tc)/tm)   # decaimiento segun gusinov
            T0 = np.sqrt(T0**2+0.5*(3.6*2.4*3.6e6)*tm*T0*T1*(1-np.exp(-2*(age-tc)/tm)))
            T1 = B**2/T0 
        elif age < 0.0:
            T0 = 0.0
            T1 = 0.0 
        ##########################################################    

        if T1 != 0.0:
            logQcrit = 2.0*np.log10(T0)-16.52 + 1 + np.random.normal(0,0.5/2)    # deathline segun BING ZHAN
            if np.log10(T1) > logQcrit and not np.isinf(T1):        
                self.P = T0        # s              
                self.Q = T1        # adim    
                
                # brillo promedio mas grande permitido valor del atnf 26100 mJy
                self.L = float(np.exp(log_Brillo(self.P,p2[0],p2[1],p2[2])))      
                while self.L > 27: self.L = float(np.exp(log_Brillo(self.P,p2[0],p2[1],p2[2])))   
    
                self.W = float(np.exp(log_W(self.P,p3[0],p3[1],p3[2])))        # s 
                         
                self.dist = np.random.randint(0,3)  # entero 0,1,2--aleatorio               
                if self.dist == 0: self.sigma = np.random.uniform(1,100)     #normal
                elif self.dist == 1: self.sigma = np.random.uniform(1.01,100) #lognormal
                elif self.dist == 2: self.sigma = np.random.uniform(0,3)     #powerlaw

            else:   # esta debajo de la dead line 
                self.P = 0.0         # s              
                self.Q = 0.0         # adim   
                self.dist = 'D'      # invisible   
                self.L = 0.0         # Jk kpc^2
                self.W = 0.0         # s 

        else:   # el pulsar no ha nacido
            self.P = 0.0         # s              
            self.Q = 0.0         # adim   
            self.dist = 'D'      # invisible   
            self.L = 0.0         # Jk kpc^2
            self.W = 0.0         # s 

    def observeIt(self, t):
        maxBurst = self.maxburst * self.L
        if self.dist == 'D':# invisible
            pulses = [0]    # si es invisible no emite pulsos
        else:        
            maxT = int(t*3600/self.P)   # numero de iteraciones que depende del periodo
            pulses = [0]*maxT           # estoy generando una lista de longitud maxT llena de ceros
        
        if self.dist == 0: # distribucion normal
            i=0
            while i < maxT:
                x = np.random.normal(self.L,self.sigma) 
                if(x>0) and x/self.L <= maxBurst: 
                    pulses[i] = x 
                    i = i + 1                  

        if self.dist == 1: # distribucion lognormal
            i=0
            while i < maxT:
                z = pulses[i] = np.random.lognormal(np.log(self.L)-0.5*np.log(self.sigma)**2,np.log(self.sigma))  
                if(z/self.L <= maxBurst):  #se reducen los picos hasta 10^5 el promedio
                    pulses[i] = z 
                    i = i + 1  

        if self.dist == 2: # distribucion de powerlaw
            i=0
            while i < maxT:
                y = self.L*self.sigma/(self.sigma+1) + np.random.pareto(self.sigma) 
                if(y/self.L <= maxBurst):  #se reducen los picos hasta 10^5 el promedio
                    pulses[i] = y 
                    i = i + 1   
  
        self.pulses = pulses
        return pulses

    def pulseIt(self):
        maxBurst = self.maxburst*self.L
        if self.dist == 'D': #invisible
            return 0     
          
        if self.dist == 0: #normal
            x = np.random.normal(self.L,self.sigma)
            while x < 0 or x/self.L > maxBurst:
                x = np.random.normal(self.L,self.sigma)       
            return x

        if self.dist == 1: #lognormal
            z = np.random.lognormal(np.log(self.L)-0.5*np.log(self.sigma)**2,np.log(self.sigma))
            while z/self.L > maxBurst:  #se reducen los picos hasta 10^5 el promedio
                z = np.random.lognormal(np.log(self.L)-0.5*np.log(self.sigma)**2,np.log(self.sigma))  
            return z

        if self.dist == 2: #powerlaw
            y = self.L*self.sigma/(self.sigma+1) + np.random.pareto(self.sigma)
            while y/self.L > maxBurst:  #se reducen los picos hasta 10^5 el promedio
                y = self.L*self.sigma/(self.sigma+1) + np.random.pareto(self.sigma) 
            return y
            
    def showIt(self): # Exporta los datos del objeto pulsar a una lista, para extraer datos
        return self.name, self.P, self.Q, self.L, self.W, self.dist, self.sigma

###########################################################################################
#                                DEFINICIONES PRELIMINARES
###########################################################################################

def log_Brillo(x,D,E,sigmaL):
    return D*x + E + np.random.normal(0,sigmaL)

def log_W(x,F,G,sigmaW):
    #return 1.06*np.log(F*np.exp(x)**0.9+G)+ np.random.normal(0,sigmaW)
    return G*x + F + np.random.normal(0,sigmaW)

###########################################################################################
#           CLASE genericModel EN EL VAN A ESTAR LOS PULSARES
###########################################################################################

class genericModel(object):
    def __init__(self, name, nTot, distance, p1, p2, p3, p4, age):
        self.name = name            # solo para id
        self.nTot = int(nTot)       # Cantidad de objetos
        self.distance = float(distance)   # distancia de la galaxia
        self.p1 = p1    # modelo para distribucion de P0 
        self.p2 = p2    # modelo para distribucion de L(P)
        self.p3 = p3    # modelo para distribucion de W(P)
        self.p4 = p4    # modelo para distribucion de P1      
        self.pulsar = [0]*self.nTot    # Espacio para los pulsares del tamano requerido para no fragmentar 
      
        for i in range(self.nTot):
            x = np.random.normal(p1[1],p1[2])   # aplicacion del modelo de P0
            q = np.random.normal(p4[1],p4[2])   # aplicacion del modelo de P1      
            self.pulsar[i]=(genericPulsar(i,x,q,age,self.p2,self.p3))   # crea los pulsares

        del x, q, p1, p2, p3, p4, name, age, nTot, distance  

    def graficar(self): # devuelve los valores de cada pulsar individual
        p=[]        
        q=[]
        l=[]
        w=[]
        b=[]
        D=[]
        sigma=[]

        for i in range(self.nTot):  # BASE 10  
            D.append(self.pulsar[i].dist)   
            p.append(np.log10(self.pulsar[i].P))
            l.append(np.log10(self.pulsar[i].L))
            w.append(np.log10(self.pulsar[i].W))
            q.append(np.log10(self.pulsar[i].Q))
            sigma.append(self.pulsar[i].sigma)
            b.append(np.log10(3.2e19*(self.pulsar[i].P*self.pulsar[i].Q)**(1./2.)))
        return p, q, l, w, b, D, sigma  

###########################################################################################
#           CLASE SFHModel Se considera la sfh como un grupo de burst
###########################################################################################

class SFHModel(object):
    def __init__(self, name, distance, p1, p2, p3, p4, SFH):
        # LA SFH SE INTRODUCE COMO UNA LISTA DE (EDAD, NUMERO ESPERADO DE PULSARES)

        self.name = name            # solo para id
        self.distance = float(distance)   # distancia de la galaxia
        self.p1 = p1  
        self.p2 = p2  
        self.p3 = p3  
        self.p4 = p4         
        self.pulsar = [0]*N    # Espacio para los pulsares del tamano requerido para no fragmentar 

        age = SFH[0]        # edad del burst
        cuantos = SFH[1]    # numero de pulsares esperados
        self.nTot = array(cuantos).sum()    # numero total de pulsares
        J = 0
        N = 0
        for i in range(len(cuantos)): N = N + cuantos[i]

        for i in range(len(cuantos)):
            for j in range(int(cuantos[i])):
                x = np.random.normal(p1[1],p1[2])
                q = np.random.normal(p4[1],p4[2])         
                self.pulsar[j+J]=(genericPulsar(j+J ,x,q,age[i],self.p2,self.p3))   # crea el pulsar
            J = J + int(cuantos[i])
            j=0
        
        del q, x, p1, p2, p3, p4, name, SFH, cuantos, distance  

    def graficar(self):
        p=[]        
        q=[]
        l=[]
        w=[]
        b=[]
        D=[]
        sigma=[]
        for i in range(self.nTot):  # BASE 10  
            D.append(self.pulsar[i].dist)   
            p.append(np.log10(self.pulsar[i].P))
            l.append(np.log10(self.pulsar[i].L))
            w.append(np.log10(self.pulsar[i].W))
            q.append(np.log10(self.pulsar[i].Q))
            sigma.append(self.pulsar[i].sigma)
            b.append(np.log10(3.2e19*(self.pulsar[i].P*self.pulsar[i].Q)**(1./2.)))
        return p, q, l, w, b, D, sigma


###########################################################################################
#           CLASE CompletModel Se considera la sfh, la imf  y la metalicidad 
###########################################################################################

class CompletModel(object):
    def __init__(self, name, distance, p1, p2, p3, p4, SFH, index):
        # LA SFH SE INTRODUCE COMO UNA LISTA DE (EDAD, SFR, Z, ancho_de_burst)

        age = SFH[0]        # edad de cada burst (yr)
        sfr = SFH[1]        # star formation rate (m_sol/yr)
        metal = SFH[2]      # [Z]
        width = SFH[3]      # ancho del burst (yr)
        Mtotal = 0.0        # masa estelar total M_star (m_sol)
        cuantos = [0]*len(sfr)  # numero de bursts
        Mmin = [0]*len(sfr)     # masa minima para un SNeII (m_sol)
        Mmax = [0]*len(sfr)     # masa minima para un SNeII (m_sol)
        Mu = 100.0  # masa maxima en la galaxia (m_sun)  

        self.index = index
        self.name = name                # solo para id
        self.distance = float(distance) # distancia de la galaxia
        self.p1 = p1  # igual que el anterior
        self.p2 = p2  
        self.p3 = p3  
        self.p4 = p4         

        #### parte anulada calculo INCORRECTO
        """        
        Z0 = 3.55e-5
        Z1 = 1.404e-4
        Z2 = 1.3614e-3
        Z3 = -1.699

        for i in range(len(sfr)):
            m = np.log10(Z0*10**metal[i]/(Z1*10**metal[i]+Z2))-Z3 # reescalar Z en cuncion de [Fe/H]
            Mmin[i] = 9.40858+1.14548*m+3.96e-1*m**2+2.96e-2*m**3-8.79e-3*m**4-1.96e-3*m**5-1.12e-4*m**6
            Mtotal = Mtotal + sfr[i]*width[i] 
       """##### 

        for i in range(len(sfr)):
            if metal[i] > 0.0: Mmax[i] = Mu
            else: Mmax[i] = 25.0
            m = metal[i]
            Mmin[i] = 9.40858+1.14548*m+3.96e-1*m**2+2.96e-2*m**3-8.79e-3*m**4-1.96e-3*m**5-1.12e-4*m**6
            Mtotal = Mtotal + sfr[i]*width[i] 

        a0 = np.random.normal(0.3,0.7/2)
        a1 = np.random.normal(1.3,0.5/2)
        a2 = np.random.normal(2.3,0.3/2)

        C0 = (0.08**(2-a0)-0.01**(2-a0))/(2.-a0)+(0.5**(2-a1)-0.08**(2-a1))*0.08**(a1-a0)/(2.-a1)+(1.0**(2-a2)-0.5**(2-a2))*0.08**(a1-a0)*0.5**(a2-a1)/(2.-a2) # c_0 parametro de la imf
        C1 =  0.08**(a1-a0)*0.5**(a2-a1)*1.**(index-a2)    # c_1 parametro de la imf

        K3 =C1/(C0+C1/(2-index)*(Mu**(2-index)-1))    # parametro de la imf 
        
        self.C0 =C0
        self.C1 =C1

        for i in range(len(sfr)):
            cuantos[i] = int(1./6.*K3*sfr[i]*width[i]*(Mmax[i]**(1-index)-Mmin[i]**(1-index))/(1-index))   # aplicacion del modelo de la imf

        self.nTot = array(cuantos).sum()    # numero total de pulsares
        J = 0
        N = 0
        for i in range(len(cuantos)): N = N + cuantos[i]
        self.pulsar = [0]*N    # Espacio para los pulsares del tamano requerido para no fragmentar 

        # imprime los resultados del modelo imf-sfh-metalicidad
        print "age, sfr, metal, width, Mmin, cuantos"         
        for i in range(len(cuantos)):
            print age[i], sfr[i], metal[i], width[i], Mmin[i], cuantos[i] 
        #print C0,C1,K3
        print 'Masa total: ', Mtotal, 'En total: ', self.nTot, ' pulsares', 'Indice: ', index
        
        for i in range(len(cuantos)):
            for j in range(int(cuantos[i])):
                x = np.random.normal(p1[1],p1[2])   # aplicacion del modelo para P0
                q = np.random.normal(p4[1],p4[2])   # aplicacion del modelo para P1
                age_i = age[i]+np.random.uniform(-width[i]/2,width[i]/2)    # distribuye el nacimiento en el ancho del burst
                self.pulsar[j+J]=(genericPulsar(j+J ,x,q,age_i,self.p2,self.p3))    # crea el pulsar   
            J = J + int(cuantos[i])
            j=0

    def graficar(self): 
        p=[]        
        q=[]
        l=[]
        w=[]
        b=[]
        D=[]
        sigma=[]
        for i in range(self.nTot):  # BASE 10  
            D.append(self.pulsar[i].dist)   
            p.append(np.log10(self.pulsar[i].P))
            l.append(np.log10(self.pulsar[i].L))
            w.append(np.log10(self.pulsar[i].W))
            q.append(np.log10(self.pulsar[i].Q))
            sigma.append(self.pulsar[i].sigma)
            b.append(np.log10(3.2e19*(self.pulsar[i].P*self.pulsar[i].Q)**(1./2.)))
        return p, q, l, w, b, D, sigma

    def TotalPulsares(self):
        return self.nTot
    
    def showParametros(self):
        return self.index, self.C0, self.C1, self.nTot

###########################################################################################
#               CLASE gernericTelescope AQUI ES DONDE SE VEN LOS PULSARES
###########################################################################################

class genericTelescope(object):
    def __init__(self, name, G, beta, np, deltaF, Tsys, f_obs):
        """Input: Caracteristicas fisicas del RT 
        """
        self.name = int(name)       # entero solo para id
        self.G = float(G)           # K/Jk 
        self.beta = float(beta)     # adim
        self.np = float(np)         # adim
        self.deltaF = float(deltaF) # Hz (s^-1)
        self.Tsys = float(Tsys)     # K
        self.f_obs = float(f_obs)   # MHz frecuencia de observacion nuevo para chaPulin_3.6

    def periodicDetect(self, SN, galaxy, time):
        d = galaxy.distance
        Detected = []   # pulsares detectados
        for i in range(len(galaxy.pulsar)): 
            spectral_index = galaxy.pulsar[i].spectral_index 
            tipo = galaxy.pulsar[i].dist
            if  tipo == 'D':
                Detected.append(0)
            else:            
                p = galaxy.pulsar[i].P   
                w = galaxy.pulsar[i].W
                l = galaxy.pulsar[i].L*(self.f_obs/400)**spectral_index  # correccion por banda de observacion
                if (p > w):        
                    Lmin = d**2*self.beta*SN*self.Tsys/(self.G*np.sqrt(self.np * self.deltaF \
                        * time * 3600))*np.sqrt(w/(p-w))    # ecuacion de la S/N    
                    if (Lmin < l): Detected.append(1)  
                    else: Detected.append(0)
                else: Detected.append(0)              

        return Detected

    def pulseDetect(self, SN, galaxy, time):
        d = galaxy.distance 
        giantPulses = []
        xx = 0  

        for i in range(len(galaxy.pulsar)):
            spectral_index = galaxy.pulsar[i].spectral_index 
            tipo = galaxy.pulsar[i].dist
            if  tipo == 'D':
                giantPulses.append(0)
            else:   
                p = galaxy.pulsar[i].P
                w = galaxy.pulsar[i].W
                l = galaxy.pulsar[i].L*(self.f_obs/400)**spectral_index  # correccion por banda de observacion
                Lmin = d**2*SN*self.Tsys/(self.G*np.sqrt(self.np*self.deltaF*w))    # ecuacion de la S/N 
                # mejora para liberar memoria
                if galaxy.pulsar[i].maxburst*l < Lmin:
                    giantPulses.append(0)
                else:
                    j = 0
                    nn = int(3600*time/galaxy.pulsar[i].P)
                    while xx == 0 and j < nn: 
                        if Lmin < galaxy.pulsar[i].pulseIt(): xx = 1 
                        j = j+1
                    giantPulses.append(xx)
                    xx = 0
            
        return giantPulses

    def pulsarDetect(self, SN, galaxy, time):
        Periodicos = self.periodicDetect(SN, galaxy, time)
        Pulsos = self.pulseDetect(SN, galaxy, time) 
        Total = [0]*len(Pulsos)
        for i in range(len(Pulsos)):        
            if Periodicos[i] ==1 or Pulsos[i]==1: Total[i] = 1
        p, q, l, w, b, D, sigma = galaxy.graficar()
        return p, q, l, w, b, D, sigma, Periodicos, Pulsos, Total         
  
###########################################################################################
#                       FUNCIONES UTILITARIAS
###########################################################################################
    
def pulseHist(pulsar,time,ID): # grafica el histograma de pulsos
    pulsar.observeIt(time)

    plt.clf()
    plt.suptitle(r'Tipo '+ repr(pulsar.dist)+', $L_{media}$ = ' + repr(pulsar.L)+ ' [Jy kpc$^2$]')  
    plt.hist(pulsar.pulses,30,color='black',histtype='step',lw='2') 
    #, histtype='step', histtype='bar',log='True',normed='False',log='True'
    plt.xlabel(r'$L$ [Jy kpc$^{2}$]')
    plt.ylabel(r'$N_{Pulsos}$')
    plt.xscale('linear')    
    plt.yscale('symlog') 
    plt.grid(True)
                  
    plt.savefig('out/img/pulsar_'+ repr(pulsar.name)+'_'+ID+'.png')  
    plt.clf()      
