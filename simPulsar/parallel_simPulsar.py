#!/usr/bin/env python
""" Tercera implementacion del codigo
    Tercera modularizacion
    Primera paralelizacion
    Datos obtenidos de atnfParameters
    Primer paralelizacion completa.
    chaPulin8.8
"""

import math, sys, time
import numpy as np
from scipy import array
import simPulsar 
import pp
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})

ppservers = ()
#ppservers = ("10.0.0.1",)

###########################################################################################
#                                   rutinaInicial()
###########################################################################################

def rutinaInicial():
    resumen = 'out/out_simPulsar.dat'
    ss = file(resumen, "w") # archivo de resumen
    print >>ss, "# Archivo\tPeriodic\tPulsos\tAmbos"
    ss.close()
    return "El archivo " + repr(resumen) + " ha sido creado."

    resumengal = 'sim_parametros.gal'
    sa = file(resumengal, "w") # archivo de resumen
    print >>sa, "#\\index\tC0\tC1\tN_p"
    sa.close()
    return "El archivo " + repr(resumen) + " ha sido creado."

###########################################################################################
#                                   Iteracion()
###########################################################################################

def Iteracion():
    return simPulsar.main()

###########################################################################################
#                                   PARALELIZACION
###########################################################################################

if len(sys.argv) > 1:
    ncpus = int(sys.argv[1])
    job_server = pp.Server(ncpus, ppservers=ppservers)
else:
    job_server = pp.Server(ppservers=ppservers)

print "Iniciando pp con", job_server.get_ncpus(), "esclavos"
inicializar = job_server.submit(rutinaInicial)
print inicializar()

start_time = time.time()

inputs = 1000
jobs = []
for i in range(inputs):
    jobs.append(job_server.submit(Iteracion,(),(),("simPulsar",),group='sim'))

job_server.wait('sim')

###########################################################################################
#          AQUI VOY A CALCULAR LOS PROMEDIOS Y DESVIACIONES
###########################################################################################

A = [0.]*inputs
B = [0.]*inputs 
status = 'VOID' 
totA = array([0.]*4)
totB = array([0.]*4)
dA = array([0.]*4)
dB = array([0.]*4)

for i in range(inputs):
    if jobs[i]() is not None: status, A[i], B[i] = jobs[i]()
    else: status, A[i], B[i] = 'no', [0.]*4, [0.]*4
    print 'Tarea ', i, status
    totA = totA + A[i] 
    totB = totB + B[i] 

for j in range(4):
    for i in range(inputs):
        dA[j] = dA[j] + (totA[j]/inputs - A[i][j])**2
        dB[j] = dB[j] + (totB[j]/inputs - B[i][j])**2   

devA = np.sqrt(array(dA)/(inputs-1))
devB = np.sqrt(array(dB)/(inputs-1))

totA = array(totA)
totB = array(totB)

print "Tiempo de ejecucion para ", inputs, " tareas: ", time.time() - start_time, "s"

###########################################################################################
#          AQUI VOY A GUARDAR LOS DATOS DE RESUMEN EN EL ARCHIVO DE RESUMEN
###########################################################################################
  
labels = ['normal', 'log-norm', 'p (<2)', 'p (>2)']

resumen = 'out/out_simPulsar.dat'
ss = file(resumen, "a") # archivo de resumen
print >>ss, "-->> Estadisticas <<--"

print >>ss, "Busqueda por periodicidad", totA.sum()/inputs, "(Gran total)"
totA = totA/inputs
print >>ss, "Promedios:"
for i in range(4): print >>ss, "\t", labels[i], "\t", totA[i] ,"\pm", devA[i]

print >>ss, "Busqueda por pulsos", totB.sum()/inputs, "(Gran total)"
totB = totB/inputs
print >>ss, "Promedios:"
for i in range(4): print >>ss, "\t", labels[i], "\t", totB[i] ,"\pm", devB[i]

ss.close()

###########################################################################################
#          AQUI VOY A GRAFICAR EN BARRAS LOS RESULTADOS
###########################################################################################

ind = np.arange(4)
maxj = max(totA+totB)
width = 0.35

p1 = plt.bar(ind, totA, width, color='red', yerr=devA, alpha=1,error_kw=dict(elinewidth=1,ecolor='black'))
p2 = plt.bar(ind+width, totB, width, color='orange', yerr=devB, alpha=1,error_kw=dict(elinewidth=1,ecolor='black')) #, bottom=totA,)

plt.ylabel(r'$N_{p}$')
#plt.title(r"Detectabilidad seg\'un distribuci\'on de pulsos")
plt.xticks(ind+width/2., (r'Normal',r'Log-normal',r'Potencia ($\alpha < 2.0$)',r'Potencia ($\alpha > 2.0$)') )
#plt.yticks(np.arange(0,maxj+10,10))
plt.legend( (p1[0], p2[0]), ('Periodicidad', 'Pulsos') )
plt.grid(True)
plt.savefig('out/barras_detectados.pdf')
plt.clf()

###########################################################################################
#                           AQUI VOY A TERMINAR
###########################################################################################

job_server.print_stats()
