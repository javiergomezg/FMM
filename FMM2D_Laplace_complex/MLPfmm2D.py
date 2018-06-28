#se importan las librerias necesarias
import numpy as np
from numpy import log
import math
from time import time
from matplotlib import pyplot as plt
from matplotlib import rcParams, cm
import Classdef
import implement
import inter_type
import Datacontrol
import pruebas

from matplotlib import pyplot as plt
from matplotlib import rcParams, cm

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8


###########################################################
# Se generan las particulas fuente y objetivo

#Nup = 100

#print Nup

#Ns = Nup
#Nt = Nup

#pond_source = 0.5
#Ns = int(pond_source * Nup)
#Nt = Nup - Ns

#source = Sources(Ns)
#target = Target(Nt)
#source.z[0] = 0.6 + 0.7*1j
#target.z[0] = 0.4 + 0.4*1j


#print source.z

###########################################################
# Se establece la masa que tendran las particulas fuente.
# Inicialmente sera una masa igual para cada particula.
#m = 1.0
#source.m[:] = m/Ns

#print source.z
#print 
#print abs(target.z[0]-source.z[0])
#print
#print source.m
#print 
#print target.z
#print

#fig = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
#plottarget = plt.plot(target.z.real, target.z.imag, 'ro',linewidth = 3, label='target')
#plotsources = plt.plot(source.z.real, source.z.imag, 'bo',linewidth = 3, label='target')
#plt.grid(True)
#plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#plt.axis([0,1,0,1])
#plt.show()

###########################################################
# Potencial teorico
#target_aux = Target(10)
#target_aux.z = target.z[:10]
#print target_aux.z
#t_i = time()
#phi_teorico = potencial_teorico(source,target)
#t_f = time()    
#t_teo = t_f-t_i
#print
#print 't_teo =', t_teo
#print
#print 'phi_teo'
#print phi_teorico

##############################################################
####Prueba de funciones#######################################
#zord(2)
#print
#zord(3)
#zord(4)


#ti = time()
#prueba = cord_to_index(2, 5, 3)
#ti = time()
#index_to_cord(prueba, 3)
#tf = time()
#print prueba, tf-ti

##############################################################
####Potencial FMM#############################################

#n_crit = 200
#Lx,Ly = 1.0, 1.0
#p = 6
#level_start = 4

#ti = time()
#phi_FMM,_,_,_, tpre, tcalc = potencial_FMM(source,target, n_crit, Lx, Ly, level_start, p)
#tf = time()
#print 't_fmm =' ,tf-ti, lv
#print 't_precalc = ', tpre
#print 't_calc = ', tcalc

#print 'phi_fmm'
#print phi_FMM

#error = max(abs(phi_teorico - phi_FMM))

#print
#print phi_teorico
#print
#print phi_FMM
#print
#print abs(phi_teorico-phi_FMM)
#print

#print 'error = ', error
#print 

#for i in range(n_box):
#    print i, box[i].neigh

#print target.z, target.box

#print 15/6, 15%6

#for i in range(len(box)):
#    print i, box[i].center, box[i].sources

#yL = box[8].center
#print box[8].center
#yM = box[10].center
#print box[10].center
#M = box[10].multipole
#print box[10].multipole

#print abs(box[8].center-box[10].center)


#fact = factorial_gen(p)
#print fact

#L_prueba = M2L(M,fact,yM,yL)
#print 
#print L_prueba, len(L_prueba)
#del L_prueba

##############################################################################
############Test##############################################################


###########Prueba de complejidad algoritmica###############################
#N_test = [100, 1000, 10000, 50000, 100000, 1000000]#, 10000000]###
#test_complex(N_test)


############################################################################
##########Error Vs Orden del multipolo######################################

#Nup = 10000
#l_start = 4
#p_max = 40
#p_test = np.arange(2, p_max+1, 1)
#porder_test(Nup, p_test, l_start)

###########Prueba de tiempo##########################
#print
N_test = [100, 1000, 10000, 50000, 100000, 1000000]#, 10000000]###
#test_time(N_test)

