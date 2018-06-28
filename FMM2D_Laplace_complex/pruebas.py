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

from matplotlib import pyplot as plt
from matplotlib import rcParams, cm

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8

#################Pruebas para determinar complejidad algoritmica###########################
def test_complex(N_test):
    #N_test = [100, 1000]#, 10000, 50000, 100000, 500000]
    ttest_teo, ttest_FMM = np.zeros(len(N_test)), np.zeros(len(N_test))
    error = np.zeros(len(N_test))
    ttest_FMM = np.zeros(len(N_test))
    
    l_start = 4    
    for i in range(len(N_test)):
        Ns, Nt = N_test[i], N_test[i]
        source, target = Sources(Ns), Target(Nt)
        m = 1.0
        source.m[:] = m/Ns
        n_crit = 200
        Lx,Ly = 1.0, 1.0
        p = 4
            
        t_i = time()
        phi_teo = potencial_teorico(source,target)
        t_f = time()    
        ttest_teo[i] = t_f-t_i
        
        ti = time()
        phi_FMM,_,level,_,_,_,_,_,_,_,_,_= potencial_FMM(source,target, n_crit, Lx, Ly, l_start, p)
        tf = time()
        ttest_FMM[i] = tf-ti
        
        error[i] = max(abs(phi_teo - phi_FMM))

        l_start = level        
        
        print N_test[i], ttest_teo[i], ttest_FMM[i], error[i], level
        
        del source, target
        del phi_teo, phi_FMM
        
    
    fig = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
    plottteo = plt.plot(N_test,ttest_teo , 'ko-',linewidth = 3, label = 'teorico')
    plottfmm = plt.plot(N_test,ttest_FMM , 'ko--',linewidth = 3, label = 'FMM')
    plt.grid(True)
    plt.legend(loc = 2)
    plt.title("Tiempo de calculo V/S N. particulas", fontsize = 14)
    plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('N particulas', fontsize = 12)
    plt.ylabel('t[s]', fontsize = 12)
    plt.savefig('complexB.png')
    plt.show()
    
    fig = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
    ploterr = plt.plot(N_test,error , 'ko-',linewidth = 3, label = 'teorico')
    #plottfmm = plt.plot(N_test,ttest_FMM , 'ko--',linewidth = 3, label = 'FMM')
    plt.grid(True)
    plt.legend(loc = 2)
    plt.title("Error relativo V/S N. particulas", fontsize = 14)
    plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('N particulas', fontsize = 12)
    plt.ylabel('Error ', fontsize = 12)
    plt.savefig('errorB.png')
    plt.show()

############Pruebas de tiempo de calculo FMM solo##############################################
def test_time(N_test):
    #N_test = [100, 1000]#, 10000, 50000, 100000, 500000]
    ttest_FMM = np.zeros(len(N_test))
    ttest_FMMprec = np.zeros(len(N_test))
    ttest_FMMcalc = np.zeros(len(N_test))
    t_ref = np.zeros(len(N_test), dtype = float)
    ttest_P2M = np.zeros(len(N_test))
    ttest_M2M = np.zeros(len(N_test))
    ttest_M2L = np.zeros(len(N_test))
    ttest_L2L = np.zeros(len(N_test))
    ttest_L2P = np.zeros(len(N_test))
    ttest_P2P = np.zeros(len(N_test))


    
    l_start = 4    
    for i in range(len(N_test)):
        Ns, Nt = N_test[i], N_test[i]
        source, target = Sources(Ns), Target(Nt)
        m = 1.0
        source.m[:] = m/Ns
        n_crit = 200
        Lx,Ly = 1.0, 1.0
        p = 4
            
        ti = time()
        phi_FMM,_,level,_,tpre,tcalc,TP2M, TM2M, TL2L, TM2L, TL2P, TP2P= potencial_FMM(source,target, n_crit, Lx, Ly, l_start, p)
        tf = time()
        ttest_FMM[i] = tf-ti
        ttest_FMMprec[i] = tpre
        ttest_FMMcalc[i] = tcalc
        ttest_P2M[i] = TP2M
        ttest_M2M[i] = TM2M
        ttest_M2L[i] = TM2L
        ttest_L2L[i] = TL2L
        ttest_L2P[i] = TL2P
        ttest_P2P[i] = TP2P
        
        l_start = level        
        
        print 'N_part = ',N_test[i]
        print 't_FMM = ', ttest_FMM[i]
        print 't_precalc = ',ttest_FMMprec[i]
        print 't_calc = ', ttest_FMMcalc[i]
        print 'level = ', level
        print 't_P2M = ', TP2M
        print 't_M2M = ', TM2M
        print 't_M2L = ', TM2L
        print 't_L2L = ', TL2L
        print 't_L2P = ', TL2P
        print 't_P2P = ', TP2P
        print
        
        del source, target
        del phi_FMM

    m  = ((ttest_FMM[-1]-ttest_FMM[1])/(N_test[-1]-N_test[1]))
    a, b = N_test[1], ttest_FMM[1]

    for i in range(len(t_ref)):
        t_ref[i] = m*(N_test[i] - a) + b 
        
        
    
    fig = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
    plottfmm = plt.loglog(N_test,ttest_FMM , 'ko--',linewidth = 3, label = 'FMM')
    plottref = plt.loglog(N_test,t_ref , 'k+-',linewidth = 3, label = 'ref')
    plt.grid(True)
    plt.legend(loc = 2)
    plt.title("Tiempo de calculo V/S N. particulas", fontsize = 14)
    plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    #plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel('N particulas', fontsize = 12)
    plt.ylabel('t[s]', fontsize = 12)
    plt.savefig('timecalc.png')
    plt.show()
    
    fig = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
    plottfmm = plt.loglog(N_test,ttest_FMM , 'ko--',linewidth = 3, label = 'FMM_total')
    plottref = plt.loglog(N_test,ttest_FMMprec , 'b+-',linewidth = 3, label = 'precalculo')
    plottref = plt.loglog(N_test,ttest_FMMcalc , 'r+-',linewidth = 3, label = 'calculo')
    plt.grid(True)
    plt.legend(loc = 2)
    plt.title("Tiempo de calculo V/S N. particulas", fontsize = 14)
    plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('N particulas', fontsize = 12)
    plt.ylabel('t[s]', fontsize = 12)
    plt.savefig('timedetalle.png')
    plt.show()
        
    fig = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
    plottfmm = plt.loglog(N_test,ttest_FMM , 'ko--',linewidth = 3, label = 'FMM_total')
    plottref = plt.loglog(N_test,ttest_P2M , 'b+-',linewidth = 3, label = 'P2M')
    plottref = plt.loglog(N_test,ttest_M2M , 'r+-',linewidth = 3, label = 'M2M')
    plottref = plt.loglog(N_test,ttest_M2L , 'y+-',linewidth = 3, label = 'M2L')
    plottref = plt.loglog(N_test,ttest_L2L , 'g+-',linewidth = 3, label = 'L2L')
    plottref = plt.loglog(N_test,ttest_L2P , 'r*--',linewidth = 3, label = 'L2P')
    plottref = plt.loglog(N_test,ttest_P2P , 'b*--',linewidth = 3, label = 'P2P')
    plt.grid(True)
    plt.legend(loc = 2)
    plt.title("Tiempo de calculo V/S N. particulas", fontsize = 14)
    plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('N particulas', fontsize = 12)
    plt.ylabel('t[s]', fontsize = 12)
    plt.savefig('tinteract.png')
    plt.show()
        

########################Pruebas para determinar el error respecto al orden de P################
def porder_test(Nup, p_test, l_start):    
    Ns, Nt = Nup, Nup
    source, target = Sources(Ns), Target(Nt)
    m = 1.0
    source.m[:] = m/Ns
    n_crit = 200
    Lx,Ly = 1.0, 1.0
    #l_start = 5
    ttest_FMM = np.zeros(len(p_test))
    error = np.zeros(len(p_test))
    
    #print p_test
    t_i = time()
    phi_teo = potencial_teorico(source,target)
    t_f = time()    
    ttest_teo = t_f-t_i
    
    for i in range(len(p_test)):
        #t_i = time()
        #phi_teo = potencial_teorico(source,target)
        #t_f = time()    
        #ttest_teo[i] = t_f-t_i
        
        ti = time()
        phi_FMM,_,_,_,_,_,_,_,_,_,_,_ = potencial_FMM(source,target, n_crit, Lx, Ly, l_start, p_test[i])
        tf = time()
        ttest_FMM[i] = tf-ti
        
        error[i] = max(abs(phi_teo - phi_FMM))
        
        print p_test[i], ttest_FMM[i], error[i]
        
        del phi_FMM
     
    fig = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
    ploterr = plt.plot(p_test,error , 'ko-',linewidth = 3, label = 'teorico')
    #plottfmm = plt.plot(N_test,ttest_FMM , 'ko--',linewidth = 3, label = 'FMM')
    plt.grid(True)
    plt.legend(loc = 2)
    plt.title("Error relativo V/S Orden p del multipolo", fontsize = 14)
    plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    #plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('orden p del multipolo', fontsize = 12)
    plt.ylabel('Error ', fontsize = 12)
    plt.savefig('errorp.png')
    plt.show()
    
    fig2 = plt.figure(num = None, figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
    plottime = plt.plot(p_test,ttest_FMM , 'ko-',linewidth = 3, label = 'teorico')
    #plottfmm = plt.plot(N_test,ttest_FMM , 'ko--',linewidth = 3, label = 'FMM')
    plt.grid(True)
    plt.legend(loc = 2)
    plt.title("Tiempo de calculo V/S Orden p del multipolo", fontsize = 14)
    plt.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    #plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('orden p del multipolo', fontsize = 12)
    plt.ylabel('Tiempo [s] ', fontsize = 12)
    plt.savefig('timep.png')
    plt.show()
