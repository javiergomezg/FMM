import numpy as np
import time
import Classdef
import Datacontrol
import inter_type


def potencial_teorico(source, target):
    Nt = len(target.z)
    phi_teorico = np.zeros(Nt, dtype = float)
    for i in range(Nt):
        phi_teorico[i] = P2P(source.m, source.z, target.z[i])
    return phi_teorico


def potencial_FMM(source, target, n_crit = 10, Lx = 1.0, Ly = 1.0, lv_start = 4, p = 4):
    ti = time()
    nbp, level = level_max(source.z, n_crit, Lx, Ly, lv_start)
    box, n_box = box_make(level, p)
    box = box_center(box, n_box, Lx, Ly)
    box = parent_child(box, n_box, level, Lx, Ly)
    box = clustering(box, nbp, level, source, Lx, Ly)
    target = box_target(target, level, Lx, Ly)
    factorial = factorial_gen(p)
    box = neigh_loc(box, n_box)
    tf = time()
    t_precalc = tf-ti    
    #print factorial
    ti = time()
    box, TP2M, TM2M = multipole_calc(box, level, factorial)
    box, TL2L, TM2L = local_calc(box,level, n_box, factorial)
    phi_aprox, TL2P, TP2P = phi_calc(target, box, factorial)
    phi_FMM = phi_aprox.real
    tf = time()
    t_calc = tf-ti    
    #print factorial
    return phi_FMM, n_box, level, box, t_precalc, t_calc, TP2M, TM2M, TL2L, TM2L, TL2P, TP2P