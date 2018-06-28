import numpy as np
import inter_type

def zord(ls = 4):
    nx, ny = 2**ls, 2**ls
    nl = nx*ny
    #print nl
    #print
    aux  = np.zeros([nx, ny], dtype=int)
    #print aux
    #print
    cont = ((4**ls)-1)/3 
    #print cont
    for i in range(nx):
        for j in range(ny):
            aux[i,j] = cord_to_index(i,j,ls)
            
            #print i, j, ib[2:], jb[2:], nb 
    print aux
    
    
def cord_to_index(xi,yj,lv):
    ib, jb = bin(xi), bin(yj)
    nb = max(len(ib)-2, len(jb)-2)
    while(abs(len(ib)-len(jb)) > 0):
        if len(ib)-2 < nb:
            ib = '0b'+'0'+ib[2:]
        if len(jb)-2 < nb:
            jb = '0b'+'0'+jb[2:]
    auxbin = ib[2]+jb[2]
    for k in range(1,nb):
        auxbin = auxbin + ib[2+k]+jb[2+k]
    #print auxbin
    return int(auxbin, 2) + ((4**lv)-1)/3
    
def index_to_cord(index, lv):
    indaux = index - ((4**lv)-1)/3
    indexbin = bin(indaux)
    if (len(indexbin) - 2)%2 > 0:
        indexbin = '0b'+'0'+indexbin[2:]
    nb = (len(indexbin) - 2)/2
    #print indexbin, nb
    ixb , jyb = indexbin[2], indexbin[3]
    for i in range(1,nb):
        ixb , jyb = ixb + indexbin[2 + 2*i], jyb + indexbin[2*i + 3]
    ix, jy = int(ixb, 2), int(jyb, 2)
    #print ix, jy
    return ix, jy
    

def level_max(S, n_crit = 10, Lx = 1.0, Ly = 1.0, ls = 4):
    n_p, level, acepted = len(S), ls - 1, False
    while acepted ==  False:
        level +=1
        n_box = 4**level
        dx, dy = Lx/(2**level), Ly/(2**level)
        n_b_part = np.zeros(n_box, dtype = int)
        acepted = True
        for i in range(n_p):
            kx, ky = int(S[i].real/dx), int(S[i].imag/dy)
            k = cord_to_index(ky, kx, level) - ((4**level)-1)/3
            n_b_part[k] +=1
            if n_b_part[k] > n_crit:
                acepted = False
                break
    #print n_b_part
    return n_b_part, level

def box_make(level=4, p = 1):
    n_box = (4**(level+1) - 1)/3
    box = np.empty(n_box, dtype = object)
    lv = 0
    j = 0
    for i in range(n_box):
        box[i] = Celula(lv, p)
        j+=1
        if j == 4**lv:
            lv += 1
            j = 0
    return box, n_box

def box_center(b, n_box, Lx = 1.0, Ly = 1.0):
    box = np.copy(b)
    for i in range(n_box):
        l_aux = box[i].nivel
        dx, dy = Lx/(2**l_aux), Ly/(2**l_aux)
        jy, ix = index_to_cord(i, l_aux)
        box[i].center = 0.5*(2*ix + 1)*dx + 0.5*(2*jy + 1)*dy*1j 
        #print i, l_aux, box[i].center
    return box

def parent_child(b, n_box,lv,Lx = 1.0, Ly = 1.0):
    box = np.copy(b)
    for i in range(n_box - 4**lv):
        l_aux = box[i].nivel
        dx, dy = Lx/(2**l_aux), Ly/(2**l_aux)
        #print i
        #print l_aux, dx, dy
        #print box[i].center
        #print box[i].center + 0.5*(dx + dy*1j)
        
        j0x = int((box[i].center.real - 0.25*dx)/(0.5*dx))
        i0y = int((box[i].center.imag - 0.25*dy)/(0.5*dy))
        #print i0x, j0y
        i0 = cord_to_index(i0y,j0x,l_aux + 1)
        
        j1x = int((box[i].center.real + 0.25*dx)/(0.5*dx))
        i1y = int((box[i].center.imag - 0.25*dy)/(0.5*dy))
        #print  i1x, j1y
        i1 = cord_to_index(i1y,j1x,l_aux + 1)
        
        j2x = int((box[i].center.real - 0.25*dx)/(0.5*dx))
        i2y = int((box[i].center.imag + 0.25*dy)/(0.5*dy))
        #print  i2x, j2y
        i2 = cord_to_index(i2y,j2x,l_aux + 1)
        
        j3x = int((box[i].center.real + 0.25*dx)/(0.5*dx))
        i3y = int((box[i].center.imag + 0.25*dy)/(0.5*dy))
        #print  i3x, j3y
        i3 = cord_to_index(i3y,j3x,l_aux + 1)
        
        box[i].child = [i0, i1, i2, i3]

        box[i0].parent,box[i1].parent,box[i2].parent,box[i3].parent=i,i,i,i
        
        #print i, box[i].child, box[i].parent
    return box
    
def clustering(b, nbp, lv, source, Lx=1.0, Ly=1.0):
    box = np.copy(b)
    cont = (4**lv - 1)/3
    for i in range(len(nbp)):
        box[i + cont].sources = np.zeros(nbp[i], dtype=complex)
        box[i + cont].sourcesm = np.zeros(nbp[i], dtype=float)
        #print i+cont, nbp[i], box[i + cont].sources
    j_aux = np.zeros(len(nbp), dtype=int)
    dx, dy = Lx/(2**lv), Ly/(2**lv)
    for i in range(len(source.z)):
        jx, iy = int(source.z[i].real/dx), int(source.z[i].imag/dy)
        k = cord_to_index(iy, jx, lv)
        box[k].sources[j_aux[k - cont]] = source.z[i]
        box[k].sourcesm[j_aux[k - cont]] = source.m[i]
        j_aux[k - cont] += 1
    return box
    
def box_target(t, lv, Lx=1.0, Ly=1.0):
    target = t
    dx, dy = Lx/(2**lv), Ly/(2**lv)
    for i in range(len(target.z)):
        jx, iy = int(target.z[i].real/dx), int(target.z[i].imag/dy)
        target.box[i] = cord_to_index(iy, jx, lv)
    return target
    
def factorial_gen(p):
    factorial = np.zeros(p + 1, dtype=float)
    factorial[0] = 1

    for i in range(1,p+1):
        factorial[i] = i*factorial[i-1]

    return factorial
    
def multipole_calc(b, lv, fact):
    box = np.copy(b)
    cont = (4**lv - 1)/3
    p = len(fact)
    #print fact
    ti = time()
    for i in range(cont, len(box)):
        #print i
        box[i].multipole = P2M(box[i].sourcesm, fact,
                            box[i].sources, box[i].center)
        #print i, box[i].multipole, len(box[i].sources)
        #print box[i].sources
        #print box[i].sourcesm
        #print box[i].center
    tf = time()
    TP2M = tf - ti
        
    ti = time()
    for i in range(cont-1, 4 , -1):
        #print i, box[i].child
        box[i].multipole = np.zeros(p+1)
        Mc0 = box[box[i].child[0]].multipole
        #print Mc0
        Mc1 = box[box[i].child[1]].multipole
        #print Mc1
        Mc2 = box[box[i].child[2]].multipole
        #print Mc2
        Mc3 = box[box[i].child[3]].multipole
        #print Mc3
        
        
        M0 = M2M(Mc0, fact, box[i].center,box[box[i].child[0]].center)
        #print M0
        M1 = M2M(Mc1, fact, box[i].center,box[box[i].child[1]].center)
        #print M1
        M2 = M2M(Mc2, fact, box[i].center,box[box[i].child[2]].center)
        #print M2
        M3 = M2M(Mc3, fact, box[i].center,box[box[i].child[3]].center)
        #print M3
        
        box[i].multipole = M0+M1+M2+M3
        #print i, box[i].multipole
    tf = time()
    TM2M = tf - ti
        
    return box, TP2M, TM2M
    
    

def neigh_loc(b, n_box):
    box = np.copy(b)
    for i in range(1, n_box):
        lv_aux = box[i].nivel
        ipt, jpt = index_to_cord(i, lv_aux)
        
        ingh = np.empty(8, dtype=int)
        jngh = np.empty(8, dtype=int)
        
        ingh[0], jngh[0] = ipt - 1, jpt
        ingh[1], jngh[1] = ipt - 1, jpt + 1
        ingh[2], jngh[2] = ipt    , jpt + 1
        ingh[3], jngh[3] = ipt + 1, jpt + 1
        ingh[4], jngh[4] = ipt + 1, jpt
        ingh[5], jngh[5] = ipt + 1, jpt - 1
        ingh[6], jngh[6] = ipt    , jpt - 1
        ingh[7], jngh[7] = ipt - 1, jpt - 1
        
        #print i, ingh, jngh, 2**lv_aux
        
        n_neigh = 0
        
        for j in range(8):
            if ingh[j] >= 0 and jngh[j] >= 0:
                if ingh[j] < 2**lv_aux and jngh[j] < 2**lv_aux:
                    n_neigh +=1
        aux_neigh = np.zeros(n_neigh, dtype=int)
        #print i, n_neigh
        
        l = 0
        for j in range(8):
            if ingh[j] >= 0 and jngh[j] >= 0:
                if ingh[j] < 2**lv_aux and jngh[j] < 2**lv_aux:
                    aux_neigh[l] = cord_to_index(ingh[j],jngh[j],lv_aux)
                    l +=1
        box[i].neigh = aux_neigh
        #print i, box[i].neigh
    return box

def local_calc(b, lv, n_box, fact):
    p = len(fact)
    box = np.copy(b)
    #stop = (4**lv - 1)/3
    TL2L = 0.0
    TM2L = 0.0

    for i in range(5, n_box):
        parent = box[i].parent
        neighp = box[parent].neigh
        #print i, parent, neighp
        interM2L = np.zeros(len(neighp)*4, dtype = int)
        L_aux = np.zeros(p, dtype=complex)
        #print i 
        #print interM2L
        #print box[parent].local
        ti = time()
        L_aux += L2L(box[parent].local, fact, box[parent].center, box[i].center)
        tf = time()
        TL2L_aux = tf - ti
        TL2L += TL2L_aux
        #print i, parent
        #print box[parent].local
        #print L_aux
        ti = time()
        for j in range(len(neighp)):
            #print interM2L[j:j+4]
            #print neighp[j]
            #print box[neighp[j]].child
            interM2L[4*j:4*j+4] = box[neighp[j]].child
        #print i, interM2L
        for j in range(len(interM2L)):
            cond = True
            for k in range(len(box[i].neigh)):
                if interM2L[j] - box[i].neigh[k] == 0:
                    cond = False
            if cond == True:
                #print i, interM2L[j]
                L_aux += M2L(box[interM2L[j]].multipole,fact,box[interM2L[j]].center,box[i].center)
                
        box[i].local = L_aux
        tf = time()
        TM2L_aux = tf - ti
        TM2L += TM2L_aux
        #print i, box[i].local
    del L_aux, cond, interM2L, p, neighp, parent
    return box, TL2L, TM2L
    
def phi_calc(T, b, fact):
    target = T
    box = b
    Nt = len(target.z)
    phi_aprox = np.zeros(Nt, dtype = complex)    
    TL2P = 0.
    TP2P = 0.
    for i in range(Nt):
        ti = time()
        phi_aprox[i] = L2P(box[target.box[i]].local,fact, target.z[i],box[target.box[i]].center)
        tf = time()
        TL2P_aux = tf - ti
        TL2P += TL2P_aux
        
        ti = time()
        phi_aprox[i] += P2P(box[target.box[i]].sourcesm,box[target.box[i]].sources, target.z[i])
        neigh = box[target.box[i]].neigh
        for j in range(len(neigh)):
            phi_aprox[i] += P2P(box[neigh[j]].sourcesm,box[neigh[j]].sources, target.z[i])
        tf = time()
        TP2P_aux = tf - ti
        TP2P += TP2P_aux
    #print phi_aprox
    return phi_aprox, TL2P, TP2P

    

        
        