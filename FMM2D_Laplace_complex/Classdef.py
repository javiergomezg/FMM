import numpy as np

class Sources():
    def __init__(self,Ns, Lx = 1.0, Ly = 1.0):
        self.z = Lx*np.random.rand(Ns) + Ly*np.random.rand(Ns)*1j
        self.m = np.zeros(Ns, dtype = float)

class Target():
    def __init__(self, Nt, Lx = 1.0, Ly = 1.0):
        self.z = Lx*np.random.rand(Nt) + Ly*np.random.rand(Nt)*1j
        self.box = np.zeros(Nt, dtype = int)

class Celula():
    def __init__(self, lv, p = 1):
        self.nivel = lv
        self.parent = 0
        self.child = np.zeros(4, dtype = int)
        self.center = 0. + 0.*1j
        self.neigh = np.zeros(0, dtype = int)
        self.sources = np.empty(0, dtype = complex)
        self.sourcesm = np.empty(0, dtype = float)
        self.multipole = np.zeros(p+1, dtype = complex)
        self.local = np.zeros(p+1, dtype = complex)