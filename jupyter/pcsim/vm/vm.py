import ROOT
import os

## import the utility library
import pcsim.util
from pcsim.util import tf1_wrap

## compile the associated root macro that contains the actual VM info if needed
ROOT.gROOT.ProcessLine('.L %s.cc+' % os.path.splitext(__file__)[0])

## wrap around the TF1s for more intuitive usage
class epsilon_Q2(tf1_wrap):
    def factory(self, lim): 
        return ROOT.vm.epsilon_Q2(lim[0], lim[1])
    def pardef(self):
        return {'W': 0, 'Ebeam': 1}
class dsigma_dt_vm_brodsky(tf1_wrap):
    def factory(self, lim): 
        return ROOT.vm.dsigma_dt_vm_brodsky(lim[0], lim[1])
    def pardef(self):
        return {'c2g': 0, 'c3g': 1, 'b': 2, 'Mv': 3, 'W': 4}
class dsigma_dexp_bt_vm_brodsky(tf1_wrap):
    def factory(self, lim): 
        return ROOT.vm.dsigma_dexp_bt_vm_brodsky(lim[0], lim[1])
    def pardef(self):
        return {'c2g': 0, 'c3g': 1, 'b': 2, 'Mv': 3, 'W': 4}
class sigma_vm_brodsky_W(tf1_wrap):
    def factory(self, lim): 
        return ROOT.vm.sigma_vm_brodsky(lim[0], lim[1])
    def pardef(self):
        return {'c2g': 0, 'c3g': 1, 'b': 2, 'Mv': 3}
class R_vm_martynov(tf1_wrap):
    def factory(self, lim):
        return ROOT.vm.R_vm_martynov(lim[0], lim[1])
    def pardef(self):
        return {'Mv': 0, 'R_c': 1, 'R_n': 2}
class dipole_ff_vm(tf1_wrap):
    def factory(self, lim):
        return ROOT.vm.dipole_ff_vm(lim[0], lim[1])
    def pardef(self):
        return {'Mv': 0, 'dipole_n': 1}
class sigma_vm_brodsky_W_Q2(tf1_wrap):
    def factory(self, lim):
        return ROOT.vm.sigma_vm_brodsky_W_Q2(lim[0], lim[1])
    def pardef(self):
        return {'c2g': 0,
                'c3g': 1,
                'b'  : 2,
                'Mv' : 3,
                'Q2' : 4,
                'W'  : 5,
                'Ebeam': 6,
                'dipole_n': 7,
                'R_c': 8,
                'R_n': 9}
class dsigma_dt_vm_brodsky_Q2W(tf1_wrap):
    def factory(self, lim):
        return ROOT.vm.sigma_vm_brodsky_W_Q2(lim[0], lim[1])
    def pardef(self):
        return {'c2g': 0,
                'c3g': 1,
                'b'  : 2,
                'Mv' : 3,
                'Q2' : 4,
                'W'  : 5,
                'Ebeam': 6,
                'dipole_n': 7,
                'R_c': 8,
                'R_n': 9}
class r00_04(tf1_wrap):
    def factory(self, lim):
        return ROOT.vm.r00_04(lim[0], lim[1])
    def pardef(self):
        return {'W': 0, 'Mv': 1, 'Ebeam': 2, 'R_c': 3, 'R_n': 4}
class ctheta_schc(tf1_wrap):
    def factory(self, lim):
        return ROOT.vm.ctheta_cshc(lim[0], lim[1])
    def pardef(self):
        return {'Q2': 0, 'W': 1, 'Mv': 2, 'Ebeam': 3, 'R_n': 4, 'R_c': 5}
class dsigma_dctheta_vm_brodsky_Q2W(tf1_wrap):
    def factory(self, lim):
        return ROOT.vm.dsigma_dctheta_vm_brodsky_Q2W(lim[0], lim[1])
    def pardef(self):
        return {'c2g': 0,
                'c3g': 1,
                'b'  : 2,
                'Mv' : 3,
                'Q2' : 4,
                'W'  : 5,
                'Ebeam': 6,
                'dipole_n': 7,
                'R_c': 8,
                'R_n': 9}
