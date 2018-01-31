import ROOT, pcsim.util

## load the TF1 factories
ROOT.gROOT.ProcessLine('#include <pcsim/root/physics/vm.hh>')

vm_impl = ROOT.pcsim.root.physics.vm

## import the utility tf1 wrapper
from pcsim.util import range_type, tf1_wrap

## wrap around the TF1s for more intuitive usage
class epsilon_Q2(tf1_wrap):
    def factory(self, lim): 
        return vm_impl.epsilon_Q2(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'W': 0, 'Ebeam': 1}
class dsigma_dt_vm_brodsky(tf1_wrap):
    def factory(self, lim): 
        return vm_impl.dsigma_dt_vm_brodsky(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'c2g': 0, 'c3g': 1, 'b': 2, 'Mv': 3, 'W': 4}
class dsigma_dexp_bt_vm_brodsky(tf1_wrap):
    def factory(self, lim): 
        return vm_impl.dsigma_dexp_bt_vm_brodsky(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'c2g': 0, 'c3g': 1, 'b': 2, 'Mv': 3, 'W': 4}
class sigma_vm_brodsky_W(tf1_wrap):
    def factory(self, lim):
        return vm_impl.sigma_vm_brodsky_W(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'c2g' : 0, 'c3g' : 1, 'b' : 2, 'Mv' : 3}
class R_vm_martynov(tf1_wrap):
    def factory(self, lim):
        return vm_impl.R_vm_martynov(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'Mv': 0, 'R_c': 1, 'R_n': 2}
class dipole_ff_vm(tf1_wrap):
    def factory(self, lim):
        return vm_impl.dipole_ff_vm(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'Mv': 0, 'dipole_n': 1}
class sigma_vm_brodsky_W_Q2(tf1_wrap):
    def factory(self, lim):
        return vm_impl.sigma_vm_brodsky_W_Q2(range_type(lim[0], lim[1]))
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
        return vm_impl.sigma_vm_brodsky_W_Q2(range_type(lim[0], lim[1]))
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
        return vm_impl.r00_04(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'W': 0, 'Mv': 1, 'Ebeam': 2, 'R_c': 3, 'R_n': 4}
class ctheta_schc(tf1_wrap):
    def factory(self, lim):
        return vm_impl.ctheta_cshc(range_type(lim[0], lim[1]))
    def pardef(self):
        return {'Q2': 0, 'W': 1, 'Mv': 2, 'Ebeam': 3, 'R_n': 4, 'R_c': 5}
class dsigma_dctheta_vm_brodsky_Q2W(tf1_wrap):
    def factory(self, lim):
        return vm_impl.dsigma_dctheta_vm_brodsky_Q2W(range_type(lim[0], lim[1]))
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
