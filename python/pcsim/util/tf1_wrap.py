import ROOT, pcsim.util
ROOT.gROOT.ProcessLine('#include <pcsim/core/interval.hh>')

range_type = ROOT.pcsim.interval(ROOT.double)

## wrap around the TF1s for more intuitive usage
class tf1_wrap:
    '''A simple tf1 wrapper that allows for a more pythonic TF1 construction.'''
    def __init__(self, limits, **kwargs):
        '''Create a function object with pre-defined limits.

        limits: a min-max pair
        kwargs: Optional arguments are set as function parameters if recognized.
        '''
        self.f = self.factory(limits)
        self.set_options(**kwargs)

    def factory():
        '''(dummy) Returns a TF1.'''
        return None
    def pardef():
        '''(dummy) Return a dictionary of the allowed arguments.'''

    def set_options(self, **kwargs):
        '''General function to configure the underlying TF1.'''
        known_pars = self.pardef()
        for par in kwargs:
            if not par in known_pars:
                raise KeyError('%s: Unkown parameter name: ' % (self.__class__, key))
            val = kwargs[par]
            self.f.SetParameter(known_pars[par], val)
