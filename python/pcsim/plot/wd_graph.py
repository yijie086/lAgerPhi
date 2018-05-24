import ROOT
import numpy as np

class wd_graph:
    options = {
        'color': lambda s, val : s.config(mcolor=val, lcolor=val),
        'mcolor': lambda s, val : s.graph.SetMarkerColor(val),
        'mstyle': lambda s, val : s.graph.SetMarkerStyle(val),
        'lcolor': lambda s, val : s.graph.SetLineColor(val),
        'lwidth': lambda s, val : s.graph.SetLineWidth(val),
        'name': lambda s, val : setattr(s, 'name', val),
        'msize': lambda s, val : s.graph.SetMarkerSize(val)
    }

    def __init__(self, x, y, ex=None, ey=None, **kwargs):
        '''Create a graph from existing data.

        parameters: 
            * x: array of x values
            * y: array of y values
            * ex: x errors, can be an array, or {min: [], max: []} dict (Opt., D: None)
            * ey: y errors, can be an array, or {min: [], max: []} dict (Opt., D: None)
        kwargs: additional options, see config() for mode info
        '''
        self.x = np.array(x)
        self.y = np.array(y)
        self.ex = self.init_errors(ex)
        self.ey = self.init_errors(ey)
        self.graph = ROOT.TGraphAsymmErrors(
            len(self.x),
            self.x, self.y, 
            np.array(self.ex['min']), np.array(self.ex['max']), 
            np.array(self.ey['min']), np.array(self.ey['max']))
        self.config(**kwargs)
    
    def init_errors(self, e):
        if e is None:
            return {'min': np.zeros(len(self.x)), 'max': np.zeros(len(self.x))}
        elif not isinstance(e, dict):
            return {'min': np.array(e), 'max': np.array(e)}
        return e


    def config(self, **kwargs):
        '''Configure the underlying histogram.

        See wd_graph.options for the available options.
        '''
        for val in kwargs:
            if val in wd_graph.options:
                wd_graph.options[val](self, kwargs[val])
            else:
                print("%s: unknown option: %s" % (__name__, val))
