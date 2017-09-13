import ROOT
import numpy

class wd_graph:
    options = {
        'color': lambda (s, val) : s.config(mcolor=val, lcolor=val),
        'mcolor': lambda (s, val) : s.graph.SetMarkerColor(val),
        'lcolor': lambda (s, val) : s.graph.SetLineColor(val),
        'lwidth': lambda (s, val) : s.graph.SetLineWidth(val),
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
        self.ex = init_errors(ex)
        self.ey = init_errors(ey)
        self.graph = ROOT.TGraphAsymmErrors(
            len(self.x),
            self.x, self.y, 
            self.ex['min'], self.ex['max'], 
            self.ey['min'], self.ey['max'])
        config(**kwargs)
    
    def init_errors(self, e):
        if e is None:
            return {'min': zeros(len(self.x)), 'max': zeros(len(self.x))}
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
                print "%s: unknown option: %s" % (__name__, val)
