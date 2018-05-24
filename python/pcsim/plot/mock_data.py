import ROOT

class mock_data:
    options = {
        'xtitle': lambda s_val : s_val[0].graph.GetXaxis().SetTitle(s_val[1]),
        'ytitle': lambda s_val1 : s_val1[0].graph.GetYaxis().SetTitle(s_val1[1]),
        'color': lambda s_val2 : s_val2[0].config(mcolor=s_val2[1], lcolor=s_val2[1]),
        'mcolor': lambda s_val3 : s_val3[0].graph.SetMarkerColor(s_val3[1]),
        'lcolor': lambda s_val4 : s_val4[0].graph.SetLineColor(s_val4[1]),
        'lwidth': lambda s_val5 : s_val5[0].graph.SetLineWidth(s_val5[1]),
        'title': lambda s_val6 : s_val6[0].graph.SetTitle(s_val6[1])
    }
    rng = ROOT.TRandom3(0)

    def __init__(self, name, hacc, ftheo, scale=1., random=False, 
                 offset=1., minevt = 1.1, **kwargs):
        '''Create mock data from a histogram of accepted events.

        parameters: 
            * name: the histogram name
            * hacc: histogram of accepted events
            * ftheo: theory TF1 that describes our mock data
            * scale: scale factor to be applied to hacc (D: 1.)
            * random: use Poisson statistics for realistic mock data (D: False)
            * offset: offset factor to be multiplied with the mock data (D: 1.)
            * minevt: minimum number required to be considered a measurement (D: 1.1)
        kwargs: additional options, see config() for mode info
        '''
        self.name = name
        self.hacc = hacc.Clone()
        self.hacc.SetName('%s-acc' % name)
        self.graph = hacc.Clone()
        self.graph.SetName(name)
        self.config(**kwargs)
        ## calculate the actual mock data
        for i in range(1, hacc.GetNbinsX() + 1):
            x = hacc.GetBinCenter(i)
            y = ftheo.Eval(x)
            dy = 0.
            N = hacc.GetBinContent(i) * scale
            Nold = N
            if random:
                N = mock_data.rng.Poisson(N)
            ## do we have a reasonable amount of events
            if (N > 1.1):
                y *= N / Nold
                dy = y / ROOT.TMath.Sqrt(N)
                y *= offset
            else:
                y = 0.
                dy = 0.
                N = 0.
            self.graph.SetBinContent(i, y)
            self.graph.SetBinError(i, dy)
            self.hacc.SetBinContent(i, N)

    def config(self, **kwargs):
        '''Configure the underlying histogram.

        See mock_data.options for the available options.
        '''
        for val in kwargs:
            if val in mock_data.options:
                mock_data.options[val]((self, kwargs[val]))
            else:
                print("%s: unknown option: %s" % (__name__, val))
