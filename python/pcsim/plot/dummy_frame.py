import ROOT

class dummy_frame:
    dummy_cnt = 0
    options = {
        'xtitle': lambda s, val: s.graph.GetXaxis().SetTitle(val),
        'ytitle': lambda s, val: s.graph.GetYaxis().SetTitle(val),
        'ylim': lambda s, val: s.graph.GetYaxis().SetRangeUser(*val),
        'title': lambda s, val: s.graph.SetTitle(val),
        'logx': lambda s, val: s.set_logx(val),
        'logy': lambda s, val: s.set_logy(val)}
    def __init__(self, title='', href=None, xlim=(0., 1.), logx=False, logy=False, **kwargs):
        '''Create a frame (dummy histogram).

        Arguments:
            Title: graph title (D: '')
            href: reference histogram to start from (D: None)
            xlim: x limits, ignored when using a reference histogram (D: 0., 1.)
            logx: log scale in x? (D: False)
            logy: log scale in y? (D: False)
        kwargs:
            see dummy_frame.config
        '''
        if href is not None:
            self.graph = href.Clone()
            self.graph.Reset()
        else:
            self.graph = ROOT.TH1F('dummy_%i' % dummy_frame.dummy_cnt, title, 1000, *xlim)
            dummy_frame.dummy_cnt += 1
        self.graph.SetLineColor(ROOT.kWhite)
        self.graph.GetXaxis().CenterTitle()
        self.graph.GetYaxis().CenterTitle()
        self.logx = logx
        self.logy = logy
        self.config(**kwargs)

    def config(self, **kwargs):
        '''Configure the underlying histogram.

        See dummy_frame.options for the available options.
        '''
        for val in kwargs:
            if val in dummy_frame.options:
                dummy_frame.options[val](self, kwargs[val])
            else:
                print "%s: unknown option: %s" % (__name__, val)
    def Draw(self):
        ROOT.gPad.SetLogx(self.logx)
        ROOT.gPad.SetLogy(self.logy)
        self.graph.Draw()

    def set_logx(self, val):
        self.logx = val
    def set_logy(self, val):
        self.logy = val
