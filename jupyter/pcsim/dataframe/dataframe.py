import pcsim.util, ROOT, os

ROOT.gROOT.ProcessLine('#include "%s/dataframe.hh"' % os.path.dirname(__file__))

def make_histo1D(df, href, vname, wname = ''):
    make_histo1D.dummy_cnt += 1
    href.SetName('%s_%i' % (href.GetName(), make_histo1D.dummy_cnt))
    return ROOT.dataframe.make_histo1D(type(df))(df, href, vname, wname)
make_histo1D.dummy_cnt = 0

def make_histo2D(df, href, v1name, v2name, wname = ''):
    make_histo2D.dummy_cnt += 1
    href.SetName('%s_%i' % (href.GetName(), make_histo2D.dummy_cnt))
    return ROOT.dataframe.make_histo2D(type(df))(df, href, v1name, v2name, wname)
make_histo2D.dummy_cnt = 0
