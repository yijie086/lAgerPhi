import pcsim.util, ROOT, os

## compile the associated root macro that contains the actual VM info if needed
ROOT.gROOT.ProcessLine('.L %s.cc+' % os.path.splitext(__file__)[0])

from ROOT.dataframe import anl_vm
