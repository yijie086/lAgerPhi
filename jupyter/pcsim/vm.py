import ROOT
import os
import setup_root

## compile the root macro that contains the actual VM info
ROOT.gROOT.ProcessLine('.L %s.cc+' % os.path.splitext(__file__)[0])

## wrap around the TF1s for more intuitive usage
