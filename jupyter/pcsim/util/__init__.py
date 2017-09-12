## /usr/bin/python

'''utility routines.'''

## initialize PCSIM/ROOT
import ROOT
if not hasattr(ROOT, "pcsim"):
    from setup_root import setup_root
    setup_root()

## import the sub-modules
from tf1_wrap import *
from classdict import *
from path import *
