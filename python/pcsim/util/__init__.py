## /usr/bin/python

'''utility routines.'''

## initialize PCSIM/ROOT
import ROOT, setup_root
import setup_root
if not setup_root.root_loaded():
    setup_root.setup_root()

## import the sub-modules
from tf1_wrap import *
from classdict import *
from path import *
