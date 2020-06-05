import ROOT
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser('Read lager output')

parser.add_argument('root_file', help='path to the root file.')

args = parser.parse_args()

f = ROOT.TFile.Open(args.root_file)
tree = f.lAger
tree.GetEntry(0)

for par in tree.particles:
    print('Name = {} Status = {:d} Energy = {:.2f} GeV Ek = {:.4f} GeV'.format(par.GetName(), par.GetStatusCode(), par.Energy(), par.Ek()))
