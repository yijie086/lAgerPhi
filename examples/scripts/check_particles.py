import ROOT
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser('Read lager output')

parser.add_argument('root_file', help='path to the root file.')

args = parser.parse_args()

f = ROOT.TFile.Open(args.root_file)
tree = f.lAger

kcols = ['Q2', 'W', 't']
pcols = ['E', 'Ek', 'Px', 'Py', 'Pz']
pars = { 4: 'e', 5: 'jpsi', 6: 'he4' }

nev = tree.GetEntries()
kbuf = np.zeros(shape=(nev, len(kcols)))
pbufs = {ip: np.zeros(shape=(nev, len(pcols))) for ip in list(pars)}

for iev in np.arange(nev):
    tree.GetEntry(iev)
    kbuf[iev] = (tree.Q2, tree.W, tree.t)
    for ip in list(pars):
        par = tree.particles[ip]
        pbufs[ip][iev] = (par.Energy(), par.Ek(), par.Px(), par.Py(), par.Pz())


res = pd.DataFrame(columns=kcols + [p + '_' + col for _, p in pars.items() for col in pcols ],
                   data=np.concatenate([kbuf] + [pbuf for _, pbuf in pbufs.items()], axis=1))

res.to_csv('particles.csv', index=False)
