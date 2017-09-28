import pcsim.util, ROOT

ROOT.gROOT.ProcessLine('#include <pcsim/root/physics/physics.hh>')

M_EL = ROOT.pcsim.root.physics.M_EL
M_P = ROOT.pcsim.root.physics.M_P
M_JPSI = ROOT.pcsim.root.physics.M_JPSI
M_UPSILON = ROOT.pcsim.root.physics.M_UPSILON

THRESHOLD_JPSI = M_P + M_JPSI;
THRESHOLD_UPSILON = M_P + M_UPSILON;
