#ifndef OLEKSII_BH_IMPL_LOADED
#define OLEKSII_BH_IMPL_LOADED

namespace oleksii_bh_impl {

// calculate the cross section
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame
double calc_xsec(const double t, const double q2, const double Egamma,
                 const double thetaCM, const double phiCM);
} // namespace oleksii_bh_impl

#endif
