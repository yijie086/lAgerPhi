#ifndef OLEKSII_JPSI_IMPL_LOADED
#define OLEKSII_JPSI_IMPL_LOADED

namespace oleksii_jpsi_impl {

std::complex<double> fT1(const double q2, const double Egamma);

double TT(const double t, const double q2, const double Egamma,
          const double thetaCM, const double phiCM);

double Re_jpsi_T(const double q2, const double Egamma, const double T_0 = 0);
double Im_jpsi_T(const double q2, const double Egamma, const double T_0 = 0);
double Re_jpsi_fT(const double q2, const double Egamma, const double t,
                  const double T_0 = 0);
double Im_jpsi_fT(const double q2, const double Egamma, const double t,
                  const double T_0 = 0);

// calculate the cross section
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame
double calc_xsec(const double t, const double q2, const double Egamma,
                 const double thetaCM, const double phiCM,
                 const double T_0 = 0);
} // namespace oleksii_jpsi_impl

#endif
