#ifndef PCSIM_GEN_JPSI_LOADED
#define PCSIM_GEN_JPSI_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/gen/bremsstrahlung.hh>
#include <pcsim/physics/jpsi.hh>

namespace pcsim {
namespace gen {

struct jpsi_event {
  // cross section and integrated photon flux
  double xsec = 0;
  double flux = 0;
  // scattering variables
  double s = 0;
  double t = 0;
  double tmin = 0;
  double tmax = 0;
  double W = 0;
  // particles
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector recoil;
  TLorentzVector jpsi;
  TLorentzVector positron;
  TLorentzVector electron;
  // is this a good event?
  bool good = false;
};

class jpsi : generator<jpsi, jpsi_event> {
public:
  using base_type = generator<jpsi, jpsi_event>;

  jpsi(const ptree& settings, const string_path& path,
       std::shared_ptr<TRandom> r);

  jpsi_event gen_impl(const photon_beam& photon);

private:
  const physics::jpsi_xsec xsec_;
  const double me_;  // electron mass
  const double Mjp_; // J/Psi pole mass
  const double Mp_;  // proton mass
  const double Wjp_; // J/Psi decay width
  const double ctheta_min_; // min and max cos(theta) that can be reached (equal
  const double ctheta_max_; // to -1 and +1)
};

} // gen
} // pcsim

#endif
