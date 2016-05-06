#ifndef PCSIM_GEN_JPSI_LOADED
#define PCSIM_GEN_JPSI_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/physics/bremsstrahlung.hh>
#include <pcsim/physics/jpsi.hh>

namespace pcsim {
namespace gen {

struct jpsi_event {
  double volume = 0;  // phase space volume used for event generation
  size_t ntrials = 0; // number of trials before this accept was accepted
  double weight = 1;  // event weight (default=1)
  // scattering variables
  double s = 0;
  double W = 0;
  double t = 0;
  // particles
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector recoil;
  TLorentzVector jpsi;
  TLorentzVector positron;
  TLorentzVector electron;
};

class jpsi : public generator<jpsi, jpsi_event> {
public:
  using base_type = generator<jpsi, jpsi_event>;

  jpsi(const configuration& conf, const string_path& path,
       std::shared_ptr<TRandom> r);

  jpsi_event gen_impl();

private:
  // utility function to calculate the maximum cross section and t-window size
  double calc_max_xsec() const;
  interval<double> calc_t_range(double Egamma) const;

  // cross sections
  const physics::bremsstrahlung brems_;
  const physics::jpsi_xsec xsec_;

  // phase space limits
  const interval<double> s_range_;
  const interval<double> t_range_;

  const double xsec_max_; // the maximum cross section
};

} // gen
} // pcsim

#endif
