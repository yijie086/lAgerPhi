#ifndef PCSIM_GEN_EXCLUSIVE_LOADED
#define PCSIM_GEN_EXCLUSIVE_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <memory>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam.hh>
#include <pcsim/gen/photon.hh>

namespace pcsim {
namespace gen {

// describes photon + target ==> a + b events
// vm is the primary produced particle, and recoil is the recoil target
// (possible excited)
struct exclusive_data : generator_data {
  double t;
  particle vm;     // vector meson
  particle recoil; // recoil
  exclusive_data(double xs, double ps) : generator_data(xs, ps) {}
  exclusive_data() = default;
};

// base class for exclusive events
class exclusive : public generator<exclusive_data, photon_data, particle> {
public:
  exclusive(const configuration& conf, const string_path& path,
            std::shared_ptr<TRandom> r)
      : generator{conf, path, std::move(r)} {}
  exclusive_data generate(const photon_data& photon,
                                  const particle& target) = 0;
};

class brodsky_tchannel : public exclusive {
public:
  brodsky_tchannel(const configuration& conf, const string_path& path,
                   std::shared_ptr<TRandom> r);

  virtual double max_cross_section() const { return max_; }
  virtual exclusive_data generate(const photon_data& photon,
                                  const particle& target);

protected:
  double dsigma_dt(const double W2, const double t, const double Mt) const;
  double dsigma_dexp_bt(const double W2, const double Mt) const;

private:
  double calc_max(const configuration& cf) const;
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf) const;

  // calc t range
  interval<double> t_range(const double W2, const double Q2,
                           const double Mt) const;
  interval<double> exp_bt_range(const double W2, const double Q2,
                                const double Mt) const;

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;
  const double threshold2_;

  // cross section settings
  const double b_;   // b-parameters
  const double c2g_; // constant for 2-gluon amplitude

  // t range and cross section maxima
  // note: we actually throw flat in exp(bt), not in t
  const interval<double> max_t_range_;
  const interval<double> max_exp_bt_range_;
  const double max_;
};





} // gen
} // pcsim

#endif
