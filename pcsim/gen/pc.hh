#ifndef PCSIM_GEN_PC_LOADED
#define PCSIM_GEN_PC_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/gen/jpsi.hh>
#include <pcsim/physics/bremsstrahlung.hh>
#include <pcsim/physics/pc.hh>

namespace pcsim {
namespace gen {

class pc : public generator<pc, jpsi_event> {
public:
  using base_type = generator<pc, jpsi_event>;

  pc(const configuration& conf, const string_path& path,
       std::shared_ptr<TRandom> r);

  jpsi_event gen_impl();

private:
  // utility functions
  interval<double> calc_s_range() const;

  // cross sections
  const physics::bremsstrahlung brems_;
  const physics::pc_xsec xsec_;
  const double xsec_max_; // the maximum cross section

  // Phase space limits
  const interval<double> s_range_; // s range
};

} // gen
} // pcsim

#endif
