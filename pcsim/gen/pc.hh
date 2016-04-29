#ifndef PCSIM_GEN_PC_LOADED
#define PCSIM_GEN_PC_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/gen/bremsstrahlung.hh>
#include <pcsim/gen/jpsi.hh>
#include <pcsim/physics/pc.hh>

namespace pcsim {
namespace gen {

class pc : public generator<pc, jpsi_event> {
public:
  using base_type = generator<pc, jpsi_event>;

  pc(const ptree& settings, const string_path& path,
       std::shared_ptr<TRandom> r);

  jpsi_event gen_impl(const photon_beam& photon);

private:
  const physics::pc_xsec xsec_;
  const double me_;  // electron mass
  const double Mjp_; // J/Psi pole mass
  const double Mp_;  // proton mass
  const double Wjp_; // J/Psi decay width
  const double Bje_; // J/Psi --> e+e- branching ratio
  const double ctheta_min_; // min and max cos(theta) that can be reached (equal
  const double ctheta_max_; // to -1 and +1)
};

} // gen
} // pcsim

#endif
