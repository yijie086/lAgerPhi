#ifndef PCSIM_PROC_RECONSTRUCTION_LP_GAMMA_LOADED
#define PCSIM_PROC_RECONSTRUCTION_LP_GAMMA_LOADED

#include <pcsim/gen/lp_gamma_event.hh>
#include <pcsim/proc/reconstruction/reconstruction.hh>

namespace pcsim {
namespace reconstruction {
// =============================================================================
// RECONSTRUCTION ALL UNSTABLE PARTICLES IN A GAMMA_P EVENT
//
// also handles chained reconstructions
//
// Supported channels:
//  * Pc according to Wang
//  * e+e- reconstruction of VMs
//
// Note: adds event weight to account for branching ratios when not simulating
// the full reconstruction width
// =============================================================================

class lp_gamma : public reconstruction<lp_gamma_event> {
public:
  using base_type = reconstruction<lp_gamma_event>;
  lp_gamma(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}
  virtual void process(lp_gamma_event& e) const;

private:
};

} // namespace reconstruction
} // namespace pcsim

#endif
