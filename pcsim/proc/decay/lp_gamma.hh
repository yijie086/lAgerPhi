#ifndef PCSIM_PROC_DECAY_LP_GAMMA_LOADED
#define PCSIM_PROC_DECAY_LP_GAMMA_LOADED

#include <pcsim/gen/lp_gamma_event.hh>
#include <pcsim/physics/decay.hh>
#include <pcsim/proc/decay/decay.hh>

namespace pcsim {
namespace decay {
// =============================================================================
// DECAY ALL UNSTABLE PARTICLES IN A GAMMA_P EVENT
//
// also handles chained decays
//
// Supported channels:
//  * Pc according to Wang
//  * e+e- decay of VMs
//
// Note: adds event weight to account for branching ratios when not simulating
// the full decay width
// =============================================================================

class lp_gamma : public decay<lp_gamma_event> {
public:
  using base_type = decay<lp_gamma_event>;
  lp_gamma(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}
  virtual void process(lp_gamma_event& e) const;

private:
  void quarkonium_schc(lp_gamma_event& e, const int index) const;
  void pentaquark_wang(lp_gamma_event& e, const int index) const;
};

} // namespace decay
} // namespace pcsim

#endif
