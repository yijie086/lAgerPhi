#ifndef LIEGE_PROC_DECAY_LP_GAMMA_LOADED
#define LIEGE_PROC_DECAY_LP_GAMMA_LOADED

#include <liege/core/particle.hh>
#include <liege/gen/lp_gamma_event.hh>
#include <liege/physics/decay.hh>
#include <liege/proc/decay/decay.hh>

namespace liege {
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
class radiative_decay_vm {
public:
  radiative_decay_vm();
  void process(lp_gamma_event& e, const int vm_index);
};

class lp_gamma : public decay<lp_gamma_event> {
public:
  using base_type = decay<lp_gamma_event>;
  lp_gamma(const configuration&, const string_path&,
           std::shared_ptr<TRandom> r);
  virtual void process(lp_gamma_event& e) const;

private:
  void quarkonium_schc(lp_gamma_event& e, const int index) const;
  void pentaquark_wang(lp_gamma_event& e, const int index) const;

  const particle vm_decay_lplus_;
  const particle vm_decay_lminus_;
  const double vm_decay_br_;
  std::unique_ptr<radiative_decay_vm> radiative_decay_;
};

} // namespace decay
} // namespace liege

#endif
