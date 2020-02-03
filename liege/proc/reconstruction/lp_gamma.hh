#ifndef LIEGE_PROC_RECONSTRUCTION_LP_GAMMA_LOADED
#define LIEGE_PROC_RECONSTRUCTION_LP_GAMMA_LOADED

#include <liege/gen/lp_gamma_event.hh>
#include <liege/proc/reconstruction/reconstruction.hh>

namespace liege {
namespace reconstruction {
// =============================================================================
// RECONSTRUCT ADDITIONAL PARTICLES IN A GAMMA_P EVENT
//
// also handles trigger-level cuts
//
// Note: sets the event weight to zero for events that don't pass the cut
//       (will then be removed by the event_generator)
// =============================================================================

class lp_gamma : public reconstruction<lp_gamma_event> {
public:
  using base_type = reconstruction<lp_gamma_event>;
  lp_gamma(const configuration& conf, const string_path& path,
           std::shared_ptr<TRandom> r);
  virtual void process(lp_gamma_event& e) const;

private:
  const bool require_leading_{false};
  const bool veto_leading_{false};
  const bool require_scat_{false};
  const bool veto_scat_{false};
};

} // namespace reconstruction
} // namespace liege

#endif
