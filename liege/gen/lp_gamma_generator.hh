#ifndef LIEGE_GEN_LP_GAMMA_GENERATOR2_LOADED
#define LIEGE_GEN_LP_GAMMA_GENERATOR2_LOADED

#include <memory>
#include <liege/core/generator.hh>
#include <liege/gen/beam/generator.hh>
#include <liege/gen/lp_gamma/generator.hh>
#include <liege/gen/lp_gamma_event.hh>
#include <liege/proc/decay/lp_gamma.hh>
#include <liege/proc/detector/detector.hh>
#include <liege/proc/reconstruction/lp_gamma.hh>

namespace liege {

class lp_gamma_generator
    : public event_generator<lp_gamma_event, lp_gamma_data> {
public:
  using base_type = event_generator<lp_gamma_event, lp_gamma_data>;

  lp_gamma_generator(const configuration& cf, const string_path& path,
                     std::shared_ptr<TRandom> r);

protected:
  virtual lp_gamma_data generate_initial() const;
  virtual void build_event(lp_gamma_event& e) const;

private:
  // initial state generators
  std::shared_ptr<beam::vertex_generator> vertex_gen_;
  std::shared_ptr<beam::primary_generator> lepton_gen_;
  std::shared_ptr<beam::primary_generator> proton_gen_;
  std::shared_ptr<beam::photon_generator> photon_gen_;
  // event processors
  std::shared_ptr<decay::lp_gamma> decay_proc_;
  std::shared_ptr<detector::detector> detector_proc_;
  std::shared_ptr<reconstruction::lp_gamma> rc_proc_;
};

} // namespace liege

#endif
