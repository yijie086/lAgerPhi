#ifndef PCSIM_ROOT_DATAFRAME_LP_GAMMA_LOADED
#define PCSIM_ROOT_DATAFRAME_LP_GAMMA_LOADED

#include "dataframe.hh"

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <string>
#include <vector>

namespace pcsim {
namespace root {
namespace dataframe {

class lp_gamma : public custom_dataframe {
public:
  lp_gamma(std::string_view fname, const double scale = 1.)
      : custom_dataframe{{"lp_gamma_event", fname}, scale} {}

protected:
  // utility functions to get certain 4-vectors from the event
  static TLorentzVector get_vector(const TClonesArray& particles,
                                   const int16_t index) {
    TLorentzVector v;
    static_cast<TParticle*>(particles.At(index))->Momentum(v);
    return v;
  }
  template <int ChildIndex>
  static TLorenzVector get_child_vector(const TClonesArray& particles,
                                        const int16_t parent_index) {
    const int16_t index = static_cast<TParticle*>(particles.At(parent_index))
                              ->GetDaughter(ChildIndex);
    return get_vector(particles, index);
  }
};
} // namespace dataframe
} // namespace root
} // namespace pcsim

#endif
