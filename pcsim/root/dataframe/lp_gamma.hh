#ifndef PCSIM_ROOT_DATAFRAME_LP_GAMMA_LOADED
#define PCSIM_ROOT_DATAFRAME_LP_GAMMA_LOADED

#include "dataframe.hh"

#include <Math/Vector4D.h>
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
  using vector_type = ROOT::Math::XYZTVector;

  lp_gamma(std::string_view fname, const double scale = 1.)
      : custom_dataframe{{"lp_gamma_event", fname}, scale} {}
  lp_gamma(const std::vector<std::string>& fnames, const double scale = 1.)
      : custom_dataframe{{"lp_gamma_event", fnames}, scale} {}

  // utility functions to get certain 4-vectors from the event
  static vector_type vector(const TClonesArray& particles,
                            const int16_t index) {
    TLorentzVector v;
    static_cast<TParticle*>(particles.At(index))->Momentum(v);
    return {v.X(), v.Y(), v.Z(), v.T()};
  }
  template <int ChildIndex>
  static vector_type child_vector(const TClonesArray& particles,
                                  const int16_t parent_index) {
    const int16_t index = static_cast<TParticle*>(particles.At(parent_index))
                              ->GetDaughter(ChildIndex);
    return vector(particles, index);
  }
};
} // namespace dataframe
} // namespace root
} // namespace pcsim

#endif
