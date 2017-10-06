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
  // data structures
  using particle_type = TParticle;
  using particles_type = TClonesArray;
  using particle_crefs_type =
      custom_dataframe::particle_crefs_type<particle_type>;
  // mc_particle and rc_particle are stored in the exact same structure
  using mc_particle_type = particle_type;
  using mc_particles_type = particles_type;
  using mc_particle_crefs_type = particle_crefs_type;
  using rc_particle_type = particle_type;
  using rc_particles_type = particles_type;
  using rc_particle_crefs_type = particle_crefs_type;

  // column identifiers
  constexpr static const char* MC_PARTS = "particles";
  constexpr static const char* RC_PARTS = "rc_particles";

  // tree name
  constexpr static const char* TREE_NAME = "lp_gamma_event";

  // constructors
  lp_gamma(std::string_view fname, const double scale = 1.)
      : custom_dataframe{{TREE_NAME, fname}, scale} {}
  lp_gamma(const std::vector<std::string>& fnames, const double scale = 1.)
      : custom_dataframe{{TREE_NAME, fnames}, scale} {}

  // utility functions to extract data from the event
  static ROOT::Math::XYZTVector vector(const particles_type& parts,
                                       const int16_t index) {
    return vector(*static_cast<particle_type*>(parts.At(index)));
  }
  template <int16_t index>
  ROOT::Math::XYZTVector vector(const particles_type& parts) {
    return vector(parts, index);
  }
  static ROOT::Math::XYZTVector vector(const particle_type& part) {
    TLorentzVector v;
    part.Momentum(v);
    return {v.X(), v.Y(), v.Z(), v.T()};
  }
  template <int ChildIndex>
  static ROOT::Math::XYZTVector child_vector(const particles_type& parts,
                                             const int16_t parent_index) {
    const int16_t index = static_cast<particle_type*>(parts.At(parent_index))
                              ->GetDaughter(ChildIndex);
    return vector(parts, index);
  }
  template <int pid> static size_t pid_count(const particles_type& parts) {
    size_t cnt = 0;
    const size_t N = parts.GetSize();
    for (size_t i = 0; i < N; ++i) {
      const auto& part = *static_cast<TParticle*>(parts.At(i));
      if (part.GetWeight() == 1. && pid == part.GetPdgCode()) {
        cnt += 1;
      }
    }
    return cnt;
  }
  template <int pid>
  static particle_crefs_type select_pid(const particles_type& parts) {
    particle_crefs_type selected;
    const size_t N = parts.GetSize();
    for (size_t i = 0; i < N; ++i) {
      const auto& part = *static_cast<TParticle*>(parts.At(i));
      if (part.GetWeight() == 1. && pid == part.GetPdgCode()) {
        selected.push_back(part);
      }
    }
    return selected;
  }
};
} // namespace dataframe
} // namespace root
} // namespace pcsim

#endif
