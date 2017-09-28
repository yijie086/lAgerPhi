#ifndef PCSIM_ROOT_DATAFRAME_ANL_LOADED
#define PCSIM_ROOT_DATAFRAME_ANL_LOADED

#include <pcsim/root/config.hh>

#ifdef SUPPORT_ANL

#include "dataframe.hh"
#include <Math/Vector4D.h>
#include <algorithm>
#include <functional>
#include <lcio2/MCParticleData.h>
#include <lcio2/ReconstructedParticleData.h>
#include <vector>

namespace pcsim {
namespace root {
namespace dataframe {

class anl : public custom_dataframe {
public:
  using particles_type = std::vector<lcio2::ReconstructedParticleData>;
  using particle_crefs_type = std::vector<
      std::reference_wrapper<const lcio2::ReconstructedParticleData>>;
  using mc_particles_type = std::vector<lcio2::MCParticleData>;
  using mc_particle_crefs_type =
      std::vector<std::reference_wrapper<const lcio2::MCParticleData>>;
  using vector_type = ROOT::Math::PxPyPzMVector;

  anl(const std::string_view& fname, const double scale = 1.)
      : custom_dataframe({"events", fname}, scale) {}
  anl(const std::vector<std::string>& fnames, const double scale = 1.)
      : custom_dataframe({"events", fnames}, scale) {}

  template <class ParticleVector>
  static vector_type vector(const ParticleVector& parts, const int index) {
    const auto v3 = parts[index].momentum;
    const auto m = parts[index].mass;
    return {v3[0], v3[1], v3[2], m};
  }
  static size_t pid_count(const particles_type& parts, const int pid) {
    return std::count_if(parts.begin(), parts.end(),
                         [=](const auto& part) { return part.type == pid; });
  }
  static particle_crefs_type select_pid(const particles_type& parts,
                                        const int pid) {
    particle_crefs_type selected;
    for (const auto& part : parts) {
      if (part.type == pid) {
        selected.push_back(part);
      }
    }
    return selected;
  }
};

} // namespace dataframe
} // namespace root
} // namespace pcsim

#endif // SUPPORT_ANL

#endif
