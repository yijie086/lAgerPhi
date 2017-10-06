#ifndef PCSIM_ROOT_DATAFRAME_ANL_LOADED
#define PCSIM_ROOT_DATAFRAME_ANL_LOADED

#include <pcsim/root/config.hh>

#ifdef SUPPORT_ANL

#include "dataframe.hh"
#include <Math/Vector4D.h>
#include <TRandom3.h>
#include <algorithm>
#include <functional>
#include <lcio2/MCParticleData.h>
#include <lcio2/ReconstructedParticleData.h>
#include <pcsim/core/particle.hh>
#include <vector>

namespace pcsim {
namespace root {
namespace dataframe {

class anl : public custom_dataframe {
public:
  // data structures
  using mc_particle_type = lcio2::MCParticleData;
  using mc_particles_type = std::vector<mc_particle_type>;
  using mc_particle_crefs_type =
      custom_dataframe::particle_crefs_type<mc_particle_type>;
  using rc_particle_type = lcio2::ReconstructedParticleData;
  using rc_particles_type = std::vector<rc_particle_type>;
  using rc_particle_crefs_type =
      custom_dataframe::particle_crefs_type<rc_particle_type>;

  // column identifiers
  constexpr static const char* MC_PARTS = "MCParticle";
  constexpr static const char* RC_PARTS = "PandoraPFOCollection";

  // tree name
  constexpr static const char* TREE_NAME = "events";

  // constructors
  anl(const std::string_view& fname, const double scale = 1.)
      : custom_dataframe({TREE_NAME, fname}, scale) {}
  anl(const std::vector<std::string>& fnames, const double scale = 1.)
      : custom_dataframe({TREE_NAME, fnames}, scale) {}

  // utility functions to extract data from the event
  template <size_t index>
  static ROOT::Math::PxPyPzMVector vector(const rc_particles_type& parts) {
    return vector(parts[index]);
  }
  static ROOT::Math::PxPyPzMVector vector(const rc_particle_type& part) {
    const auto& p = part.momentum;
    const auto mass = part.mass;
    return {p[0], p[1], p[2], mass};
  }
  template <size_t index>
  static ROOT::Math::XYZTVector vector(const mc_particles_type& parts) {
    return vector(parts[index]);
  }
  static ROOT::Math::XYZTVector vector(const mc_particle_type& part) {
    const auto pid = static_cast<pcsim::pdg_id>(part.pdg);
    const auto& p = part.momentum;
    // this will automatically initiate the particle with its pole mass (be
    // careful with unstable particles and virtual particles!
    const particle tmp{pid, {p[0], p[1], p[2]}};
    return tmp.p();
  }

  template <int pid> static size_t pid_count(const rc_particles_type& parts) {
    return std::count_if(
        parts.begin(), parts.end(),
        [](const rc_particle_type& part) { return part.type == pid; });
  }
  template <int pid>
  static rc_particle_crefs_type select_pid(const rc_particles_type& parts) {
    rc_particle_crefs_type selected;
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
