#ifndef PCSIM_GEN_BEAM_PRIMARY_GEN_LOADED
#define PCSIM_GEN_BEAM_PRIMARY_GEN_LOADED

#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam/beam.hh>
#include <pcsim/gen/beam/primary.hh>

namespace pcsim {
namespace beam {

// =============================================================================
// beam::beam
//
// constant 4-vector beam
// =============================================================================
class beam : public generator {
public:
  primary_gen(const configuration& cf, const string_path& path,
              std::shared_ptr<TRandom> r)
      : generator{std::move(r)}
      , beam_{static_cast<pdg_id>(cf.get<int>(path / "particle_type")),
              cf.get_vector3<particle::XYZVector>(path / "dir"),
              cf.get<double>(path / "energy"), particle::status_code::BEAM} {
    LOG_INFO("beam::primary_gen", "type: " + std::string(beam_.pdg()->GetName()));
    LOG_INFO("beam::primary_gen",
             "energy [GeV]: " + std::to_string(beam_.energy()));
  }

  virtual primary generate() { return {beam_}; }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  const particle beam_;
};

} // namespace beam
} // namespace pcsim

#endif
