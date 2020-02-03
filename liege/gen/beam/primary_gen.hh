#ifndef LIEGE_GEN_BEAM_PRIMARY_GEN_LOADED
#define LIEGE_GEN_BEAM_PRIMARY_GEN_LOADED

#include <liege/core/factory.hh>
#include <liege/core/generator.hh>
#include <liege/core/particle.hh>
#include <liege/gen/beam/generator.hh>
#include <liege/gen/beam/primary.hh>

namespace liege {
namespace beam {

// =============================================================================
// beam::beam
//
// constant 4-vector beam
// =============================================================================
class beam : public primary_generator {
public:
  beam(const configuration& cf, const string_path& path,
       std::shared_ptr<TRandom> r)
      : primary_generator{std::move(r)}
      , beam_{static_cast<pdg_id>(cf.get<int>(path / "particle_type")),
              cf.get_vector3<particle::XYZVector>(path / "dir"),
              cf.get<double>(path / "energy"), particle::status_code::BEAM} {
    LOG_INFO("beam::primary_gen",
             "type: " + std::string(beam_.pdg()->GetName()));
    LOG_INFO("beam::primary_gen",
             "energy [GeV]: " + std::to_string(beam_.energy()));
  }

  virtual primary generate(const particle::XYZTVector& vertex) {
    particle ret = beam_;
    ret.vertex() = vertex;
    return {ret};
  }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  const particle beam_;
};

} // namespace beam
} // namespace liege

#endif
