#ifndef PCSIM_GEN_BEAM_VERTEX_GEN_LOADED
#define PCSIM_GEN_BEAM_VERTEX_GEN_LOADED

#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam/generator.hh>

namespace pcsim {
namespace beam {

// =============================================================================
// beam::beam
//
// constant 4-vector beam
// =============================================================================
class origin_vertex : public vertex_generator {
public:
  origin_vertex(const configuration& cf, const string_path& path,
                std::shared_ptr<TRandom> r)
      : vertex_generator{std::move(r)}, vertex_{0, 0, 0, 0} {
    LOG_INFO("beam::origin_vertex", "Vertex generator initialized");
  }

  virtual particle::XYZTVector generate() { return vertex_; }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  const particle::XYZTVector vertex_;
};

class linear_vertex : public vertex_generator {
public:
  linear_vertex(const configuration& cf, const string_path& path,
                std::shared_ptr<TRandom> r)
      : vertex_generator{std::move(r)}
      , range_{cf.get_range<double>(path / "range")} {
    LOG_INFO("beam::linear_vertex", "Vertex generator initialized");
    LOG_INFO("beam::linear_vertex", "Vertex range [cm]: [" +
                                        std::to_string(range_.min) + ", " +
                                        std::to_string(range_.max) + "]");
  }

  virtual particle::XYZTVector generate() {
    const double vz = rng()->Uniform(range_.min, range_.max);
    LOG_JUNK2("linear_vertex", "Vertex position [cm]: " + std::to_string(vz));
    return {0., 0., vz, 0};
  }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  const interval<double> range_;
};

} // namespace beam
} // namespace pcsim

#endif
