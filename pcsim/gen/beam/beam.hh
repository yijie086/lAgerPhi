#ifndef PCSIM_BEAM_BEAM_LOADED
#define PCSIM_BEAM_BEAM_LOADED

#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>

namespace pcsim {
namespace beam {

// =============================================================================
// beam::data
//
// generic beam data
// =============================================================================
class data : public generator_data {
public:
  data() = default;
  data(const particle& part) : beam_{part} {}
  data(const particle& part, const double xs)
      : generator_data{xs}, beam_{part} {}
  data(const double xs) : generator_data{xs} {}

  const particle& beam() const { return beam_; }

private:
  particle beam_;
};

// =============================================================================
// beam::primary
//
// constant 4-vector beam
// =============================================================================
class primary : public generator<particle> {
public:
  static factory<primary> factory;

  primary(const configuration& cf, const string_path& path,
       std::shared_ptr<TRandom> r)
      : generator{std::move(r)}
      , beam_{static_cast<pdg_id>(cf.get<int>(path / "type")),
              cf.get_vector3<particle::XYZVector>(path / "dir"),
              cf.get<double>(path / "energy")} {
    LOG_INFO("beam::primary", "type: " + std::string(beam_.pdg->GetName()));
    LOG_INFO("beam::primary", "energy [GeV]: " + std::to_string(beam_.mom.E()));
  }

  virtual particle generate() { return beam_; }

protected:
  const particle beam_;
};

} // namespace beam
} // namespace pcsim

#endif
