#ifndef PCSIM_GEN_BEAM_LOADED
#define PCSIM_GEN_BEAM_LOADED

#include <TLorentzVector.h>
#include <TVector3.h>
#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>

namespace pcsim {
namespace gen {

// =============================================================================
// simple beam along a constant four-vector
// =============================================================================
class beam : public generator<particle> {
public:
  static factory<beam> factory;

  beam(const configuration& cf, const string_path& path,
       std::shared_ptr<TRandom> r)
      : generator{cf, path, std::move(r)}
      , beam_{static_cast<pdg_id>(conf().get<int>("type")),
              conf().get_vector3<TVector3>("dir"),
              conf().get<double>("energy")} {}

  virtual particle generate() { return beam_; }

private:
  const particle beam_;
};

}
}

#endif
