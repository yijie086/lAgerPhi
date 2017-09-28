#include "dataframe.hh"
#include <ROOT/TDataFrame.hxx>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <lcio2/MCParticleData.h>
#include <lcio2/ReconstructedParticleData.h>
#include <string>
#include <vector>

#ifdef __CLING__
#pragma link C++ nestedclass;
#pragma link C++ class const lcio2::ReconstructedParticleData+;
#pragma link C++ class std::vector<const lcio2::ReconstructedParticleData*>+;
#pragma link C++ class std::tuple<double, std::vector<const lcio2::ReconstructedParticleData*>, std::vector<const lcio2::ReconstructedParticleData*>>+;
#pragma link C++ class std::vector<std::tuple<double, std::vector<const lcio2::ReconstructedParticleData*>, std::vector<const lcio2::ReconstructedParticleData*>>>+;
#endif

using MCParticles    = std::vector<lcio2::MCParticleData>;
using ReconParticles = std::vector<lcio2::ReconstructedParticleData>;
using ReconParticlePtrs = std::vector<const lcio2::ReconstructedParticleData*>;
using PairRecon = std::tuple<double, std::vector<const lcio2::ReconstructedParticleData*>, std::vector<const lcio2::ReconstructedParticleData*>>;

namespace dataframe {

class anl_vm : public custom_dataframe {
public:
  anl_vm(std::string_view fname,
           const std::vector<std::string>& custom_branches,
           const double scale = 1.)
      : custom_dataframe{{"events", fname}, scale} {
    init(custom_branches);
  }
  anl_vm(std::string_view fname, const double scale = 1.)
      : anl_vm(fname, {}, scale) {}
  virtual ~anl_vm() {}

protected:
  virtual void init(const std::vector<std::string>& custom_branches) {
    def_interface_type df = interface();
    df = df.Define("eta_e", eta_e_truth, {"MCParticle"});
    interface() = df;
  }

private:
  static double eta_e_truth(const MCParticles& parts0) {
    auto ePrime = FourVector(parts0[3].momentum, parts0[3].mass);
    return ePrime.Eta();
  };

  static TLorentzVector FourVector(const std::array<double, 3>& p, double M) {
    double ptot2 = 0;
    for (auto p_x : p)
      ptot2 += (p_x * p_x);
    return TLorentzVector{p[0], p[1], p[2], TMath::Sqrt(ptot2 + M * M)};
  }
};
} // namespace dataframe
