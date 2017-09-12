#include "dataframe.hh"
#include <ROOT/TDataFrame.hxx>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TParticle.h>

namespace dataframe {

class pcsim_vm : public custom_dataframe {
public:
  pcsim_vm(std::string_view fname)
      : custom_dataframe{{"lp_gamma_event", fname}} {
    init();
  }
  virtual ~pcsim_vm() {}

protected:
  virtual col_interface_type custom_defines(TDataFrame& df) {
    // vectors
    auto ndf =
        df.Define("scat", get_vector, {"particles", "scat_index"})
            .Define("beam", get_vector, {"particles", "beam_index"})
            .Define("target", get_vector, {"particles", "target_index"})
            .Define("photon", get_vector, {"particles", "photon_index"})
            .Define("recoil", get_vector, {"particles", "recoil_index"})
            .Define("leading", get_vector, {"particles", "leading_index"})
            .Define("lplus", get_first_child, {"particles", "leading_index"})
            .Define("lminus", get_second_child, {"particles", "leading_index"});
    return ndf;
  }

  static TLorentzVector get_vector(const TClonesArray& particles,
                                   const int16_t index) {
    TLorentzVector v;
    static_cast<TParticle*>(particles.At(index))->Momentum(v);
    return v;
  }
  static TLorentzVector get_child(const TClonesArray& particles, int16_t index,
                                  const int child_index) {
    index =
        static_cast<TParticle*>(particles.At(index))->GetDaughter(child_index);
    return get_vector(particles, index);
  }
  static TLorentzVector get_first_child(const TClonesArray& particles, int16_t index) {
    return get_child(particles, index, 0);
  }
  static TLorentzVector get_second_child(const TClonesArray& particles, int16_t index) {
    return get_child(particles, index, 1);
  }
};
} // namespace dataframe
