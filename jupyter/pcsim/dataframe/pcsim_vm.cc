#include "dataframe.hh"
#include <ROOT/TDataFrame.hxx>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <string>
#include <vector>

namespace dataframe {

class pcsim_vm : public custom_dataframe {
public:
  pcsim_vm(std::string_view fname,
           const std::vector<std::string>& custom_branches,
           const double scale = 1.)
      : custom_dataframe{{"lp_gamma_event", fname}, scale} {
    init(custom_branches);
  }
  pcsim_vm(std::string_view fname, const double scale = 1.)
      : pcsim_vm(fname, {}, scale) {}
  virtual ~pcsim_vm() {}

protected:
  virtual void init(const std::vector<std::string>& custom_branches) {
    def_interface_type df = interface();
    for (const auto& branch : custom_branches) {
      if (branch == "scat" || branch == "beam" || branch == "target" ||
          branch == "photon" || branch == "recoil" || branch == "leading") {
        df = define_vector(df, branch);
      } else if (branch == "lplus") {
        df = define_child_vector(df, branch, "leading", 0);
      } else if (branch == "lminus") {
        df = define_child_vector(df, branch, "leading", 1);
      } else if (branch == "eta_scat" || branch == "eta_recoil") {
        df = define_eta(df, branch);
      } else if (branch == "eta_lplus") {
        df = define_child_eta(df, branch, "leading", 0);
      } else if (branch == "eta_lminus") {
        df = define_child_eta(df, branch, "leading", 1);
      } else if (branch == "theta_scat" || branch == "theta_recoil") {
        df = define_theta(df, branch);
      } else if (branch == "theta_lplus") {
        df = define_child_theta(df, branch, "leading", 0);
      } else if (branch == "theta_lminus") {
        df = define_child_theta(df, branch, "leading", 1);
      } else {
        std::cerr << "UNKNOWN BRANCH REQUEST: " << branch << std::endl;
      }
    }
    interface() = df;
  }

private:
  static def_interface_type define_vector(def_interface_type& df,
                                          const std::string& branch) {
    return df.Define(branch, get_vector, {"particles", branch + "_index"});
  }
  static def_interface_type define_eta(def_interface_type& df, const std::string& branch) {
    const std::string particle = branch.substr(4, branch.size());
    return df.Define(branch, get_eta, {"particles", particle + "_index"});
  }
  static def_interface_type define_theta(def_interface_type& df, const std::string& branch) {
    const std::string particle = branch.substr(6, branch.size());
    return df.Define(branch, get_theta, {"particles", particle + "_index"});
  }

  static def_interface_type define_child_vector(def_interface_type& df,
                                                const std::string& branch,
                                                const std::string& parent,
                                                const int child_index) {
    return df.Define(
        branch,
        [=](const TClonesArray& particles, const int16_t parent_index) {
          return get_child_vector(particles, parent_index, child_index);
        },
        {"particles", parent + "_index"});
  }
  static def_interface_type define_child_eta(def_interface_type& df,
                                             const std::string& branch,
                                             const std::string& parent,
                                             const int child_index) {
    return df.Define(
        branch,
        [=](const TClonesArray& particles, const int16_t parent_index) {
          return get_child_eta(particles, parent_index, child_index);
        },
        {"particles", parent + "_index"});
  }
  static def_interface_type define_child_theta(def_interface_type& df,
                                               const std::string& branch,
                                               const std::string& parent,
                                               const int child_index) {
    return df.Define(
        branch,
        [=](const TClonesArray& particles, const int16_t parent_index) {
          return get_child_theta(particles, parent_index, child_index);
        },
        {"particles", parent + "_index"});
  }

  static TLorentzVector get_vector(const TClonesArray& particles,
                                   const int16_t index) {
    TLorentzVector v;
    static_cast<TParticle*>(particles.At(index))->Momentum(v);
    return v;
  }
  static double get_eta(const TClonesArray& particles, const int16_t index) {
    return static_cast<TParticle*>(particles.At(index))->Eta();
  }
  static double get_theta(const TClonesArray& particles, const int16_t index) {
    return static_cast<TParticle*>(particles.At(index))->Theta();
  }
  static TLorentzVector get_child_vector(const TClonesArray& particles,
                                         int16_t index, const int child_index) {
    index =
        static_cast<TParticle*>(particles.At(index))->GetDaughter(child_index);
    return get_vector(particles, index);
  }
  static double get_child_eta(const TClonesArray& particles,
                                      int16_t index, const int child_index) {
    index =
        static_cast<TParticle*>(particles.At(index))->GetDaughter(child_index);
    return get_eta(particles, index);
  }
  static double get_child_theta(const TClonesArray& particles, int16_t index,
                                const int child_index) {
    index =
        static_cast<TParticle*>(particles.At(index))->GetDaughter(child_index);
    return get_theta(particles, index);
  }
};
} // namespace dataframe
