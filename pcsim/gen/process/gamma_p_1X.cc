#include "gamma_p_1X.hh"
#include <TF1.h>
#include <TMath.h>
#include <pcsim/core/logger.hh>
#include <pcsim/physics/photon.hh>
#include <pcsim/physics/vm.hh>

namespace pcsim {
namespace process {

factory<gamma_p_1X, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    gamma_p_1X::factory;
FACTORY_REGISTER(gamma_p_1X, gamma_p_1X_qpq, "gamma_p_1X_qpq");

// =============================================================================
// Constructor for gamma_p_1X data, calculates the final state four-vector in
// the lab-frame
// =============================================================================
gamma_p_1X_data::gamma_p_1X_data(const beam::photon_data& photon,
                                 const beam::data& target, const particle& X1,
                                 const double xs, const double R)
    : generator_data{xs}
    , R_{R}
    , X_{X1.type(), (photon.beam().p() + target.beam().p())} {}

// =============================================================================
// Constructor for process::gamma_p_1X_qpq
// =============================================================================
gamma_p_1X_qpq::gamma_p_1X_qpq(const configuration& cf, const string_path& path,
                               std::shared_ptr<TRandom> r)
    : gamma_p_1X{r}
    , vm_pole_{static_cast<pdg_id>(cf.get<int>(path / "vm_type"))}
    , qpq_{static_cast<pdg_id>(cf.get<int>(path / "qpq_type"))}
    , W2_range_{calc_W2_range(cf.get<double>(path / "n_sigma"))}
    , qpq_amplitude_{cf.get<double>(path / "qpq_amplitude")}
    , qpq_coupling_{cf.get<double>(path / "qpq_couping")}
    , R_vm_c_{cf.get<double>(path / "R_vm_c")}
    , R_vm_n_{cf.get<double>(path / "R_vm_n")}
    , dipole_n_{cf.get<double>(path / "dipole_n")}
    , max_{calc_max_xsec(cf)} {
  LOG_INFO("gamma_p_1X_qpq", "amplitude: " + std::to_string(qpq_amplitude_));
  LOG_INFO("gamma_p_1X_qpq", "coupling: " + std::to_string(qpq_coupling_));
  LOG_INFO("gamma_p_1X_qpq",
           "R_vm c-parameter: " + std::to_string(R_vm_c_));
  LOG_INFO("gamma_p_1X_qpq",
           "R_vm n-parameter (power): " + std::to_string(R_vm_n_));
  LOG_INFO("gamma_p_1X_qpq",
           "'Dipole' FF power: " + std::to_string(dipole_n_));
  LOG_INFO("gamma_p_1X_qpq",
           "Q-Pentaquark: " + std::string(qpq_.pdg()->GetName()));
  LOG_INFO("gamma_p_1X_qpq",
           "VM Pole: " + std::string(vm_pole_.pdg()->GetName()));
}

gamma_p_1X_data gamma_p_1X_qpq::generate(const beam::photon_data& photon,
                                         const beam::data& target) {

  // check if we are in the correct W2 range
  if (W2_range_.excludes(photon.W2())) {
    LOG_JUNK("gamma_p_1X_qpq",
             "Event outside of W2 range - W2: " + std::to_string(photon.W2()) +
                 " outside of [" + std::to_string(W2_range_.min) + ", " +
                 std::to_string(W2_range_.max) + "]");
    return {0.};
  }

  // no generation step necessary, we just have to evaluate the cross section
  // and create the Q-Pq

  // evaluate the cross section
  const double xs_R = R(photon.Q2());
  const double xs_dipole = dipole(photon.Q2());
  const double xs_photo = sigma(photon.W2());
  const double xs = (1 + photon.epsilon() * xs_R) * xs_dipole * xs_photo;

  LOG_JUNK("gamma_p_1X_qpq",
           "xsec: " + std::to_string(xs_photo) + " < " + std::to_string(max_));
  LOG_JUNK("gamma_p_1X_qpq", "R: " + std::to_string(xs_R));
  LOG_JUNK("gamma_p_1X_qpq", "dipole: " + std::to_string(xs_dipole));

  // return a new 1X event
  return {photon, target, qpq_, xs, xs_R};
}


// =============================================================================
// gamma_p_1X_qpq::calc_max_xsec(cf)
//
// Utility function for the generator initialization
//
// max cross section as defined by the beam/target and phase-space settings
//
// Scan for cross section maximum between the W2 limits
// =============================================================================
double gamma_p_1X_qpq::calc_max_xsec(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the lepton
  // beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/dir"),
                        cf.get<double>("beam/energy")};
  const particle target{static_cast<pdg_id>(cf.get<int>("target.type")),
                        cf.get_vector3<particle::XYZVector>("target/dir"),
                        cf.get<double>("target/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W
  const double W2max = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                          (photon.p() + target.p()).M2())
                                   : (photon.p() + target.p()).M2();
  // get the minimum W
  const double W2min =
      opt_W_range ? fmin(opt_W_range->min * opt_W_range->min, target.mass2())
                  : target.mass2();
  auto func = [this](double* W2, double*) { return this->sigma(W2[0]); };
  TF1 f("qpq_xsec", func, W2min, W2max, 0.);
  return f.GetMaximum() * 1.0000001;
}

// =============================================================================
// gamma_p_1X_qpq::calc_mass2_range()
//
// Utility function to calculate the W2 range we allow for q-Pq production
// =============================================================================
interval<double> gamma_p_1X_qpq::calc_W2_range(const double n_sigma) const {
  double min = qpq_.pole_mass() - n_sigma * qpq_.width();
  double max = qpq_.pole_mass() - n_sigma * qpq_.width();
  return {min * min, max * max};
}

// =============================================================================
// gamma_p_1X_qpq::sigma()
// gamma_p_1X_qpq::R()
// gamma_p_1X_qpq::dipole()
//
// Utility functions to calculate the cross section components
// =============================================================================
double gamma_p_1X_qpq::sigma(const double W2) const {
  return qpq_coupling_ * qpq_coupling_ * qpq_amplitude_ *
         TMath::Gaus(sqrt(W2), qpq_.mass(), qpq_.width());
}
double gamma_p_1X_qpq::R(const double Q2) const {
  return physics::R_vm_martynov(Q2, vm_pole_.mass(), R_vm_c_, R_vm_n_);
}
double gamma_p_1X_qpq::dipole(const double Q2) const {
  return physics::dipole_ff_vm(Q2, vm_pole_.mass(), dipole_n_);
}

} // namespace process
} // namespace pcsim
