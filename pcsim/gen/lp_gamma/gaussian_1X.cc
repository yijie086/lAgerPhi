#include "gaussian_1X.hh"
#include <TF1.h>
#include <TMath.h>
#include <pcsim/core/logger.hh>
#include <pcsim/physics/photon.hh>
#include <pcsim/physics/vm.hh>

namespace pcsim {
namespace lp_gamma {

// =============================================================================
// Constructor for lp_gamma::gaussian_qpq
// =============================================================================
gaussian_qpq::gaussian_qpq(const configuration& cf, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : base_type{r}
    , vm_pole_{static_cast<pdg_id>(cf.get<int>(path / "vm_type"))}
    , qpq_{static_cast<pdg_id>(cf.get<int>(path / "qpq_type")),
           particle::status_code::UNSTABLE}
    , W2_range_{calc_W2_range(cf.get<double>(path / "n_sigma"))}
    , qpq_amplitude_{cf.get<double>(path / "qpq_amplitude")}
    , qpq_coupling_{cf.get<double>(path / "qpq_coupling")}
    , R_vm_c_{cf.get<double>(path / "R_vm_c")}
    , R_vm_n_{cf.get<double>(path / "R_vm_n")}
    , dipole_n_{cf.get<double>(path / "dipole_n")}
    , max_{calc_max_xsec(cf)} {
  LOG_INFO("gaussian_qpq", "amplitude: " + std::to_string(qpq_amplitude_));
  LOG_INFO("gaussian_qpq", "coupling: " + std::to_string(qpq_coupling_));
  LOG_INFO("gaussian_qpq", "R_vm c-parameter: " + std::to_string(R_vm_c_));
  LOG_INFO("gaussian_qpq",
           "R_vm n-parameter (power): " + std::to_string(R_vm_n_));
  LOG_INFO("gaussian_qpq", "'Dipole' FF power: " + std::to_string(dipole_n_));
  LOG_INFO("gaussian_qpq",
           "Q-Pentaquark: " + std::string(qpq_.pdg()->GetName()));
  LOG_INFO("gaussian_qpq",
           "VM Pole: " + std::string(vm_pole_.pdg()->GetName()));
}

lp_gamma_event gaussian_qpq::generate(const lp_gamma_data& initial) {
  const auto& gamma = initial.photon();

  // check if we are in the correct W2 range
  if (W2_range_.excludes(gamma.W2())) {
    LOG_JUNK("gaussian_qpq",
             "Event outside of W2 range - W2: " + std::to_string(gamma.W2()) +
                 " outside of [" + std::to_string(W2_range_.min) + ", " +
                 std::to_string(W2_range_.max) + "]");
    return lp_gamma_event{0.};
  }

  // no generation step necessary, we just have to evaluate the cross section
  // and create the Q-Pq

  // evaluate the cross section
  const double xs_R = R(gamma.Q2());
  const double xs_dipole = dipole(gamma.Q2());
  const double xs_photo = sigma(gamma.W2());
  const double xs = (1 + gamma.epsilon() * xs_R) * xs_dipole * xs_photo;

  LOG_JUNK("gaussian_qpq",
           "xsec: " + std::to_string(xs_photo) + " < " + std::to_string(max_));
  LOG_JUNK("gaussian_qpq", "R: " + std::to_string(xs_R));
  LOG_JUNK("gaussian_qpq", "dipole: " + std::to_string(xs_dipole));

  // create our new event
  lp_gamma_event e{initial, xs, 1., xs_R};
  e.add_daughter({qpq_.type(), e.photon().p() + e.target().p(), qpq_.status()},
                 e.photon_index(), e.target_index());

  return e;
}

// =============================================================================
// gaussian_qpq::calc_max_xsec(cf)
//
// Utility function for the generator initialization
//
// max cross section as defined by the beam/target and phase-space settings
//
// Scan for cross section maximum between the W2 limits
// =============================================================================
double gaussian_qpq::calc_max_xsec(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the lepton
  // beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/dir"),
                        cf.get<double>("beam/energy")};
  const particle target{
      static_cast<pdg_id>(cf.get<int>("target/particle_type")),
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
  LOG_INFO("gaussian_qpq",
           "Max cross section: " + std::to_string(f.GetMaximum()));
  return f.GetMaximum() * 1.0000001;
}

// =============================================================================
// gaussian_qpq::calc_mass2_range()
//
// Utility function to calculate the W2 range we allow for q-Pq production
// =============================================================================
interval<double> gaussian_qpq::calc_W2_range(const double n_sigma) const {
  double min = qpq_.pole_mass() - n_sigma * qpq_.width();
  double max = qpq_.pole_mass() - n_sigma * qpq_.width();
  return {min * min, max * max};
}

// =============================================================================
// gaussian_qpq::sigma()
// gaussian_qpq::R()
// gaussian_qpq::dipole()
//
// Utility functions to calculate the cross section components
// =============================================================================
double gaussian_qpq::sigma(const double W2) const {
  return qpq_coupling_ * qpq_coupling_ * qpq_amplitude_ *
         TMath::Gaus(sqrt(W2), qpq_.mass(), qpq_.width());
}
double gaussian_qpq::R(const double Q2) const {
  return physics::R_vm_martynov(Q2, vm_pole_.mass(), R_vm_c_, R_vm_n_);
}
double gaussian_qpq::dipole(const double Q2) const {
  return physics::dipole_ff_vm(Q2, vm_pole_.mass(), dipole_n_);
}

} // namespace lp_gamma
} // namespace pcsim
