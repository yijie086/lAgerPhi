// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
//
// This file is part of lAger.
//
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
//

#include "resonance_qpq.hh"
#include <TF1.h>
#include <TMath.h>
#include <lager/core/logger.hh>
#include <lager/gen/initial/target_gen.hh>
#include <lager/physics/photon.hh>
#include <lager/physics/vm.hh>

namespace lager {
namespace lA {

// =============================================================================
// Constructor for lA::resonance_qpq
// =============================================================================
resonance_qpq::resonance_qpq(const configuration& cf, const string_path& path,
                             std::shared_ptr<TRandom> r)
    : base_type{r}
    , vm_pole_{cf.get<std::string>(path / "vm_type")}
    , qpq_{cf.get<std::string>(path / "qpq_type"),
           particle::status_code::UNSTABLE}
    , mass_{cf.get<double>(path / "mass")}
    , width_{cf.get<double>(path / "width")}
    , amplitude_{cf.get<double>(path / "amplitude")}
    , norm_{calc_normalization()}
    , coupling_{cf.get<double>(path / "coupling")}
    , R_vm_c_{cf.get<double>(path / "R_vm_c")}
    , R_vm_n_{cf.get<double>(path / "R_vm_n")}
    , dipole_n_{cf.get<double>(path / "dipole_n")}
    , W2_range_{calc_W2_range(cf.get<double>(path / "n_sigma"))}
    , max_{calc_max_xsec(cf)} {
  LOG_INFO("resonance_qpq", "amplitude: " + std::to_string(amplitude_));
  LOG_INFO("resonance_qpq", "width: " + std::to_string(width_));
  LOG_INFO("resonance_qpq", "norm: " + std::to_string(norm_));
  LOG_INFO("resonance_qpq", "coupling: " + std::to_string(coupling_));
  LOG_INFO("resonance_qpq", "R_vm c-parameter: " + std::to_string(R_vm_c_));
  LOG_INFO("resonance_qpq",
           "R_vm n-parameter (power): " + std::to_string(R_vm_n_));
  LOG_INFO("resonance_qpq", "'Dipole' FF power: " + std::to_string(dipole_n_));
  LOG_INFO("resonance_qpq",
           "Q-Pentaquark: " + std::string(qpq_.pdg()->GetName()));
  LOG_INFO("resonance_qpq",
           "VM Pole: " + std::string(vm_pole_.pdg()->GetName()));
}

lA_event resonance_qpq::generate(const lA_data& initial) {
  const auto& gamma = initial.photon();

  // check if we are in the correct W2 range
  if (W2_range_.excludes(gamma.W2())) {
    LOG_JUNK("resonance_qpq",
             "Event outside of W2 range - W2: " + std::to_string(gamma.W2()) +
                 " outside of [" + std::to_string(W2_range_.min) + ", " +
                 std::to_string(W2_range_.max) + "]");
    return lA_event{0.};
  }

  // no generation step necessary, we just have to evaluate the cross section
  // and create the Q-Pq

  // evaluate the cross section
  const double xs_R = R(gamma.Q2());
  const double xs_dipole = dipole(gamma.Q2());
  const double xs_photo = sigma(gamma.W2());
  const double xs = (1 + gamma.epsilon() * xs_R) * xs_dipole * xs_photo;

  LOG_JUNK("resonance_qpq",
           "xsec: " + std::to_string(xs_photo) + " < " + std::to_string(max_));
  LOG_JUNK("resonance_qpq", "R: " + std::to_string(xs_R));
  LOG_JUNK("resonance_qpq", "dipole: " + std::to_string(xs_dipole));

  // create our new event
  lA_event e{initial, xs, 1., xs_R};
  int qpq_idx = e.add_daughter(
      {qpq_.type(), e.photon().p() + e.target().p(), qpq_.status()},
      e.photon_index(), e.target_index());
  e[qpq_idx].vertex() = e.photon().vertex();

  return e;
}
// =============================================================================
// resonance_qpq::calc_normalization()
//
// Utility function for the generator initialization
//
// Calculate the normalization factor needed to get a BW shape with the
// desired "amplitude".
// =============================================================================
double resonance_qpq::calc_normalization() const {
  TF1 f(
      "fraw", [=](double* W, double*) { return sigma(W[0] * W[0]); }, 0,
      mass_ - 3 * width_, mass_ + 3 * width_);
  LOG_INFO("resonance_qpq",
           "Raw BW maximum found: " + std::to_string(f.GetMaximum()));
  const double norm = amplitude_ / f.GetMaximum();
  LOG_INFO("resonance_qpq",
           "Calculated normalization to reach desired amplitude: " +
               std::to_string(norm));
  return norm;
}

// =============================================================================
// resonance_qpq::calc_max_xsec(cf)
//
// Utility function for the generator initialization
//
// max cross section as defined by the beam/target and phase-space settings
//
// Scan for cross section maximum between the W2 limits
// =============================================================================
double resonance_qpq::calc_max_xsec(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the lepton
  // beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{initial::estimated_target(cf)};
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
  LOG_INFO("resonance_qpq",
           "Max cross section: " + std::to_string(f.GetMaximum()));
  return f.GetMaximum() * 1.0000001;
}

// =============================================================================
// resonance_qpq::calc_mass2_range()
//
// Utility function to calculate the W2 range we allow for q-Pq production
// =============================================================================
interval<double> resonance_qpq::calc_W2_range(const double n_sigma) const {
  double min = mass_ - n_sigma * width_ / 2.;
  double max = mass_ + n_sigma * width_ / 2.;
  return {min * min, max * max};
}

// =============================================================================
// resonance_qpq::sigma()
// resonance_qpq::R()
// resonance_qpq::dipole()
//
// Utility functions to calculate the cross section components
// =============================================================================
double resonance_qpq::sigma(const double W2) const {
  return coupling_ * coupling_ * norm_ *
         TMath::BreitWigner(sqrt(W2), mass_, width_);
}
double resonance_qpq::R(const double Q2) const {
  return physics::R_vm_martynov(Q2, vm_pole_.mass(), R_vm_c_, R_vm_n_);
}
double resonance_qpq::dipole(const double Q2) const {
  return physics::dipole_ff_vm(Q2, vm_pole_.mass(), dipole_n_);
}

} // namespace lA
} // namespace lager
