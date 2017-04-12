#include "gamma_p_2VmX.hh"
#include <pcsim/physics/photon.hh>
#include <pcsim/physics/vm.hh>

namespace pcsim {
namespace process {

gamma_p_2VmX_data(const beam::photon_data& photon, const beam::data& target,
                  const double t, const particle& vm1, const particle& X1,
                  const double xs, const double phi, const double R,
                  const double epsilon)
    : generator_data{xs}
    , t_{t}
    , xv_{photon.x() + vm1.mass2() / (2 * target.beam().mass())}
    , Q2plusMv2_{photon.Q2() + vm1.mass2()}
    , R_{R}
    , epsilon_{epsilon}
    , vm_{vm1}
    , X_{X1} {

  // utility shortcuts
  const double W2 = photon.W2();
  const double W = sqrt(W2);
  const double Q2 = photon.Q2();
  const double y = photon.y();
  const double nu = photon.nu();
  const double Mt2 = target.beam().mass2();
  const double Mr2 = X_.mass2();
  const double Mv2 = vm_.mass2();

  // create our final state particles in the CM frame
  // energies and momenta
  const double Et_cm = (W2 + Q2 + Mt2) / (2. * W);
  const double Pt_cm = sqrt(Et_cm * Et_cm - Mt2);
  const double Er_cm = (W2 - Mv2 + Mr2) / (2. * W);
  const double Pr_cm = sqrt(Er_cm * Er_cm - Mr2);
  const double Ev_cm = (W2 + Mv2 - Mr2) / (2. * W);
  const double Pv_cm = sqrt(Ev_cm * Ev_cm - Mv2);

  // get the VM and recoil CM 4-vectors
  // NOTE:
  //  * theta is the change in angle from the initial state, don't forget the
  //    target originally was flying backwards (already has theta=pi)
  // calculate the scattering angle theta
  const double ctheta_cm =
      (t + 2 * Et_cm * Er_cm - Mt2 - Mr2) / (2 * Pt_cm * Pr_cm);
  const double theta_cm = std::acos(ctheta_cm);
  const double phi_cm = phi;

  // create the CM 4-vectors
  const particle::Polar3DVector p3_v{Pv_cm, theta_cm, phi_cm};
  vm_.p() = {p3_v.X(), p3_v.Y(), p3_v.Z(), Ev_cm};
  const particle::Polar3DVector p3_r{Pr_cm, theta_cm + TMath::Pi(), phi_cm};
  X_.p() = { p3_r.X(), p3_r.Y(), p3_r.Z(), er_cm;

  // calculate some necessary boost and rotation vectors
  particle targ = target.beam();
  particle phot = photon.beam();
  // lab frame to target rest frame (trf)
  particle::Boost boost_trf{targ.BoostToCM()};
  targ.boost(boost_trf);
  phot.boost(boost_trf);
  // describe photon in a target rest frame where the photon moves along the
  // z-axis
  particle::XYZTVector p_phot_z{0, 0, phot.momentum(), phot.energy()};
  // boost to CM frame from this rotated trf
  particle::boost boost_cm{(targ.p() + p_phot_z).BoostToCM()};
  // 1. CM -> rotated TRF
  vm_.boost(-boost_cm);
  X_.boost(-boost_cm);
  // 2. rotated TRF -> TRF
  vm_.rotate_uz(phot.p());
  X_.rotate_uz(phot.p());
  // 3. TRF -> lab
  vm.boost(-boost_trf);
  x_.boost(-boost_trf);
  // all done!
}

} // namespace process
} // namespace pcsim
