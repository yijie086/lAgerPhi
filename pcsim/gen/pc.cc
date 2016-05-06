#include "pc.hh"
#include <pcsim/core/assert.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/physics/constants.hh>
#include <pcsim/physics/decay.hh>
#include <pcsim/physics/pdg.hh>

namespace pcsim {
namespace gen {

pc::pc(const configuration& conf, const string_path& path,
       std::shared_ptr<TRandom> r)
    : base_type{conf, path, "s-channel Charmed Pentaquark Generator",
                std::move(r)}
    , brems_{conf, path / "photon_beam"}
    , xsec_{conf, path / "xsec"}
    , xsec_max_{xsec_.max() * brems_.max()}
    , s_range_{brems_.calc_s_range()} {
  LOG_INFO("pc", "s range: [" + std::to_string(s_range_.min) + ", " +
                     std::to_string(s_range_.max) + "] GeV^2");
}

// generate a good e-p->gamma,p->Pc->J/Psi,p->e+e-,p event
jpsi_event pc::gen_impl() {

  jpsi_event ev;

  // generate a new phase space point using the cross section
  double xsec = 0; // cross section for the generated phase space point
  double test = 0; // test variable used by accept-reject
  do {
    ev.s = rng()->Uniform(s_range_.min, s_range_.max);
    ev.W = std::sqrt(ev.s);

    // get the cross section and the corresponding bremsstrahlung intensity
    xsec = xsec_(ev.s) * brems_(ev.s);

    // record the phase space volume
    ev.volume = s_range_.width() * xsec_max_;

    // record this attempt to generate a new event
    ev.ntrials += 1;

    // accept or reject!
    test = rng()->Uniform(0, xsec_max_);
    tassert(xsec <= xsec_max_, "Invalid cross section maximum!");

  } while (test > xsec);

  // setup the initial state
  const double Egamma =
      (ev.s - physics::M2_PROTON) * physics::ONE_OVER_2M_PROTON;
  ev.beam.SetPxPyPzE(0, 0, Egamma, Egamma);
  ev.target.SetPxPyPzE(0, 0, 0, physics::M_PROTON);

  // reaction in the photon-proton CM frame
  const auto cm = (ev.beam + ev.target);
  const auto beta_cm = cm.BoostVector();

  // create our Pc and boost to the lab frame
  TLorentzVector pc{0, 0, 0, ev.W};
  pc.Boost(beta_cm);

  // get the decay proton and J/Psi
  physics::decay_pc_iso(rng(), pc, ev.recoil, ev.jpsi);

  // get t from the J/Psi
  ev.t = (ev.jpsi - ev.beam).M2();

  // get the decay leptons
  physics::decay_jpsi_lepton(rng(), ev.jpsi, physics::M_ELECTRON, ev.positron,
                             ev.electron);

  // all done!
  return ev;
}

} // gen
} // pcsim
