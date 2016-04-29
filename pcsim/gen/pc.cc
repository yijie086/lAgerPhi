#include "pc.hh"
#include <pcsim/physics/decay.hh>
#include <pcsim/physics/pdg.hh>

namespace pcsim {
namespace gen {

pc::pc(const ptree& settings, const string_path& path,
       std::shared_ptr<TRandom> r)
    : base_type{settings, path, "s-channel Charmed Pentaquark Generator",
                std::move(r)}
    , xsec_{settings, path / "xsec"}
    , me_{physics::PDG_ELECTRON.Mass()}
    , Mjp_{physics::PDG_JPSI.Mass()}
    , Mp_{physics::PDG_PROTON.Mass()}
    , Wjp_{physics::PDG_JPSI_WIDTH}
    , Bje_{physics::PDG_JPSI_BRANCHING_ELEC}
    , ctheta_min_{-1.}
    , ctheta_max_{1.} {}

jpsi_event pc::gen_impl(const photon_beam& photon) {
  jpsi_event ev;

  // setup the initial state
  ev.beam.SetPxPyPzE(0, 0, photon.energy, photon.energy);
  ev.target.SetPxPyPzE(0, 0, 0, Mp_);

  // reaction in the photon-proton CM frame
  const auto cm = (ev.beam + ev.target);
  const auto beta_cm = cm.BoostVector();
  ev.s = cm.M2();
  ev.W = std::sqrt(ev.s);

  // get the cross section
  ev.xsec = xsec_(ev.W);
  // bail if the cross section is zero
  if (!(ev.xsec > 0)) {
    ev.good = false;
    return ev;
  }
  ev.weight = ev.xsec;
  ev.branching = Bje_;

  // create our Pc and boost to the lab frame
  TLorentzVector pc{0, 0, 0, ev.W};
  pc.Boost(beta_cm);

  // get the decay proton and J/Psi
  physics::decay_pc_iso(rng(), pc, ev.recoil, ev.jpsi);

  // get the J/Psi decay leptons
  physics::decay_jpsi_lepton(rng(), ev.jpsi, me_, ev.positron, ev.electron);

  // all done!
  ev.good = true;
  return ev;
}

} // gen
} // pcsim
