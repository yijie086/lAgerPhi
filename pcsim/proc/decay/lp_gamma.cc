#include "lp_gamma.hh"
#include <cmath>
#include <pcsim/core/particle.hh>
#include <pcsim/core/pdg.hh>
#include <pcsim/core/stringify.hh>

#include <HepMC/GenEvent.h>
#include <photospp/Photos.h>
#include <photospp/PhotosHepMCEvent.h>

namespace pcsim {
namespace decay {

lp_gamma::lp_gamma(const configuration& conf, const string_path& path,
                   std::shared_ptr<TRandom> r)
    : lp_gamma::base_type{std::move(r)}
    , vm_decay_lplus_{static_cast<pdg_id>(
          -abs(conf.get<int>(path / "vm_decay_lepton_type", 11)))}
    , vm_decay_lminus_{static_cast<pdg_id>(
          abs(conf.get<int>(path / "vm_decay_lepton_type", 11)))}
    , vm_decay_br_{conf.get<double>(path / "vm_branching_ratio", 1)} {
  LOG_INFO("decay", "VM decays into " + vm_decay_lplus_.name() +
                        vm_decay_lminus_.name());
  LOG_INFO("decay",
           "VM Branching ratio set to: " + std::to_string(vm_decay_br_));
  auto do_radiative_decay_vm =
      conf.get<bool>(path / "do_radiative_decay_vm", false);
  if (do_radiative_decay_vm) {
    LOG_INFO("decay", "Simulating radiative decay for VM particles");
    radiative_decay_ = std::make_unique<radiative_decay_vm>();
  }
}

void lp_gamma::process(lp_gamma_event& e) const {
  for (int i = 0; i < e.size(); ++i) {
    LOG_JUNK2("decay::lp_gamma",
              "Considering decay for particle " + e[i].name());
    // we won't decay stable particles
    if (e[i].stable()) {
      LOG_JUNK2("decay::lp_gamma", "Particle does not need to be decayed");
      continue;
    }
    // SCHC leptonic e+e- decay of vms
    if (e[i].type() == pdg_id::J_psi || e[i].type() == pdg_id::upsilon ||
        e[i].type() == pdg_id::phi) {
      LOG_JUNK2("decay::lp_gamma", "Performing SCHC decay for VMs");
      quarkonium_schc(e, i);
    }
    // Pc decay according to wang et al.
    else if (e[i].type() == pdg_id::Pc_4450_52p ||
             e[i].type() == pdg_id::Pc_4380_52p ||
             e[i].type() == pdg_id::Pc_4450_52m ||
             e[i].type() == pdg_id::Pc_4380_52m ||
             e[i].type() == pdg_id::Pc_4450_32p ||
             e[i].type() == pdg_id::Pc_4380_32p ||
             e[i].type() == pdg_id::Pc_4450_32m ||
             e[i].type() == pdg_id::Pc_4380_32m) {
      LOG_JUNK2("decay::lp_gamma", "Pc decay according to Wang et. al.");
      pentaquark_wang(e, i);
    } else {
      LOG_DEBUG("decay::lp_gamma", "Unstable particle " + e[i].name() +
                                       ", but no decay path implemented");
    }
  }
  // that's all
}

void lp_gamma::quarkonium_schc(lp_gamma_event& e, const int i) const {
  // electron or muon BR only
  e.update_weight(vm_decay_br_);
  std::pair<particle, particle> decay_products{{vm_decay_lplus_.type()},
                                               {vm_decay_lminus_.type()}};
  std::pair<particle, particle> decay_products_cm{
      {vm_decay_lplus_.type(), particle::status_code::INFO_PARENT_CM},
      {vm_decay_lminus_.type(), particle::status_code::INFO_PARENT_CM}};
  const double epsilon_R = e.epsilon() * e.R();
  const double r04 = epsilon_R / (1 + epsilon_R);
  const double phi = rng()->Uniform(0., TMath::TwoPi());
  const double ctheta = rand_f(
      {-1, 1},
      [=](const double ctheta) {
        return ((1. + r04) + (1. - 3. * r04) * ctheta * ctheta);
      },
      2.001);
  const double theta = acos(ctheta);
  physics::decay_2body(e[i], theta, phi, decay_products, decay_products_cm);
  // add the decay particles
  e.add_daughter(decay_products, i);
  if (radiative_decay_) {
    radiative_decay_->process(e, i);
  }

  // add the SCHC info in the VM helicity frame
  decay_products_cm.first.add_parent(i);
  decay_products_cm.second.add_parent(i);
  e.add_particle(decay_products_cm);

  // mark the vm as decayed
  e[i].update_status(particle::status_code::DECAYED_SCHC);
}
radiative_decay_vm::radiative_decay_vm() {
  Photospp::Photos::initialize();
  // 0.5MeV cutoff for a 3GeV J/psi
  Photospp::Photos::setInfraredCutOff(.0005 / 3.);
}
void radiative_decay_vm::process(lp_gamma_event& e, const int vm_index) {
  const std::pair<int, int> decay_index = {e[vm_index].daughter_begin(),
                                           e[vm_index].daughter_begin() + 1};
  HepMC::GenEvent evt(20, 1);
  evt.use_units(HepMC::Units::GEV, HepMC::Units::CM);
  auto vx = std::make_shared<HepMC::GenVertex>();
  evt.add_vertex(vx.get());
  vx->add_particle_in(new HepMC::GenParticle(
      HepMC::FourVector(e[vm_index].p().X(), e[vm_index].p().Y(),
                        e[vm_index].p().Z(), e[vm_index].p().E()),
      e[vm_index].type<int>(), 31));
  vx->add_particle_out(new HepMC::GenParticle(
      HepMC::FourVector(
          e[decay_index.first].p().X(), e[decay_index.first].p().Y(),
          e[decay_index.first].p().Z(), e[decay_index.first].p().E()),
      e[decay_index.first].type<int>(), 31));
  vx->add_particle_out(new HepMC::GenParticle(
      HepMC::FourVector(
          e[decay_index.second].p().X(), e[decay_index.second].p().Y(),
          e[decay_index.second].p().Z(), e[decay_index.second].p().E()),
      e[decay_index.second].type<int>(), 31));
  Photospp::PhotosHepMCEvent photos_event(&evt);
  photos_event.process();
  // did we radiate one (or more) photons?
  if (evt.particles_size() > 3) {
    // then we need to update our event record and add the photon
    auto it = evt.particles_begin();
    ++it; // skip the VM
    // next up is the first decay lepton
    e[decay_index.first].p() = (*it)->momentum();
    ++it;
    // second decay lepton
    e[decay_index.second].p() = (*it)->momentum();
    ++it;
    // now store the photons
    for (; it != evt.particles_end(); ++it) {
      e.add_daughter({static_cast<pdg_id>((*it)->pdg_id()),
                      particle::XYZTVector((*it)->momentum())},
                     vm_index);
    }
  }
}
void lp_gamma::pentaquark_wang(lp_gamma_event& e, const int i) const {
  std::pair<particle, particle> decay_products{
      {pdg_id::J_psi, particle::status_code::UNSTABLE_SCHC}, {pdg_id::p}};
  std::pair<particle, particle> decay_products_cm{
      {pdg_id::J_psi, particle::status_code::INFO_PARENT_CM},
      {pdg_id::p, particle::status_code::INFO_PARENT_CM}};
  const double phi = rng()->Uniform(0., TMath::TwoPi());
  double ctheta = -1;
  if (e[i].type() == pdg_id::Pc_4450_52p ||
      e[i].type() == pdg_id::Pc_4380_52p) {
    // result from a pol6 fit to a digitized version of figure 6c from
    // PRD92-034022(2015)
    ctheta = rand_f(
        {-1, 1},
        [](const double x) {
          const double x2 = x * x;
          const double x3 = x2 * x;
          const double x4 = x3 * x;
          const double x5 = x4 * x;
          const double x6 = x5 * x;
          return .149211 - 0.194418 * x - 0.563191 * x2 + 0.374024 * x3 +
                 0.658942 * x4 + 0.110057 * x5 + 0.0931712 * x6;
        },
        0.63);
  } else if (e[i].type() == pdg_id::Pc_4450_52m ||
             e[i].type() == pdg_id::Pc_4380_52m) {
    // result from a pol7 fit to a digitized version of figure 5c from
    // PRD92-034022(2015)
    ctheta = rand_f(
        {-1, 1},
        [](const double x) {
          const double x2 = x * x;
          const double x3 = x2 * x;
          const double x4 = x3 * x;
          const double x5 = x4 * x;
          const double x6 = x5 * x;
          const double x7 = x6 * x;
          return 1.31241 - 1.19802 * x + 1.58351 * x2 + 17.1514 * x3 +
                 20.8306 * x4 - 4.43848 * x5 + 2.67151 * x6 + 6.06378 * x7;
        },
        44.06);
  } else if (e[i].type() == pdg_id::Pc_4450_32p ||
             e[i].type() == pdg_id::Pc_4380_32p) {
    // result from a expo fit to a digitized version of figure 5b from
    // PRD92-034022(2015)
    ctheta = rand_f(
        {-1, 1}, [](const double x) { return exp(-5.944 - x); }, 0.00713);
  } else if (e[i].type() == pdg_id::Pc_4450_32m ||
             e[i].type() == pdg_id::Pc_4380_32m) {
    // result from a pol2 fit to a digitized version of figure 6b from
    // PRD92-034022(2015)
    ctheta = rand_f(
        {-1, 1},
        [](const double x) {
          const double x2 = x * x;
          return 0.00845846 - 0.0128146 * x + 0.00526053 * x2;
        },
        0.0266);
  }
  const double theta = acos(ctheta);
  physics::decay_2body(e[i], theta, phi, decay_products, decay_products_cm);
  // add the decay particles as leading/recoil if not already set
  if (e.leading_index() < 0) {
    e.add_leading(decay_products.first, i);
    e.add_recoil(decay_products.second, i);
  } else {
    e.add_daughter(decay_products, i);
  }

  // also add the Pc CM info particles
  decay_products_cm.first.add_parent(i);
  decay_products_cm.second.add_parent(i);
  e.add_particle(decay_products_cm);

  // mark the Pc as decayed
  e[i].update_status(particle::status_code::DECAYED);
} // namespace decay

} // namespace decay
} // namespace pcsim
