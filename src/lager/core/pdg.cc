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

#include <TDatabasePDG.h>
#include <lager/core/pdg.hh>
#include <memory>

namespace lager {

// unnamed namespace
namespace {

// wrapper around ROOTs PDG database that automatically adds info on common
// nuclei (therefor avoiding the need for a custom database)
// implemented using a single persistent static database
// note: nuclear masses in GeV from
// http://hyperphysics.phy-astr.gsu.edu/hbase/pertab
class pdg_handler {
public:
  pdg_handler() : db_{std::make_unique<TDatabasePDG>()} {
    // some nuclei
    db_->AddParticle("H-2", "Deuteron", 1.875613, true, 0, 3, "nucleus",
                     static_cast<int32_t>(pdg_id::H2));
    db_->AddParticle("H-3", "Triton", 2.808921, true, 0, 3, "nucleus",
                     static_cast<int32_t>(pdg_id::H3));
    db_->AddParticle("He-3", "Helium-3", 2.808391, true, 0, 6, "nucleus",
                     static_cast<int32_t>(pdg_id::He3));
    db_->AddParticle("He-4", "Helium-4", 3.727379, true, 0, 6, "nucleus",
                     static_cast<int32_t>(pdg_id::He4));
    db_->AddParticle("C-12", "Carbon-12", 11.1750, true, 0, 18, "nucleus",
                     static_cast<int32_t>(pdg_id::C12));
    db_->AddParticle("N-14", "Nitrogen-14", 13.0403, true, 0, 21, "nucleus",
                     static_cast<int32_t>(pdg_id::N14));
    db_->AddParticle("Al-27", "Aluminum-27", 25.1166, true, 0, 21, "nucleus",
                     static_cast<int32_t>(pdg_id::Al27));
    // pentaquarks (mass and width not used internally as this changes for each
    // assumption)
    db_->AddParticle("Pc_wang_52p", "P_{c} (5/2+)", 0, false, 0, 3,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_wang_52p));
    db_->AddParticle("Pc_wang_52m", "P_{c} (5/2-)", 0, false, 0, 3,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_wang_52m));
    db_->AddParticle("Pc_wang_32p", "P_{c} (3/2+)", 0, false, 0, 3,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_wang_32p));
    db_->AddParticle("Pc_wang_32m", "P_{c} (3/2-)", 0, false, 0, 3,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_wang_32m));
    db_->AddParticle("Pc_iso_52p", "P_{c} (5/2+)", 0, false, 0, 3, "pentaquark",
                     static_cast<int32_t>(pdg_id::Pc_iso_52p));
    db_->AddParticle("Pc_iso_52m", "P_{c} (5/2-)", 0, false, 0, 3, "pentaquark",
                     static_cast<int32_t>(pdg_id::Pc_iso_52m));
    db_->AddParticle("Pc_iso_32p", "P_{c} (3/2+)", 0, false, 0, 3, "pentaquark",
                     static_cast<int32_t>(pdg_id::Pc_iso_32p));
    db_->AddParticle("Pc_iso_32m", "P_{c} (3/2-)", 0, false, 0, 3, "pentaquark",
                     static_cast<int32_t>(pdg_id::Pc_iso_32m));
    // Unknown
    db_->AddParticle("Unknown", "Unknown", 0, true, 0, 0, "Unknown",
                     static_cast<int32_t>(pdg_id::unknown));
  }

  TParticlePDG* find(const pdg_id id) {
    return db_->GetParticle(static_cast<int32_t>(id));
  }
  TParticlePDG* find(const std::string& name) {
    return db_->GetParticle(name.c_str());
  }

private:
  std::unique_ptr<TDatabasePDG> db_;
};

// global PDG handler
pdg_handler glb_pdg;

} // unnamed namespace

TParticlePDG* pdg_particle(const pdg_id id) { return glb_pdg.find(id); }
TParticlePDG* pdg_particle(const std::string& name) {
  return glb_pdg.find(name);
}

} // namespace lager
