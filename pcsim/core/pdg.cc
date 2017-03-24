#include <pcsim/core/pdg.hh>
#include <TDatabasePDG.h>
#include <memory>

namespace pcsim {

// unnamed namespace
namespace {

// wrapper around ROOTs PDG database that automatically adds info on common
// nuclei (therefor avoiding the need for a custom database)
// implemented using a single persistent static database
// note: nuclear masses in GeV from
// http://hyperphysics.phy-astr.gsu.edu/hbase/pertab
class pdg_handler {
public:
  pdg_handler() : db_{make_unique<TDatabasePDG>{}} {
    // some nuclei
    db_->AddParticle("H-2", "Deuteron", 1.875613, true, 0, 1, "nucleus",
                     static_cast<int32_t>(pdg_id::H2));
    db_->AddParticle("H-3", "Triton", 2.808921, true, 0, 1, "nucleus",
                     static_cast<int32_t>(pdg_id::H3));
    db_->AddParticle("He-3", "Helium-3", 2.808391, true, 0, 1, "nucleus",
                     static_cast<int32_t>(pdg_id::He3));
    db_->AddParticle("He-4", "Helium-4", 3.727379, true, 0, 1, "nucleus",
                     static_cast<int32_t>(pdg_id::He4));
    db_->AddParticle("C-12", "Carbon-12", 11.1750, true, 0, 1, "nucleus",
                     static_cast<int32_t>(pdg_id::C12));
    db_->AddParticle("N-14", "Nitrogen-14", 13.0403, true, 0, 1, "nucleus",
                     static_cast<int32_t>(pdg_id::N14));
    // pentaquarks
    // Pc(4450)
    db_->AddParticle("Pc_4450_52p", "P_{c} (4450, 5/2+)", 4.450, false, 0.039,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4450_52p));
    db_->AddParticle("Pc_4450_52m", "P_{c} (4450, 5/2-)", 4.450, false, 0.039,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4450_52m));
    db_->AddParticle("Pc_4450_32p", "P_{c} (4450, 3/2+)", 4.450, false, 0.039,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4450_32p));
    db_->AddParticle("Pc_4450_32m", "P_{c} (4450, 3/2-)", 4.450, false, 0.039,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4450_32m));
    // Pc(4380)
    db_->AddParticle("Pc_4380_32m", "P_{c} (4380, 3/2-)", 4.380, false, 0.205,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4380_32m));
    db_->AddParticle("Pc_4380_32p", "P_{c} (4380, 3/2+)", 4.380, false, 0.205,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4380_32p));
    db_->AddParticle("Pc_4380_52m", "P_{c} (4380, 5/2-)", 4.380, false, 0.205,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4380_52m));
    db_->AddParticle("Pc_4380_52p", "P_{c} (4380, 5/2+)", 4.380, false, 0.205,
                     "pentaquark", static_cast<int32_t>(pdg_id::Pc_4380_52p));
  }
  
  TParticlePDG* find(const pdg_id id) {
    return db_->GetParticle(static_cast<int32_t>(id));
  }

private:
  std::unique_ptr<TDatabasePDG> db_;
};

// global PDG handler
pdg_handler glb_pdg;

} // unnamed namespace

TParticlePDG* pdg_particle(const pdg_id id) { return glb_pdg->find(id); }

} // namespace pcsim
