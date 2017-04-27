#ifndef PCSIM_GEN_LP_GAMMA_GENERATOR_LOADED
#define PCSIM_GEN_LP_GAMMA_GENERATOR_LOADED

#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/gen/lp_gamma_event.hh>

// =============================================================================
// main include file for lp_gamma generators
// =============================================================================
 
namespace pcsim {
namespace lp_gamma {

// =============================================================================
// lp_gamma generator type
// =============================================================================
using generator = process_generator<lp_gamma_event, lp_gamma_data>;

} // namespace beam
} // namespace pcsim
