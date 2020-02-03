#ifndef LIEGE_GEN_LP_GAMMA_GENERATOR_LOADED
#define LIEGE_GEN_LP_GAMMA_GENERATOR_LOADED

#include <liege/core/factory.hh>
#include <liege/core/generator.hh>
#include <liege/gen/lp_gamma_event.hh>

// =============================================================================
// main include file for lp_gamma generators
// =============================================================================
 
namespace liege {
namespace lp_gamma {

// =============================================================================
// lp_gamma generator type
// =============================================================================
using generator = liege::process_generator<lp_gamma_event, lp_gamma_data>;

} // namespace beam
} // namespace liege

#endif
