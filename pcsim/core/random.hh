#ifndef PCSIM_CORE_RANDOM_LOADED
#define PCSIM_CORE_RANDOM_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/interval.hh>

namespace pcsim {

// =============================================================================
// Generate a random value along an arbitrary distribution (has to be positive
// in the interval; normalization is not needed)
//
// Impelemented using an accept-reject algorithm
// =============================================================================
template <class Func1D>
inline random_func(std::shared_ptr<TRandom> rng, const interval<double>& range,
                   Func1D f, const double fmax) {
  double x = 0;
  double test = -1;
  double fx = -1;
  do {
    x = rng->Uniform(range.min, range.max);
    test = rng->Uniform(0, fmax);
    fx = f(x);
    // ensure a proper fmax
    tassert(fx <= fmax, "fmax set too small in random_func call");
    // if the test fails, try again
  } while (test > fx);

  // all done
  return x;
}

#endif
