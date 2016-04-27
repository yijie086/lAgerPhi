#ifndef PCSIM_CORE_ACCEPT_REJECT_LOADED
#define PCSIM_CORE_ACCEPT_REJECT_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/interval.hh>

namespace pcsim {

// simple 1D accept reject function (TODO: generalize to N-dimensions, use
// auto-increasing fmax)
template <class Func1D>
double accept_reject_1D(std::shared_ptr<TRandom> rng,
                        const interval<double>& range, Func1D f, double fmax) {
  const double x = rng->Uniform(range.min, range.max);
  const double test = rng->Uniform(0, fmax);
  const double fx = f(x);
  // if the test fails, try until we accept a value
  if (test > fx) {
    return accept_reject_1D(rng, range, f, fmax);
  }
  // ensure a proper fmax
  tassert(fx <= fmax, "fmax set too small!!!");

  return x;
}

} // pcsim

#endif
