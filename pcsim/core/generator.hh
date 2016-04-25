#ifndef PCSIM_GENERATOR_LOADED
#define PCSIM_GENERATOR_LOADED

#include <memory>
#include <TRandom.h>

namespace pcsim {

// Base class for all generators. Using a recursive template implementation for
// maximum runtime performace.
template <class Derived>
class generator {
public:
  using event_type = typename Derived::event_type;

  template <class... Input> event_type generate(const Input&... input) {
    return derived().gen_impl(input...);
  }
  template <class... Input> event_type operator()(const Input&... input) {
    return generate(input);
  }

protected:
  std::shared_ptr<TRandom> rng_;

private:
  Derived& derived() { return static_cast<Derived>(*this); }
  const Derived& derived() const { return static_cast<const Derived>(*this); }
};

}

#endif
