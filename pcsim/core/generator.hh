#ifndef PCSIM_CORE_GENERATOR_LOADED
#define PCSIM_CORE_GENERATOR_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/interval.hh>

namespace pcsim {

// =============================================================================
// Base class for all generators
//
// Owns a shared pointer to the random generator, and handles the optional
// configuration options
//
// Note: Data should be a struct with at least a cross_section data member
// =============================================================================
template <class Data, class... Input> class generator : public configurable {
public:
  using data_type = Data;

  generator(const configuration& conf, const string_path& path,
            std::shared_ptr<TRandom> r)
      : configurable{conf, path}, rng_{std::move(r)} {}

  virtual data_type generate(const Input&... input) = 0;
  data_type operator()(const Input&... input) { return generate(input...); }

  // contribution to the maximum cross section
  virtual double max_cross_section() const { return 1; }

protected:
  // access the RNG
  std::shared_ptr<TRandom> rng() const { return rng_; }

  // generate a random number following an arbitrary function
  // parameters:
  //  range: generation range
  //  f: arbitrary function (continous within our range)
  //  fmax: maximum of f within our range
  template <class Func1D>
  double rand_f(const interval<double>& range, Func1D f, double fmax) const {
    const double x = rng()->Uniform(range.min, range.max);
    const double test = rng()->Uniform(0, fmax);
    const double fx = f(x);
    // if the test fails, try again (accept-reject)
    if (test > fx) {
      return rand_f(range, f, fmax);
    }
    // ensure a proper fmax
    tassert(fx <= fmax, "fmax set too small in rand_f call");

    return x;
  }

private:
  std::shared_ptr<TRandom> rng_;
};

// =============================================================================
// Base class for event generators that handle the following steps:
//    * initial generation
//    * accept-reject
//    * event building
// Keeps track of the generated cross section
// =============================================================================
template <class Event> class event_generator : public generator<Event> {
public:
  using event_type = Event;
  using base_type = generator<Event>;

  event_generator(const configuration& conf, const string_path& path,
                  std::shared_ptr<TRandom> r)
      : base_type{conf, path, std::move(r)}
      , phase_space_{0}
      , ntrials_{0}
      , nevents_{0}
      , max_cross_section_{1} {}

  virtual event_type generate() {
    // generate a phase space point
    event_type event;
    double test;
    do {
      ntrials_ += 1;
      event.part.clear();
      generate_event(event);
      tassert(event.cross_section <= max_cross_section_,
              "Invalid cross section maximum");
      // accept or reject
      test = this->rng()->Uniform(0, max_cross_section_);
    } while (test > event.cross_section);
    phase_space_ += event.phase_space * max_cross_section_;
    nevents_ += 1;
    // event builder step
    build_event(event);
    // that's all
    return event;
  }

  double cross_section() const { return phase_space_ / ntrials_; }
  double nevents() const { return nevents_; }

protected:
  virtual void generate_event(event_type&) const = 0;
  virtual void build_event(event_type&) const = 0;
  void add_max_cross_section(double max) { max_cross_section_ *= max; }
  double max_cross_section() const { return max_cross_section_; }

private:
  double phase_space_;       // cumulative phase space
  size_t ntrials_;           // total number of trials
  size_t nevents_;           // total number generated events
  double max_cross_section_; // cross section maximum
};

}
#endif
