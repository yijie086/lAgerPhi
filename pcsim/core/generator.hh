#ifndef PCSIM_CORE_GENERATOR_LOADED
#define PCSIM_CORE_GENERATOR_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/interval.hh>

namespace pcsim {

class generator_data {
public:
  explicit generator_data(const double xs = 1.) : cross_section_{xs} {}

  void update_cross_section(const double xs) { cross_section_ *= xs; }
  double cross_section() const { return cross_section_; }

private:
  double cross_section_;
};

// =============================================================================
// Base class for all generators
//
// Owns a shared pointer to the random generator, and handles the optional
// configuration options
// =============================================================================
template <class Data, class... Input> class generator {
public:
  using data_type = Data;

  generator(std::shared_ptr<TRandom> r) : rng_{std::move(r)} {}

  // generate an event
  virtual data_type generate(const Input&... input) = 0;
  data_type operator()(const Input&... input) { return generate(input...); }

  // contribution to the maximum cross section and phase space volume
  // (to be implemented by child class)
  virtual double max_cross_section() const = 0;
  virtual double phase_space() const = 0;

protected:
  // access the RNG
  std::shared_ptr<TRandom> rng() const { return rng_; }

  // generate a random number following an arbitrary function
  // parameters:
  //  range: generation range
  //  f: arbitrary function (continous within our range)
  //  fmax: maximum of f within our range
  // this is useful for e.g. azimuthal distributions that don't impact the total
  // cross section estimate
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
//
// Note: Event should derive from the event class (in core/event.hh)
// =============================================================================
template <class Event>
class event_generator : public generator<Event>, public configurable {
public:
  using event_type = Event;
  using base_type = generator<Event>;

  event_generator(const configuration& conf, const string_path& path,
                  std::shared_ptr<TRandom> r)
      : base_type{std::move(r)}, configurable{conf, path} {}

  virtual event_type generate() {
    // generate a phase space point
    event_type event;
    double test;
    do {
      n_trials_ += 1;
      event = {};
      generate_event(event);
      tassert(event.cross_section() <= max_cross_section_,
              "Invalid cross section maximum");
      // accept or reject
      test = this->rng()->Uniform(0, max_cross_section_);
    } while (test > event.cross_section());
    n_events_ += 1;
    // event builder step
    build_event(event);
    // that's all
    return event;
  }

  // total cross section is given by the size of the generator box 
  // (phase_space * max_cross_section) times the fraction of accepted events
  // compared to the number of trials
  double cross_section() const {
    return phase_space_ * max_cross_section_ * n_events_ / n_trials_;
  }
  double n_events() const { return n_events_; }

  // get the maximum cross section and total phase space volume
  virtual double max_cross_section() const { return max_cross_section_; }
  virtual double phase_space() const { return phase_space_; }

protected:
  // actual event generation step, to be implemented by child class
  virtual void generate_event(event_type&) const = 0;
  // "event builder" step, to be implemented by child class
  virtual void build_event(event_type&) const = 0;

  // register a sub-process generator with this event generator.
  // This stores the relevant phase_space and max_cross_section variables with
  // the event generator
  template <class Data, class... Input>
  void add(const generator<Data, Input...>& gen) {
    max_cross_section_ *= gen.max_cross_section();
    phase_space_ *= gen.phase_space();
  }
  void update_max_cross_section(double max) { max_cross_section_ *= max; }
  void update_phase_space(const double ps) { phase_space_ *= ps; }

private:
  double phase_space_{1.};       // phase space
  size_t n_events_{0};           // total number generated events
  size_t n_trials_{0};           // total number of trials
  double max_cross_section_{1.}; // cross section maximum
};

} // namespace pcsim
#endif
