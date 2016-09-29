#ifndef PCSIM_CORE_GENERATOR_LOADED
#define PCSIM_CORE_GENERATOR_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/configuration.hh>

namespace pcsim {

// base struct to be used for generator data
struct generator_data {
  double phase_space = 1;   // phase space volume
  double cross_section = 1; // cross section
};

// Base class for all generators
//
// Owns a shared pointer to the random generator, and handles the optional
// configuration options
template <class Data, class... Input> class generator : public configurable {
public:
  using data_type = Data;

  generator(const configurable& conf, const string_path& path,
            const std::string& title, std::shared_ptr<TRandom> r)
      : configurable{conf, path}, rng_{std::move(r)} {}

  virtual data_type generate(const Input&... input) = 0;
  data_type operator()(const Input&... input) { return generate(input...); }
  std::shared_ptr<TRandom> rng() const { return rng_ }

private:
  std::shared_ptr<TRandom> rng_;
};

// Base class for event generators that handle the following steps:
//    * initial generation
//    * accept-reject
//    * event building
// Keeps track of the generated cross section
template <class Event> class event_generator : public generator<Event> {
public:
  using event_type = Event;
  using base_type = generator<Event>;

  event_generator(const configurable& conf, const string_path& path,
                  const std::string& title, std::shared_ptr<TRandom> r)
      : base_type{conf, path, title, std::move(r)}
      , phase_space_{0}
      , ntrials_{0}
      , nevents_{0}
      , max_cross_section_{0} {}

  virtual event_type generate() {
    // generate a phase space point
    event_type event;
    do {
      ntrials_ += 1;
      generate_event(event);
      tassert(event.cross_section <= max_cross_section_,
              "Invalid cross section maximum");
      // accept or reject
      double test = rng()->Uniform(0, max_cross_section_);
    } while (test > event.cross_section);
    phase_space_ += event.phase_space;
    // event builder step
    build_event(event);
    // that's all
    return event;
  }

  double cross_section() const { return phase_space_ / ntrials_; }
  double nevents() const { return nevents; }

protected:
  virtual void generate_event(event_type&) const = 0;
  virtual void build_event(event_type&) const = 0;
  void set_max_cross_section(double max) { max_cross_section_ = max; }

private:
  double phase_space_;       // cumulative phase space
  size_t ntrials_;           // total number of trials
  size_t nevents_;           // total number generated events
  double max_cross_section_; // cross section maximum
};
}

#endif
