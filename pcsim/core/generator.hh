#ifndef PCSIM_CORE_GENERATOR_LOADED
#define PCSIM_CORE_GENERATOR_LOADED

#include <TRandom.h>
#include <algorithm>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/factory.hh>
#include <pcsim/core/interval.hh>

namespace pcsim {

// =============================================================================
// Base class for all generator data
// =============================================================================
class generator_data {
public:
  explicit generator_data(const double xs = 1.) : cross_section_{xs} {}

  double cross_section() const { return cross_section_; }
  double jacobian() const { return jacobian_; }

  void update_cross_section(const double xs) { cross_section_ *= xs; }
  void update_jacobian(const double j) { jacobian_ *= j; }

private:
  double cross_section_;
  double jacobian_{1.}; // jacobian to transform to more sane coordinate
                        // system, if necessary
};

// =============================================================================
// Base class for all generators
//
// Owns a shared pointer to the random generator.
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
// Base class for all process_generators
//
// Generator input: initial reaction information
// Generator output: a valid event
//
// Note:
//    * Event should derive from the event class (in core/event.hh)
//    * InitialData should derive from generator_data
// =============================================================================
template <class Event, class InitialData>
class process_generator : public generator<Event, InitialData> {
public:
  using event_type = Event;
  using initial_type = InitialData;
  using base_type = generator<Event, InitialData>;

  static factory<process_generator, const configuration&, const string_path&,
                 std::shared_ptr<TRandom>>
      factory;

  process_generator(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}

  virtual event_type generate(const initial_type&) = 0;
};

template <class Event, class InitialData>
factory<process_generator<Event, InitialData>, const configuration&,
        const string_path&, std::shared_ptr<TRandom>>
    process_generator<Event, InitialData>::factory;

// =============================================================================
// Base class for all event_processors (detectors/decay_handlers/...)
//
// Processor input: a valid event
// The processor will modify the input event
//
// Note: Event should derive from the event class (in core/event.hh)
// =============================================================================
template <class Event> class event_processor : public generator<void> {
public:
  using event_type = Event;
  using base_type = generator<void>;

  event_processor(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}

  virtual void process(event_type&) const = 0;

private:
  virtual void generate() {}
  virtual double max_cross_section() const { return -1.; }
  virtual double phase_space() const { return -1.; }
};

// =============================================================================
// Base class for event generators that handle the following steps:
//    * generate initial state
//    * evaluate all processes (up to 10) simulataneously
//    * accept-reject each process
//    * event building for each process
// Keeps track of the generated cross section
//
// Usage:
//    * the user should derive from this class, and provide a definition of the
//      virtual member functions
//        - InitialData generate_initial();
//        - void build_event(Event&);
//      the rest of the generation process will be handled automatically
//
// Note:
//    * Event should derive from the event class (in core/event.hh)
//    * InitialData should derive from generator_data
// =============================================================================
template <class Event, class InitialData>
class event_generator : public generator<std::vector<Event>>,
                        public configurable {
public:
  using event_type = Event;
  using initial_type = InitialData;
  using base_type = generator<std::vector<Event>>;
  using process_type = process_generator<event_type, initial_type>;

  event_generator(const configuration& cf, const string_path& path,
                  std::shared_ptr<TRandom> r)
      : base_type{std::move(r)}, configurable{cf, path} {
    init_process_list();
    init_lumi(cf);
  }

  virtual std::vector<event_type> generate() {
    // buffer for the final events we return
    // this is needed because we do one additional accept-reject step where we
    // remove events to compensate for a weight smaller than 1, which can occur
    // when, e.g., the event builder only simulates one particular decay channel
    std::vector<event_type> good_event_list;
    do {

      // generate a phase space point
      std::vector<event_type> event_list;
      do {
        n_trials_ += 1;
        auto initial = generate_initial();
        // start over if we already have a bad initial state
        if (initial.cross_section() <= 0) {
          LOG_JUNK("event_generator",
                   "Initial cross section <= 0, abandoning trial cycle.");
          continue;
        }
        LOG_JUNK("event_generator",
                 "Initial cross section: " +
                     std::to_string(initial.cross_section()));
        // generate the sub_processes
        for (auto& process : process_list_) {
          LOG_JUNK(process.name, "Generating a trial event");
          // check if we need to generate an event for this process (ensure
          // correct sub-process mixing)
          if (proc_volume_ != process.vol &&
              this->rng()->Uniform(0, proc_volume_) > process.vol) {
            LOG_JUNK(process.name, "Skipping this trial cycle.");
            continue;
          }
          // generate one event
          auto event = process.gen->generate(initial);
          // go to the next process have a bad cross section, print a warning if
          // the cross section maximum was violated
          if (event.cross_section() <= 0) {
            LOG_JUNK(process.name,
                     "Cross section <= 0, skipping this trial cycle");
            continue;
          } else if (event.cross_section() > initial_max_ * process.max) {
            LOG_WARNING(
                process.name,
                "Cross section maximum exceeded, please check "
                "the cross section maximum calculation. The MC "
                "distributions will be invalid if this happens too often.");
          }
          // accept/reject this event
          LOG_JUNK(process.name,
                   "Testing accept reject for xs: " +
                       std::to_string(event.cross_section()) + " (max: " +
                       std::to_string(process.max * initial_max_) + ")");
          if (this->rng()->Uniform(0, initial_max_ * process.max) <
              event.cross_section()) {
            LOG_JUNK(process.name, "Event accepted!");
            event.update_process(process.id);
            event_list.push_back(event);
          } else {
            LOG_JUNK(process.name, "Event rejected.");
          }
        }
      } while (event_list.empty());

      for (auto& event : event_list) {
        LOG_JUNK("generator", "Processing event (process " +
                                  std::to_string(event.process()) + ")");
        build_event(event);
        if (event.weight() > 0) {
          LOG_JUNK("generator",
                   "Event accepted after event builder step (weight: " +
                       std::to_string(event.weight()) + ", reset to 1)");
          // events are weighted with the partial branching ratio, also add
          // account for events we did not consider when only doing a limited
          // set of decay channels
          n_tot_events_ += 1 / event.weight();
          event.update_stat(static_cast<size_t>(n_tot_events_),
                            cross_section());
          event.reset_weight();
          good_event_list.push_back(event);
        } else {
          LOG_JUNK("generator",
                   "Event rejected after event builder step (weight: " +
                       std::to_string(event.weight()) + ")");
        }
      }
      // ensure we actually have an event, else start over
    } while (good_event_list.empty());

    return good_event_list;
  }

  // total cross section is given by the size of the generator box
  // (phase_space * max_cross_section) times the fraction of accepted events
  // compared to the number of trials
  //
  // THis is true for all processes, as we rescaled the number of trials T2 for
  // a process with less generation volume V2 compared to the prime process by a
  // factor of V2/V1, i.e. T2 = V2/V1 * T1
  //
  // Therefore we get that
  // sigma_2 = G2 / T2 * V2 = G2 / T1 * V1
  //
  // n_tot_events is the sum of all the generated events in all subprocesses,
  // hence we obtain sigma_tot = sum_i sigma_i = sum_i G_i * V1/T1
  double cross_section() const {
    // return a safe upper boundary in case we don't have enough events yet to
    // have some kind of reasonable estimate
    if (n_tot_events_ < 50) {
      return volume_ * process_list_.size();
    }
    // the actual cross section estimate
    return volume_ * n_tot_events_ / n_trials_;
  }
  int64_t n_events() const { return n_tot_events_; }

  // calculate the number of requested events from the lumi * cross section,
  // or alternatively use the fixed number of events
  int64_t n_requested() const {
    return (n_requested_ > 0)
               ? n_requested_
               : static_cast<int64_t>(std::round(lumi_ * cross_section()));
  }

  bool finished() const { return (n_events() >= n_requested()); }

protected:
  // GENERATION STEPS
  // 1. generate the initial reaction, to be implemented by child class
  virtual initial_type generate_initial() const = 0;
  // 2. "event builder" step, to be implemented by child class
  virtual void build_event(event_type&) const = 0;

  // register an initial state  sub-generator (not a process generator) with
  // this event generator. This stores the relevant phase_space and
  // max_cross_section variables with the event generator
  template <class InitialGen>
  void register_initial(const std::shared_ptr<InitialGen>& gen) {
    LOG_DEBUG("event_generator",
              "Registering phase space and cross section max");
    tassert(gen, "Requested generator is a null pointer");
    initial_ps_ *= gen->phase_space();
    initial_max_ *= gen->max_cross_section();
    LOG_DEBUG("event_generator",
              "New initial cross section max: " + std::to_string(initial_max_));
    LOG_DEBUG("event_generator",
              "New initial phase space: " + std::to_string(initial_ps_));
    update_volume();
  }

private:
  constexpr static const int N_MAX_PROC{10}; // maximum number of sub processes;
  constexpr static const char* PROC_KEY{"process_"}; // config file key
  static std::string process_id(const int i) {
    return PROC_KEY + std::to_string(i);
  }

  // the maximum cross section and total phase space volume functions
  // don't make sense here, as they are different for each of the sub-processes
  virtual double max_cross_section() const { return -1; }
  virtual double phase_space() const { return -1; }

  // update the total generation volume
  void update_volume() { volume_ = initial_ps_ * initial_max_ * proc_volume_; }

  // initialize the process list
  // the factory will construct a new process generator for each of the
  // configuration file entries
  void init_process_list() {
    LOG_INFO("event_generator", "Initializing the process list");
    for (int i = 0; i < N_MAX_PROC; ++i) {
      const string_path path{PROC_KEY + std::to_string(i)};
      // shortcut to avoid overload of template keywords
      const configuration& cf = conf();
      // check if the process is requested
      auto type = cf.get_optional<std::string>(path / "type");
      if (type) {
        LOG_DEBUG(path.str(),
                  "Creating a new process sub-generator (" + *type + ")");
        process_list_.push_back(
            {i, FACTORY_CREATE(process_type, cf, path, this->rng())});
        LOG_DEBUG(path.str(), "Cross section max: " +
                                  std::to_string(process_list_.back().max));
        LOG_DEBUG(path.str(),
                  "Phase space: " + std::to_string(process_list_.back().ps));
        // check if we have a larger generation volume, update if needed
        const double volume = process_list_.back().vol;
        if (volume > proc_volume_) {
          LOG_DEBUG("event_generator",
                    "Updating maximum process generation volume: " +
                        std::to_string(volume));
          proc_volume_ = volume;
          update_volume();
        }
      } else {
        LOG_JUNK(path.str(), "Not requested");
      }
    }
    tassert(process_list_.size() > 0,
            "At least one process has to be specified");
  }

  // init the number of requested events, or alternatively the requested
  // integrated luminosity
  void init_lumi(const configuration& cf) {
    auto lumi = cf.get_optional<double>("lumi");
    if (lumi) {
      lumi_ = *lumi * 1e6; // conversion from fb^-1 to nb^-1
      n_requested_ = -1;
    } else {
      n_requested_ = cf.get<int>("events");
      lumi_ = -1;
    }
  }

  struct process_info {
    const int id;                      // process identifier
    const std::string name;            // process name
    double ps{0};                      // process dependent phase space
    double max{0};                     // max cross section
    double vol{0};                     // generation volume
    double n_events{0};                // number of events
    std::shared_ptr<process_type> gen; // process sub-generator
    process_info(const int id, std::shared_ptr<process_type> g)
        : id{id}
        , name{process_id(id)}
        , ps{g->phase_space()}
        , max{g->max_cross_section()}
        , vol{max * ps}
        , gen{g} {}
  };

  double initial_ps_{1.};   // initial state generator phase space
  double initial_max_{1.};  // initial state generator max cross section
  double proc_volume_{-1.}; // largest generation volume in process_list
  double volume_{1.};       // total volume

  double n_trials_{0.};                    // global trial counter
  double n_tot_events_{0};                 // total number of events
  std::vector<process_info> process_list_; // process dependent info

  int64_t n_requested_{-1}; // number of requested events
  double lumi_{-1}; // or alternatively, the requested luminosity (in fb^-1)
};

} // namespace pcsim
#endif
