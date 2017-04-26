#ifndef PCSIM_CORE_GENERATOR_LOADED
#define PCSIM_CORE_GENERATOR_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/interval.hh>

namespace pcsim {

// =============================================================================
// Base class for all generator data
// =============================================================================
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

  static factory<process_generator> factory;

  process_generator(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}

  virtual event_type generate(const initial_type&) = 0;
};

// =============================================================================
// Base class for all event_processors (detectors/decay_handlers/...)
//
// Processor input: a valid event
// The processor will modify the input event
//
// Note: Event should derive from the event class (in core/event.hh)
// =============================================================================
template <class Event>
class event_processor : public generator<void> {
public:
  using event_type = Event;
  using base_type = generator<void>;

  process_generator(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}

  virtual void process(event_type& event) = 0;

private:
  virtual void generate() {}
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
  using base_type = generator<Event>;
  using process_type = process_generator<event_type, initial_type>;

  event_generator(const configuration& cf, const string_path& path,
                  std::shared_ptr<TRandom> r)
      : base_type{std::move(r)}, configurable{cf, path} {
    init_process_list();
  }

  virtual std::vector<event_type> generate() {
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
      // generate the sub_processes
      for (auto& process : process_list_) {
        LOG_JUNK(process.name, "Generating a trial event");
        // check if we need to generate an event for this process (ensure
        // correct sub-process mixing)
        if (rng()->Uniform(0, proc_volume_) > process.vol) {
          LOG_JUNK(process.name, "Skipping this trial cycle.");
          continue;
        }
        // generate one event
        auto event = process->generate(initial);
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
        if (rng()->Uniform(0, process.max) > event.cross_section()) {
          LOG_JUNK(process.name, "Event accepted!");
          event.update_process(process.id);
          event_list.push_back(event);
        } else {
          LOG_JUNK(process.name, "Event rejected.");
        }
      }
    } while (event_list.size() == 0);
    n_events_ += event_list.size();
    // event builder step
    for (auto& event event_list) {
      event.update_mc(n_events_, cross_section());
      build_event(event);
    }
    // that's all!
    return event_list;
  }

  // total cross section is given by the size of the generator box 
  // (phase_space * max_cross_section) times the fraction of accepted events
  // compared to the number of trials
  double cross_section() const { return volume_ * n_events_ / n_trials_; }
  double n_events() const { return n_events_; }

protected:
  // GENERATION STEPS
  // 1. generate the initial reaction, to be implemented by child class
  virtual initial_type generate_initial() const = 0;
  // 2. "event builder" step, to be implemented by child class
  virtual void build_event(event_type&) const = 0;

  // register an initial state  sub-generator (not a process generator) with
  // this event generator. This stores the relevant phase_space and
  // max_cross_section variables with the event generator
  template <class Data, class... Input>
  void register_initial(const std::shared_ptr<generator<Data, Input...>>& gen) {
    LOG_DEBUG("event_generator", "Registering phase space and cross section max");
    tassert(gen, "Requested generator is a null pointer");
    initial_ps_ *= gen->phase_space();
    initial_max_ *= gen->max_cross_section();
    update_volume();
  }

private:
  constexpr static const int N_MAX_PROC{10}; // maximum number of sub processes;
  constexpr static const char* PROC_KEY{"process_"}; // config file key
  std::string process_id(const int i) { return PROC_KEY + std::to_string(i); }

  // the maximum cross section and total phase space volume functions
  // don't make sense here, as they are different for each of the sub-processes
  virtual double max_cross_section() const { return -1; }
  virtual double phase_space() const { return -1; }

  // update the total generation volume
  void update_volume() {
    volume_ = initial_ps_ * initial_max_ * proc_volume_;
  }

  // initialize the process list
  // the factory will construct a new process generator for each of the
  // configuration file entries
  void init_process_list() {
    LOG_INFO("event_generator", "Initializing the process list");
    for (int i = 0; i < N_MAX_PROC; ++i) {
      const string_path path{PROC_KEY + std::to_string(i)};
      // check if the process is requested
      auto type = conf().get_optional<std::string>(path / "type");
      if (*type) {
        LOG_DEBUG(path.str(),
                  "Creating a new process sub-generator (" + *type + ")");
        process_list_.push_back(
            {i, FACTORY_CREATE(process_type, conf(), path, rng())});
        // check if we have a larger generation volume, update if needed
        const double volume = process_list_.back().vol;
        if (volume > proc_volume_) {
          proc_volume_ = volume;
          update_volume();
        }
      } else {
        LOG_JUNK(path.str, "Not requested");
      }
    }
    tassert(process_list_.size() > 0,
            "At least one process has to be specified");
  }

  struct process_info {
    const int id;                      // process identifier
    const std::string name;            // process name
    double ps{0};                      // process dependent phase space
    double max{0};                     // max cross section
    double vol{0};                     // generation volume
    double n_events{0};                // number of events
    std::shared_ptr<process_type> gen; // process sub-generator
    process_info(std::shared_ptr<process_type> g, const int id)
        : id{id}
        , name{process_id(id)}
        , ps{g->phase_space()}
        , max{g->max_cross_section()}
        , vol{max * ps}
        , gen{g} {}
  };

  double initial_ps_{1.};  // initial state generator phase space
  double inital_max_{1.};  // initial state generator max cross section
  double proc_volume_{-1.}; // largest generation volume in process_list
  double volume_ { 1. }    // total volume

  double n_trials_{0.};     // global trial counter, applies to all processes
  double n_tot_events_{0.}; // total number of generated events
  std::vector<process_info> process_list_; // process dependent info
};

} // namespace pcsim
#endif
