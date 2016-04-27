#ifndef PCSIM_GENERATOR_LOADED
#define PCSIM_GENERATOR_LOADED

#include <memory>
#include <TRandom.h>

#include <pcsim/core/configuration.hh>
#include <pcsim/core/histogrammer.hh>


namespace pcsim {

// Base class for all generators. Using a recursive template implementation for
// maximum runtime performace.
//
// The base class owns a shared pointer to the random generator, and handles the
// optional configuration and histogramming options
template <class Derived, class Event>
class generator : public configurable {
public:
  using gen_type = Derived;
  using event_type = Event;
  using histogrammer_type = histogrammer<event_type>;
  using histo_var_type = typename histogrammer_type::histo_var_type;

  generator(const ptree& settings, const string_path& path,
            const std::string& title, std::shared_ptr<TRandom> r)
      : configurable{settings, path}
      , histos_{path, "", title}
      , rng_{std::move(r)} {}

  // add a 1D histo
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const std::string& title, const histo_var_type& var) {
    histos_.add_histo(file, name, title, var);
  }
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const histo_var_type& var) {
    add_histo(std::move(file), name, name, var);
  }
  // add a 2D histo
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const std::string& title, const histo_var_type& var_x,
                 const histo_var_type& var_y) {
    histos_.add_histo(file, name, title, var_x, var_y);
  }
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const histo_var_type& var_x, const histo_var_type& var_y) {
    add_histo(std::move(file), name, name, var_x, var_y);
  }

  template <class... Input> event_type generate(const Input&... input) {
    const auto& event = derived().gen_impl(input...);
    histos_.fill(event);
    return event;
  }
  template <class... Input> event_type operator()(const Input&... input) {
    return generate(input...);
  }

  TRandom& rng() const { return *rng_; }

private:
  gen_type& derived() { return static_cast<gen_type&>(*this); }
  const gen_type& derived() const {
    return static_cast<const gen_type&>(*this);
  }

  std::shared_ptr<TRandom> rng_;
  histogrammer_type histos_;
};

}

#endif
