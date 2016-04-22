#ifndef PYTHIA6M_UTIL_PROGRESS_METER_LOADED
#define PYTHIA6M_UTIL_PROGRESS_METER_LOADED

#include <cmath>
#include <cstdint>
#include <iostream>

// =============================================================================
// A simple console progress meter
// =============================================================================
namespace pythia6m {
class progress_meter {
public:
  constexpr static const size_t PRECISION = 1000.; // 0.1 percent
  explicit progress_meter(size_t max, size_t start_index = 0,
                          size_t precision = PRECISION)
      : max_{max}, index_{start_index}, precision_{precision} {
    std::cout << "\nProcessing " << max_ << " events..." << std::endl;
    update();
  }
  void update(size_t i = 1) {
    index_ += i;
    if (index_ > max_) {
      index_ = max_;
    }
    // update when needed
    if (!index_ || !(index_ % (max_ / precision_))) {
      double cnt = index_ * precision_ / max_;
      cnt /= (precision_ / 100.);
      char msg[10];
      sprintf(msg, "  %3.2f%%\r", cnt);
      std::cout << msg << std::flush;
    }
  }
  ~progress_meter() { std::cout << "      Done!" << std::endl; }

private:
  const size_t max_;
  size_t index_;
  const size_t precision_;
};
} // ns pythia6m

#endif
