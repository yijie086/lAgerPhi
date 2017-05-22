#ifndef PCSIM_CORE_PROGRESS_METER_LOADED
#define PCSIM_CORE_PROGRESS_METER_LOADED

#include <cmath>
#include <cstdint>
#include <iostream>

// =============================================================================
// A simple console progress meter
// =============================================================================
namespace pcsim {
class progress_meter {
public:
  constexpr static const size_t PRECISION = 10000.; // 0.01 percent
  explicit progress_meter(size_t max, size_t start_index = 0,
                          size_t precision = PRECISION)
      : max_{max}
      , index_{start_index}
      , precision_{precision > max ? max : precision} {
    std::cerr << "\nProcessing... " << std::endl;
    update();
  }
  void update(size_t i) {
    index_ = i;
    if (index_ > max_) {
      index_ = max_;
    }
    // update when needed
    // commented out because of FPE, need to fix this TODO
//    if (!index_ || !(index_ % (max_ / precision_) || !(index_ % 1000))) {
    if (!(index_ % 100)) {
      double cnt = index_ * precision_ / max_;
      cnt /= (precision_ / 100.);
      char msg[15];
      sprintf(msg, "  %3.2f%%\r", cnt);
      std::cerr << msg << std::flush;
    }
  }
  void update(size_t i, const size_t max) {
    max_ = max;
    update(i);
  }
  void update() {
    index_ += 1;
    update(index_);
  }
  ~progress_meter() { std::cerr << "      Done!" << std::endl; }

private:
  size_t max_;
  size_t index_;
  const size_t precision_;
};
} // ns pcsim

#endif
