#include <TFile.h>
#include <TTree.h>
#include <memory>
#include <pcsim/core/event.hh>
#include <string>
#include <vector>

// =============================================================================
// event_out
//
// Output for the base event record.
//
// Derive from this class for more specialized event records.
// (where the specialized event record derives from the main event class)
//    * Make sure to define your own push(your_event_type) method, and call the
//      parent push(parent_event_type) from within the method.
//    * you are responsible to create the necessary branches for your custom
//      event type, the main event branches are added by this base class
// =============================================================================

namespace pcsim {
class event_out {
public:
  event_out(std::shared_ptr<TFile> f, const std::string& name);
  ~event_out() { tree_->AutoSave(); }

  // add a new event to the event buffer, and flush the buffer to the tree
  void push(const event& e);

  TTree* tree() { return tree_; }

private:
  // clear particle portion of the event buffer
  void clear();
  // add a particle to the event buffer
  void add(const particle& part);
  // add a detected particle to the buffer
  void add_detected(const detected_particle& dp);

  void create_branches();

private:
  // file and tree
  std::shared_ptr<TFile> file_;
  TTree* tree_; // raw pointer because the TFile will have ownership of the tree

  // event data
  int32_t index_{0};
  int32_t evgen_;
  float cross_section_;
  float epsilon_;
  float R_;
  float weight_;
  int32_t process_;
  float s_;
  int16_t beam_index_;
  int16_t target_index_;

  // particle data
  int16_t n_part_ {0};
  std::vector<int32_t> type_;
  std::vector<int16_t> status_;
  std::vector<int16_t> charge_;
  std::vector<float> mass_;
  std::vector<particle::XYZTVector> p_;
  std::vector<particle::XYZTVector> vertex_;
  std::vector<int16_t> parent_first_;
  std::vector<int16_t> parent_second_;
  std::vector<int16_t> daughter_begin_;
  std::vector<int16_t> daughter_end_;
  // detected particle data
  int16_t rc_n_part_ {0};
  std::vector<particle::XYZTVector> rc_p_;
  std::vector<int32_t> rc_type_;
  std::vector<int16_t> rc_parent_;
  std::vector<int16_t> rc_status_;
};
} // namespace pcsim
