#include "dataframe.hh"
#include "vm.hh"

namespace dataframe {
class pcsim_vm : public custom_dataframe {
public:
  pcsim_vm(const std::string& fname) : custom_dataframe{fname, VM_TREE_NAME} {}

  virtual custom_dataframe::parent_type defines() {
    return bare().Define("Q4", "Q2*Q2");
  }
};
} // namespace dataframe
