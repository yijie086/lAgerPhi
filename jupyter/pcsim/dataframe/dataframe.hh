#ifndef JUPYTER_PCSIM_DATAFRAME_LOADED
#define JUPYTER_PCSIM_DATAFRAME_LOADED

#include <ROOT/TDataFrame.hxx>
#include <string>
#include <vector>
#include <TTree.h>

namespace dataframe {
// shortcuts
using ROOT::Experimental::TDataFrame;
using tdf_col_interface_type = decltype((TDataFrame{"", ""}).Define("", ""));

// 'custom dataframe' is a TInterface to a dataframe that also stores the
// original dataframe. This way we can use a constructor that automatically
// attaches several defines to a dataframe, without first having to create the
// dataframe in a seperate call
class custom_dataframe : public tdf_col_interface_type {
public:
  using parent_type = tdf_col_interface_type;
  //using column_names_type = ROOT::Detail::TDF::TDFDetail::TDFColumnNames_t;

  custom_dataframe(std::string_view treeName,
                   const std::vector<std::string>& filenamescoll/*,
                   const column_names_type& default_branches = {}*/)
      : parent_type(nullptr)
      , df_{std::move(treeName), filenamescoll/*, default_branches*/} {}
  custom_dataframe(std::string_view treeName, std::string_view filenameglob/*,
                   const column_names_type& default_branches = {}*/)
      : parent_type(nullptr)
      , df_{std::move(treeName), std::move(filenameglob)/*, default_branches*/} {}
  custom_dataframe(TTree& tree/*, const column_names_type& default_branches = {}*/)
      : parent_type(nullptr), df_{tree/*, default_branches*/} {}

  // get a handle to the bare dataframe
  TDataFrame& bare() { return df_; }
  const TDataFrame& bare() const { return df_; }

private:
  // actual defines, to be defined by the child class
  virtual parent_type defines() = 0;

  // helper function that calls define, and then attaches the result to the
  // parent tdf_col_interface_type
  void attach_defines() {
    auto new_interface = defines();
    static_cast<parent_type&>(*this) = new_interface;
  }

  // actual dataframe
  TDataFrame df_;
};
}

#endif
