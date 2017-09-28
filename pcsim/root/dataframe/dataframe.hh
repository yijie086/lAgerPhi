#ifndef PCSIM_ROOT_DATAFRAME_DATAFRAME_LOADED
#define PCSIM_ROOT_DATAFRAME_DATAFRAME_LOADED

#include <ROOT/TDataFrame.hxx>
#include <string>
#include <vector>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

namespace pcsim {
namespace root {
namespace dataframe {

// shortcuts
using ROOT::Experimental::TDataFrame;
using def_interface_type = decltype((TDataFrame{"", ""}).Define("", ""));
using fil_interface_type = decltype((TDataFrame{"", ""}).Filter("", ""));
using histo1D_type = decltype((TDataFrame{"", ""}).Histo1D(""));
using histo2D_type = decltype((TDataFrame{"", ""}).Histo2D(TH2F()));

// utility class
class dataframe_proxy {
public:
  dataframe_proxy(const TDataFrame& df) : df_{df} {}
  dataframe_proxy(TDataFrame&& df) : df_{std::move(df)} {}

  TDataFrame& bare_dataframe() { return df_; }
  const TDataFrame bare_dataframe() const { return df_; }

private:
  TDataFrame df_;
};

// a custom dataframe that automizes some of the defines/filters that we always
// want to apply to a dataset to present a more uniform interface
class custom_dataframe : public dataframe_proxy, public def_interface_type {
public:
  // scale: constant scale factor to apply to each event
  custom_dataframe(const TDataFrame& df, const double scale = 1.)
      : dataframe_proxy{df}
      , def_interface_type{
            bare_dataframe().Define("scale", [=]() { return scale; })} {}
  custom_dataframe(TDataFrame&& df, const double scale = 1.)
      : dataframe_proxy{std::move(df)}
      , def_interface_type{
            bare_dataframe().Define("scale", [=]() { return scale; })} {}

  def_interface_type& interface() {
    return static_cast<def_interface_type&>(*this);
  }
  const def_interface_type& interface() const {
    return static_cast<const def_interface_type&>(*this);
  }
};


// utility function to make a histogram starting from a template histogram. This
// is only neccesary for PyROOT, as the current python bindings do not work
// (9/26/2017)
template <class Interface>
histo1D_type make_histo1D(Interface& df, const TH1F& href,
                          std::string_view vname, std::string_view wname = "") {
  auto h2 = href;
  if (wname.size() == 0) {
    return df.Histo1D(std::move(h2), vname);
  }
  return df.Histo1D(std::move(h2), vname, wname);
}
template <class Interface>
histo2D_type make_histo2D(Interface& df, const TH2F& href,
                          std::string_view v1name, std::string_view v2name,
                          std::string_view wname = "") {
  auto h2 = href;
  if (wname.size() == 0) {
    return df.Histo2D(std::move(h2), v1name, v2name);
  }
  return df.Histo2D(std::move(h2), v1name, v2name, wname);
}

} // namespace dataframe
} // namespace root
} // namespace pcsim

#endif
