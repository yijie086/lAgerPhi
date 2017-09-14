#ifndef JUPYTER_PCSIM_DATAFRAME_LOADED
#define JUPYTER_PCSIM_DATAFRAME_LOADED

#include <ROOT/TDataFrame.hxx>
#include <string>
#include <vector>
#include <TTree.h>

namespace dataframe {
// shortcuts
using ROOT::Experimental::TDataFrame;
using def_interface_type = decltype((TDataFrame{"", ""}).Define("", ""));
using fil_interface_type = decltype((TDataFrame{"", ""}).Filter("", ""));
using histo1D_type = decltype((TDataFrame{"", ""}).Histo1D(""));
using histo2D_type = decltype((TDataFrame{"", ""}).Histo1D(""));

class dataframe_proxy {
public:
  dataframe_proxy(const TDataFrame& df) : df_{df} {}
  dataframe_proxy(TDataFrame&& df) : df_{std::move(df)} {}

  TDataFrame& bare_dataframe() { return df_; }
  const TDataFrame bare_dataframe() const { return df_; }

private:
  TDataFrame df_;
};

class custom_dataframe : public dataframe_proxy, public def_interface_type {
public:
  // scale: constant scale factor to apply to each event
  custom_dataframe(const TDataFrame& df, const double scale = 1.)
      : dataframe_proxy{df}
      , def_interface_type{bare_dataframe().Define("dummy", "index")}
      , scale_{scale} {}
  custom_dataframe(TDataFrame&& df, const double scale = 1.)
      : dataframe_proxy{std::move(df)}
      , def_interface_type{bare_dataframe().Define("dummy", "index + 1")}
      , scale_{scale} {}

protected:
  void init() {
    auto new_interface =
        custom_defines(bare_dataframe()).Define("scale", [=]() {
          return scale_;
        });

    static_cast<def_interface_type&>(*this) = new_interface;
  }
  virtual def_interface_type custom_defines(TDataFrame& df) = 0;

private:
  const double scale_;
};

template <class Interface>
histo1D_type make_histo1D(Interface& df, const TH1F& href, string_view vname,
                          string_view wname = "") {
  auto h2 = href;
  if (wname.size() == 0) {
    return df.Histo1D(std::move(h2), vname);
  }
  return df.Histo1D(std::move(h2), vname, wname);
}
#if 0
inline histo1D_type make_histo1D(def_interface_type& df, const TH1F& href,
                          const std::string& vname,
                          const std::string& wname = "") {
  auto h2 = href;
  return df.Histo1D(std::move(h2), vname, wname);
}
inline histo1D_type make_histo1D(fil_interface_type& df, const TH1F& href,
                          const std::string& vname,
                          const std::string& wname = "") {
  auto h2 = href;
  return df.Histo1D(std::move(h2), vname, wname);
}
inline histo1D_type make_histo1D(TDataFrame& df, const TH1F& href,
                          const std::string& vname,
                          const std::string& wname = "") {
  auto h2 = href;
  return df.Histo1D(std::move(h2), vname, wname);
}
#endif
}

#endif
