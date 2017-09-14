#ifndef JUPYTER_PCSIM_DATAFRAME_LOADED
#define JUPYTER_PCSIM_DATAFRAME_LOADED

#include <ROOT/TDataFrame.hxx>
#include <string>
#include <vector>
#include <TTree.h>

namespace dataframe {
// shortcuts
using ROOT::Experimental::TDataFrame;
using col_interface_type = decltype((TDataFrame{"", ""}).Define("", ""));

class dataframe_proxy {
public:
  dataframe_proxy(const TDataFrame& df) : df_{df} {}
  dataframe_proxy(TDataFrame&& df) : df_{std::move(df)} {}

  TDataFrame& bare_dataframe() { return df_; }
  const TDataFrame bare_dataframe() const { return df_; }

private:
  TDataFrame df_;
};

class custom_dataframe : public dataframe_proxy, public col_interface_type {
public:
  // scale: constant scale factor to apply to each event
  custom_dataframe(const TDataFrame& df, const double scale = 1.)
      : dataframe_proxy{df}
      , col_interface_type{bare_dataframe().Define("dummy", "index")}
      , scale_{scale} {}
  custom_dataframe(TDataFrame&& df, const double scale = 1.)
      : dataframe_proxy{std::move(df)}
      , col_interface_type{bare_dataframe().Define("dummy", "index + 1")}
      , scale_{scale} {}

protected:
  void init() {
    auto new_interface =
        custom_defines(bare_dataframe()).Define("scale", [=]() {
          return scale_;
        });

    static_cast<col_interface_type&>(*this) = new_interface;
  }
  virtual col_interface_type custom_defines(TDataFrame& df) = 0;

private:
  const double scale_;
};

}

#endif
