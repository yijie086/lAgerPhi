#include "dataframe.hh"

namespace dataframe {
histo1D_type make_histo1D(def_interface_type& df, const TH1F& href,
                          const std::string& vname,
                          const std::string& wname = "") {
  auto h2 = href;
  return df.Histo1D(std::move(h2), vname, wname);
}
histo1D_type make_histo1D(fil_interface_type& df, const TH1F& href,
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
}
