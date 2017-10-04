#include "histogrammer.hh"

namespace pcsim {
namespace root {
namespace dataframe {

namespace histogrammer_impl {
// =============================================================================
// HISTOGRAMMER_IMPL: HISTO_STYLE
// =============================================================================
histo_style::histo_style(int color, int lstyle, int lwidth)
    : lcolor_{color}
    , mcolor_{color}
    , lstyle_{lstyle}
    , mstyle_{0}
    , lwidth_{lwidth}
    , msize_{1.}
    , contour_{1000} {}
void histo_style::apply(histo1D_type& h) const {
  h->SetLineColor(lcolor_);
  h->SetMarkerColor(mcolor_);
  h->SetLineStyle(lstyle_);
  h->SetMarkerStyle(mstyle_);
  h->SetLineWidth(lwidth_);
  h->SetMarkerSize(msize_);
}
void histo_style::apply(histo2D_type& h) const {
  h->SetContour(contour_); // for 2D colz histos
}
const std::vector<histo_style> STYLE = {
    {kBlack, 0},     {kBlue + 1, 1},   {kGreen + 2, 2}, {kRed + 1, 3},
    {kViolet, 4},    {kOrange + 1, 5}, {kBlack, 2},     {kBlue + 1, 3},
    {kGreen + 2, 4}, {kRed + 1, 5},    {kViolet, 0},    {kOrange + 1, 1},
    {kBlack, 4},     {kBlue + 1, 5},   {kGreen + 2, 0}, {kRed + 1, 1},
    {kViolet, 2},    {kOrange + 1, 3}};

// =============================================================================
// HISTOGRAMMER_IMPL: CONFIGURABLE
// =============================================================================
configurable::configurable(const options_type& opts,
                           const options_type& default_opts)
    : options_{opts} {
  set_defaults(default_opts);
}
// add settings from the default if not already set
void configurable::set_defaults(const options_type& default_opts) {
  for (const auto& setting : default_opts) {
    if (!options_.count(setting.first)) {
      options_[setting.first] = setting.second;
    }
  }
}
} // namespace histogrammer_impl

// =============================================================================
// HISTOGRAMMER
// =============================================================================
histogrammer::histogrammer(const options_type& opts)
    : histogrammer_impl::configurable{opts} {}
void histogrammer::add(const histogrammer::histo1D_type& histo,
                       const options_type& plot_opts) {
  add(std::vector<histo1D_type>({histo}), plot_opts);
}
void histogrammer::add(const std::vector<histogrammer::histo1D_type>& histos,
                       const options_type& plot_opts) {
  plots1D_.push_back({histos, plot_opts, options()});
}
void histogrammer::add(const histogrammer::histo2D_type& histo,
                       options_type plot_opts) {
  // fix the rightmargin for standard colz figures
  if (!plot_opts.count("rightmargin")) {
    plot_opts["rightmargin"] = .15;
  }
  plots2D_.push_back(
      {std::vector<histo2D_type>({histo}), plot_opts, options()});
}
void histogrammer::print() {
  for (auto& plot : plots1D_) {
    plot.print();
  }
  for (auto& plot : plots2D_) {
    plot.print();
  }
}

} // namespace dataframe
} // namespace root
} // namespace pcsim
