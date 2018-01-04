#ifndef PCSIM_ROOT_DATAFRAME_HISTOGRAMMER_LOADED
#define PCSIM_ROOT_DATAFRAME_HISTOGRAMMER_LOADED

#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <boost/any.hpp>
#include <boost/filesystem.hpp>
#include <map>
#include <memory>
#include <pcsim/core/interval.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/root/dataframe/dataframe.hh>
#include <string>
#include <vector>

namespace pcsim {
namespace root {
namespace dataframe {

using range_type = interval<double>;
using options_type = std::map<std::string, boost::any>;

namespace histogrammer_impl {
// =============================================================================
// HISTOGRAMMER_IMPL: HISTO_STYLE
// =============================================================================
class histo_style {
public:
  histo_style(int color, int lstyle, int lwidth = 2);
  void apply(histo1D_type& h) const;
  void apply(histo2D_type& h) const;
  void apply(histo3D_type& h) const;

private:
  int lcolor_;
  int mcolor_;
  int lstyle_;
  int mstyle_;
  int lwidth_;
  double msize_;
  int contour_;
};
extern const std::vector<histo_style> STYLE;
// =============================================================================
// HISTOGRAMMER_IMPL: CONFIGURABLE
//
// TODO: confusing, this is not pcsim::configurable!
// =============================================================================
class configurable {
public:
  configurable(const options_type& opts, const options_type& default_opts = {});

  // add settings from the default if not already set
  void set_defaults(const options_type& default_opts);

  // get our options
  const options_type& options() const { return options_; }

protected:
  std::string generic_figure_name() const;
  std::string generic_histo_name() const;

private:
  void convert_char2string();

  options_type options_;
};

// =============================================================================
// HISTOGRAMMER_IMPL: HISTO_PROXY
// =============================================================================
template <class Histo> class histo_proxy : public configurable {
public:
  histo_proxy(const Histo& h, const options_type& opts = {});

  // initialize the histogram
  void init(); 

  // draw this histo
  void draw(const std::string& drawopt) { histo_->DrawClone(drawopt.c_str()); }
  // save this histo
  void write() { histo_->Write(); }

  // set the histo style to style i from the STYLE palette
  void set_style(const size_t i) { STYLE.at(i).apply(histo_); }

  // histo name
  std::string name() { return histo_->GetName(); }

private:
  Histo histo_;
};

// =============================================================================
// HISTOGRAMMER_IMPL: PLOT_PROXY
// =============================================================================
template <class HistoProxy> class plot_proxy : public configurable {
public:
  plot_proxy(const std::vector<HistoProxy>& h, const options_type& opts = {},
             const options_type& default_opts = {});
  
  // draw this plot
  void draw(const std::string& drawopt = "");
  void print(const std::string& drawopt = "");
  void write(const std::string& drawopt = "");

private:
  void init_histos();
  void init_canvas();
  std::vector<HistoProxy> histos_;
  std::shared_ptr<TCanvas> c_;
};
} // namespace histogrammer_impl

// =============================================================================
// HISTOGRAMMER
// =============================================================================
class histogrammer : public histogrammer_impl::configurable {
public:
  using histo1D_type = histogrammer_impl::histo_proxy<dataframe::histo1D_type>;
  using histo2D_type = histogrammer_impl::histo_proxy<dataframe::histo2D_type>;
  using histo3D_type = histogrammer_impl::histo_proxy<dataframe::histo3D_type>;
  using plot1D_type = histogrammer_impl::plot_proxy<histo1D_type>;
  using plot2D_type = histogrammer_impl::plot_proxy<histo2D_type>;
  using plot3D_type = histogrammer_impl::plot_proxy<histo3D_type>;
  using options_type = options_type;

  histogrammer(const options_type& opts = {},
               std::shared_ptr<TFile> file = nullptr,
               const std::string& tfile_dir = "");
  // print and save to file
  ~histogrammer();

  // 1D plots
  void add(const histo1D_type& histo, const options_type& plot_opts = {});
  void add(const std::vector<histo1D_type>& histos,
           const options_type& plot_opts = {});
  // 2D plots
  void add(const histo2D_type& histo, options_type plot_opts = {});
  // 3D "plots"
  void add(const histo3D_type& histo, options_type plot_opts = {});

  void print();
  void write();

private:
  std::shared_ptr<TFile> ofile_;
  const std::string tfile_dir_;
  std::vector<plot1D_type> plots1D_;
  std::vector<plot2D_type> plots2D_;
  std::vector<plot3D_type> plots3D_;
};

} // namespace dataframe
} // namespace root
} // namespace pcsim

// =============================================================================
// IMPLEMENTATION FOR: HISTOGRAMMER_IMPL: HISTO_PROXY
// =============================================================================
namespace pcsim {
namespace root {
namespace dataframe {
namespace histogrammer_impl {
template <class Histo>
histo_proxy<Histo>::histo_proxy(const Histo& h, const options_type& opts)
    : configurable{opts}, histo_{h} {
  // ensure histo has a name
  if (std::string("") == name()) {
    histo_->SetName(generic_histo_name().c_str());
  }
}

// initialize the histogram
template <class Histo> void histo_proxy<Histo>::init() {
  for (const auto& opt : options()) {
    const options_type::key_type& key = opt.first;
    const options_type::mapped_type& value = opt.second;
    if (key == "xrange") {
      auto range = boost::any_cast<range_type>(value);
      histo_->GetXaxis()->SetRangeUser(range.min, range.max);
    } else if (key == "yrange") {
      auto range = boost::any_cast<range_type>(value);
      histo_->GetYaxis()->SetRangeUser(range.min, range.max);
    } else if (key == "zrange") {
      auto range = boost::any_cast<range_type>(value);
      histo_->GetZaxis()->SetRangeUser(range.min, range.max);
    } else if (key == "fit") {
      auto func = boost::any_cast<const char*>(value);
      histo_->Fit(func);
    } else if (key == "color") {
      int color = boost::any_cast<int>(value);
      histo_->SetLineColor(color);
      histo_->SetMarkerColor(color);
    } else if (key == "lcolor") {
      int color = boost::any_cast<int>(value);
      histo_->SetLineColor(color);
    } else if (key == "mcolor") {
      int color = boost::any_cast<int>(value);
      histo_->SetMarkerColor(color);
    } else if (key == "lstyle") {
      int style = boost::any_cast<int>(value);
      histo_->SetLineStyle(style);
    } else if (key == "mstyle") {
      int style = boost::any_cast<int>(value);
      histo_->SetMarkerStyle(style);
    } else if (key == "lwidth") {
      int width = boost::any_cast<int>(value);
      histo_->SetLineWidth(width);
    } else if (key == "msize") {
      double size = boost::any_cast<double>(value);
      histo_->SetMarkerStyle(size);
    } else {
      LOG_WARNING("histogrammer", "Unknown histo option: " + key);
    }
  }
}
} // namespace histogrammer_impl
} // namespace dataframe
} // namespace root
} // namespace pcsim

// =============================================================================
// HISTOGRAMMER_IMPL: PLOT_PROXY
// =============================================================================
namespace pcsim {
namespace root {
namespace dataframe {
namespace histogrammer_impl {
template <class HistoProxy>
plot_proxy<HistoProxy>::plot_proxy(const std::vector<HistoProxy>& h,
                                   const options_type& opts,
                                   const options_type& default_opts)
    : configurable{opts, default_opts}, histos_{h} {
  // ensure dir is always set
  if (options().count("dir") == 0) {
    set_defaults({{"dir", "./"}});
  }
  // ensure name is always set, defaults to generic figure name, or histo in
  // case of 1 histo
  if (options().count("name") == 0) {
    if (histos_.size() == 1) {
      set_defaults({{"name", histos_[0].name()}});
    } else {
      set_defaults({{"name", generic_figure_name()}});
    }
  }
  options_type histo_defaults = {};
  if (options().count("histo")) {
    histo_defaults = boost::any_cast<options_type>(options().at("histo"));
  }
  for (size_t i = 0; i < histos_.size(); ++i) {
    histos_[i].set_style(i);
    histos_[i].set_defaults(histo_defaults);
  }
}

// draw this plot
template <class HistoProxy>
void plot_proxy<HistoProxy>::draw(const std::string& drawopt) {
  if (c_) {
    // already drawn
    return;
  }
  if (histos_.empty()) {
    LOG_WARNING("histogrammer", "Attempting to draw empty figure");
    return;
  }
  init_canvas();
  init_histos();
  std::string extra_drawopt = "";
  for (auto& h : histos_) {
    h.draw(drawopt + extra_drawopt);
    extra_drawopt = "same";
  }
  if (histos_.size() > 1) {
    c_->BuildLegend();
  }
}
template <class HistoProxy>
void plot_proxy<HistoProxy>::print(const std::string& drawopt) {
  draw(drawopt);
  std::string dir{boost::any_cast<std::string>(options().at("dir"))};
  if (dir.back() != '/') {
    dir += '/';
  }
  // create directory if needed
  if (!boost::filesystem::is_directory(dir)) {
    boost::filesystem::create_directories(dir);
  }
  const auto c_name = boost::any_cast<std::string>(options().at("name"));
  c_->Print((dir + c_name + ".pdf").c_str());
}
template <class HistoProxy>
void plot_proxy<HistoProxy>::write(const std::string& drawopt) {
  if (c_) {
    c_->Write();
    draw(drawopt);
  }
  for (auto& h : histos_) {
    h.write();
  }
}

template <class HistoProxy> void plot_proxy<HistoProxy>::init_histos() {
  for (auto& h : histos_) {
    h.init();
  }
}
template <class HistoProxy> void plot_proxy<HistoProxy>::init_canvas() {
  const auto c_name =
      ("plot_" + boost::any_cast<std::string>(options().at("name"))).c_str();
  c_ = std::make_shared<TCanvas>(c_name, c_name, 600, 600);
  for (const auto& opt : options()) {
    const options_type::key_type& key = opt.first;
    const options_type::mapped_type& value = opt.second;
    if (key == "logx") {
      c_->SetLogx(boost::any_cast<bool>(value));
    } else if (key == "logy") {
      c_->SetLogy(boost::any_cast<bool>(value));
    } else if (key == "logz") {
      c_->SetLogz(boost::any_cast<bool>(value));
    } else if (key == "grid") {
      c_->SetGridx(boost::any_cast<bool>(value));
      c_->SetGridy(boost::any_cast<bool>(value));
    } else if (key == "gridx") {
      c_->SetGridx(boost::any_cast<bool>(value));
    } else if (key == "gridy") {
      c_->SetGridy(boost::any_cast<bool>(value));
    } else if (key == "rightmargin") {
      c_->SetRightMargin(boost::any_cast<double>(value));
    } else if (key == "leftmargin") {
      c_->SetLeftMargin(boost::any_cast<double>(value));
    } else if (key == "topmargin") {
      c_->SetTopMargin(boost::any_cast<double>(value));
    } else if (key == "bottommargin") {
      c_->SetBottomMargin(boost::any_cast<double>(value));
    } else if (key == "optstat") {
      gStyle->SetOptStat(boost::any_cast<int>(value));
    } else if (key == "optfit") {
      gStyle->SetOptStat(boost::any_cast<int>(value));
    } else if (key == "title") {
      gStyle->SetOptTitle(boost::any_cast<bool>(value));
    } else if (key == "histo") {
      // magical histo options passed to all histograms in plot, do nothing
      // here.
    } else if (key == "dir") {
      // magical keyword that stores the (optional) directory for the plots
    } else if (key == "name") {
      // magical keywoard that stores the plot (canvas) name
    } else {
      LOG_WARNING("histogrammer", "Unknown plot option: " + key);
    }
  }
}

} // namespace histogrammer_impl
} // namespace dataframe
} // namespace root
} // namespace pcsim




#endif
