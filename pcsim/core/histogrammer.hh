#ifndef ANALYZER_CORE_HISTOGRAMMER_LOADED
#define ANALYZER_CORE_HISTOGRAMMER_LOADED

// an easy to define 1D and 2D histo that will automatically get its own
// variables from Data using a custom getter (f(Data...){return Double_t;})

#include <functional>
#include <memory>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <analyzer/core/configuration.hh>
#include <analyzer/core/interval.hh>

// =============================================================================
// helper class for more elegant histogram constructors
//
// members:
//  * name: this variable name
//  * getter: function that takes computes the variable value from (Data...)
//  * n_bins: number of histogram bins
//  * range: the [min, max] range of the histogram
// =============================================================================
namespace analyzer {
template <class... Data> struct histogram_variable {
  using getter_type = std::function<Double_t(Data...)>;

  const std::string name;
  const getter_type getter;
  const Int_t n_bins;
  const interval<Double_t> range;

  histogram_variable(const std::string& name, const getter_type& getter,
                     const Int_t n_bins, const interval<Double_t> range)
      : name{name}, getter{getter}, n_bins{n_bins}, range{range} {}
};
} // ns analyzer

// =============================================================================
// Histogrammer class to auto-generate 1D and 2D histos from data
//
// members:
//  * add_histo(shared_ptr<TFile>, name, var_x[, var_y])
//  * add_histo(shared_ptr<TFile>, name, title, var_x[, var_y])
//        ==> add 1D and 2D histos to the histogrammer
//  * fill(Data..., weight)
//        ==> fill all histograms using Data... as input
// =============================================================================
namespace analyzer {
namespace histogrammer_impl {
// utility classes needed by histogrammer
template <class... Data> class histo1D;
template <class... Data> class histo2D;
} // ns histogrammer_impl

template <class... Data> class histogrammer {
public:
  using histo1D_type = histogrammer_impl::histo1D<Data...>;
  using histo2D_type = histogrammer_impl::histo2D<Data...>;
  using histo_var_type = histogram_variable<Data...>;

  histogrammer(const string_path& path, std::string name, std::string title)
      : path_{path.str()}, name_{std::move(name)}, title_{std::move(title)} {}
  histogrammer(const std::string& name) : histogrammer("", name, name) {}
  histogrammer(const histogrammer&) = delete;
  histogrammer& operator=(const histogrammer&) = delete;

  // add 1D histos
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const std::string& title, const histo_var_type& var) {
    histos_1D_.push_back(
        std::move(histo1D_type{std::move(file), path_, format_name(name_, name),
                               format_title(title_, title), var}));
  }
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const histo_var_type& var) {
    add_histo(std::move(file), name, name, var);
  }
  // add 2D histos
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const std::string& title, const histo_var_type& var_x,
                 const histo_var_type& var_y) {
    histos_2D_.push_back(
        std::move(histo2D_type{std::move(file), path_, format_name(name_, name),
                               format_title(title_, title), var_x, var_y}));
  }
  void add_histo(std::shared_ptr<TFile> file, const std::string& name,
                 const histo_var_type& var_x, const histo_var_type& var_y) {
    add_histo(std::move(file), name, name, var_x, var_y);
  }

  // fill the histos
  void fill(Data... d, Double_t weight = 1.) {
    for (auto& h : histos_1D_) {
      h.fill(d..., weight);
    }
    for (auto& h : histos_2D_) {
      h.fill(d..., weight);
    }
  }

private:
  const std::string path_;
  const std::string name_;
  const std::string title_;
  std::vector<histo1D_type> histos_1D_;
  std::vector<histo2D_type> histos_2D_;
};

} // ns analyzer

// =============================================================================
// Implementation: histo1D and histo2D
// =============================================================================
namespace analyzer {
namespace histogrammer_impl {
template <class... Data> class histo1D
{
public:
  using var_type = histogram_variable<Data...>;
  using getter_type = typename var_type::getter_type;

  histo1D(std::shared_ptr<TFile> file, const std::string& path,
          const std::string& name, const std::string& title,
          const var_type& var)
      : file_{std::move(file)}
      , path_{path}
      , getter_{var.getter}
      , histo_{nullptr} {
    if (file_) {
      file_->cd();
      if (path.size() && !file_->GetDirectory(path.c_str())) {
        file_->mkdir(path.c_str());
      }
      file_->cd(path.c_str());
      histo_ = new TH1D{name.c_str(), title.c_str(), var.n_bins, var.range.min,
                        var.range.max};
      histo_->GetXaxis()->SetTitle(var.name.c_str());
      histo_->GetYaxis()->SetTitle("#");
      file_->cd();
    }
  }
  histo1D(std::shared_ptr<TFile> file, const std::string& name,
          const var_type& var)
      : histo1D(std::move(file), "", name, name, var) {}
  // don't allow regular copying/assigment to avoid unexpected destructor
  // actions
  histo1D(const histo1D&) = delete;
  histo1D& operator=(const histo1D&) = delete;
  // move constructor to allow for handling histo rvalues without spurious
  // destructor invocations
  histo1D(histo1D&& rhs)
      : file_{std::move(rhs.file_)}
      , path_{rhs.path_}
      , getter_{rhs.getter_}
      , histo_{rhs.histo_} {
    rhs.file_.reset();
    rhs.histo_ = nullptr;
  }

  ~histo1D() {
    if (file_ && histo_) {
      if (path_.size()) {
        file_->cd(path_.c_str());
      }
      histo_->Write();
      file_->cd();
    }
  }

  Int_t fill(Data... d, Double_t weight = 1.) {
    if (file_ && histo_) {
      return histo_->Fill(getter_(d...), weight);
    }
    return -1;
  }

  TH1D* histo() { return histo_; }

private:
  std::shared_ptr<TFile> file_;
  const std::string path_;
  const getter_type getter_;
  TH1D* histo_; // resource managed by file_
};
// 2D histo that uses a two custom 'getter's to extract the values from a Reader
template <class... Data> class histo2D {
public:
  using var_type = histogram_variable<Data...>;
  using getter_type = typename var_type::getter_type;

  histo2D(std::shared_ptr<TFile> file, const std::string& path,
          const std::string& name, const std::string& title,
          const var_type& var_x, const var_type& var_y)
      : file_{std::move(file)}
      , path_{path}
      , getter_x_{var_x.getter}
      , getter_y_{var_y.getter}
      , histo_{nullptr} {
    if (file_) {
      file_->cd();
      if (path.size() && !file_->GetDirectory(path.c_str())) {
        file_->mkdir(path.c_str());
      }
      file_->cd(path.c_str());
      histo_ = new TH2D{name.c_str(),    title.c_str(),   var_x.n_bins,
                        var_x.range.min, var_x.range.max, var_y.n_bins,
                        var_y.range.min, var_y.range.max};
      // set the hist to use the nicer colz drawing option
      histo_->SetOption("colz");
      histo_->GetXaxis()->SetTitle(var_x.name.c_str());
      histo_->GetYaxis()->SetTitle(var_y.name.c_str());
    }
  }
  histo2D(std::shared_ptr<TFile> file, const std::string& name,
          const var_type& var_x, const var_type& var_y)
      : histo2D(std::move(file), "", name, name, var_x, var_y) {}
  // don't allow regular copying/assigment to avoid unexpected destructor
  // actions
  histo2D(const histo2D&) = delete;
  histo2D& operator=(const histo2D&) = delete;
  // move constructor to allow for handling histo rvalues without spurious
  // destructor invocations
  histo2D(histo2D&& rhs)
      : file_{std::move(rhs.file_)}
      , path_{rhs.path_}
      , getter_x_{rhs.getter_x_}
      , getter_y_{rhs.getter_y_}
      , histo_{rhs.histo_} {
    rhs.file_.reset();
    rhs.histo_ = nullptr;
  }

  ~histo2D() {
    if (file_ && histo_) {
      if (path_.size()) {
        file_->cd(path_.c_str());
      }
      histo_->Write();
      file_->cd();
    }
  }

  Int_t fill(Data... d, Double_t weight = 1.) {
    if (file_ && histo_) {
      return histo_->Fill(getter_x_(d...), getter_y_(d...), weight);
    }
    return -1;
  }

  TH2D* histo() { return histo_; }

private:
  std::shared_ptr<TFile> file_;
  const std::string& path_;
  const getter_type getter_x_;
  const getter_type getter_y_;
  TH2D* histo_; // resource managed by file_
};
} // ns histogrammer_impl
} // ns analyzer

#endif
