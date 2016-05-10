#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <cassert>
#include <string>
#include <vector>

// generate theta vs p plots for different spectrometer settings

const std::string PATH = "../data";
const std::vector<std::string> SET = {"set1", "set2", "set3", "set4"};

TTree* get_tree(TFile* f) {
  TTree* t = (TTree*)f->Get("jpsi_event");
  assert(t);
  return t;
}

void plot_theta_lepton(const std::string& path, const std::string& proc,
                  const std::string& set) {

  // input files
  TFile* f2arm = new TFile(
      (path + "/pcsim." + proc + ".2arm." + set + ".run00001-1000.root")
          .c_str(),
      "read");
  TFile* f1arm = new TFile((path + "/pcsim." + proc + ".2arm." + set +
                            "-SHMSonly.run00001-10000.root")
                               .c_str(),
                           "read");
  TFile* f4pi = new TFile(
      (path + "/pcsim." + proc + ".4pi." + "run00001-1000000.root").c_str(),
      "read");
  assert(f2arm && f1arm && f4pi);

  // input trees
  TTree* t2arm = get_tree(f2arm);
  TTree* t1arm = get_tree(f1arm);
  TTree* t4pi = get_tree(f4pi);

  TCanvas* c = new TCanvas;
  gStyle->SetOptStat(0);

  // histos
  t2arm->Draw("positron.Theta() * 180. / TMath::Pi() : positron.Vect().Mag() "
              ">> h2arm_pos(50,0,11,50,0,80)",
              "", "");
  TH2F* h2arm_pos = (TH2F*) gROOT->FindObject("h2arm_pos");
  h2arm_pos->SetMarkerColor(kGreen+1);
  h2arm_pos->SetMarkerStyle(7);
  h2arm_pos->SetFillColor(kGreen+1);
  t2arm->Draw("electron.Theta() * 180. / TMath::Pi() : electron.Vect().Mag() "
              ">> h2arm_el(50,0,11,50,0,80)",
              "", "");
  TH2F* h2arm_el = (TH2F*) gROOT->FindObject("h2arm_el");
  h2arm_el->SetMarkerColor(kAzure+1);
  h2arm_el->SetMarkerStyle(7);
  h2arm_el->SetFillColor(kAzure+1);

  t1arm->Draw("positron.Theta() * 180. / TMath::Pi() : positron.Vect().Mag() "
              ">> h1arm_pos(50,0,11,50,0,80)",
              "", "");
  TH2F* h1arm_pos = (TH2F*) gROOT->FindObject("h1arm_pos");
  h1arm_pos->SetMarkerStyle(7);
  h1arm_pos->SetMarkerColor(kRed+1);
  h1arm_pos->SetFillColor(kRed+1);
  t1arm->Draw("electron.Theta() * 180. / TMath::Pi() : electron.Vect().Mag() "
              ">> h1arm_el(50,0,11,50,0,80)",
              "", "");
  TH2F* h1arm_el = (TH2F*) gROOT->FindObject("h1arm_el");
  h1arm_el->SetMarkerStyle(7);
  h1arm_el->SetMarkerColor(kMagenta+2);
  h1arm_el->SetFillColor(kMagenta+2);

  t4pi->Draw("electron.Theta() * 180. / TMath::Pi() : electron.Vect().Mag() "
             ">> h4pi_el(50,0,11,50,0,80)",
             "", "colz");
  TH2F* h4pi_el = (TH2F*) gROOT->FindObject("h4pi_el");
  h4pi_el->GetXaxis()->SetTitle("P [GeV]");
  h4pi_el->GetYaxis()->SetTitle("#theta [deg.]");
  h4pi_el->SetTitle(
      ("Lepton #theta vs P for " + set +
       (proc == "jpsi" ? " (J/#Psi t-channel)" : " (P_{c} s-channel)"))
          .c_str());

  h1arm_el->Draw("same");
  h1arm_pos->Draw("same");
  h2arm_el->Draw("same");
  h2arm_pos->Draw("same");

  TLegend* leg = new TLegend(0.4, 0.60, 0.85, 0.85);
  leg->AddEntry(h1arm_pos, "SHMS positrons", "f");
  leg->AddEntry(h1arm_el, "electrons with SHMS positrons", "f");
  leg->AddEntry(h2arm_pos, "positrons after full HMS+SHMS cuts", "f");
  leg->AddEntry(h2arm_el, "electrons after full HMS+SHMS cuts", "f");
  leg->Draw();

  c->Print(("pcsim." + set + "." + proc + ".lepton_theta.pdf").c_str());
}

void theta_lepton() {
  for (const auto& set : SET) {
//    plot_theta_lepton(PATH, "jpsi", set);
 //   plot_theta_lepton(PATH, "pc", set);
  }
  plot_theta_lepton(PATH, "jpsi", "set1");
}
