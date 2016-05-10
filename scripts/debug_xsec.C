#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSpline.h>
#include <TTree.h>
#include <cmath>
#include <iostream>

const double M_PROTON = 0.938272;                       // proton mass
const double M2_PROTON = M_PROTON * M_PROTON;           // proton mass
const double ONE_OVER_2M_PROTON = 1. / (2. * M_PROTON); // proton
const double M_JPSI = 3.09692;                          // J/Psi pole mass
const double M2_JPSI = M_JPSI * M_JPSI;                 // J/Psi

TF1* tf1;
TFile* f;
TTree* t;
TCanvas* c;

void testbs() {
  // BS spectrum
  double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  double yv[] = {1.15277, 1.04815, 0.96490, 0.90130, 0.85591, 0.82716,
                 0.81289, 0.80929, 0.80921, 0.80908, 0.80873, 0.80789,
                 0.80616, 0.80289, 0.79691, 0.78565, 0.76136, 0.64338};
  const TSpline3 brems10{"brems10", xv, yv, 18};

  auto fbrems = [=](double k) { return brems10.Eval(k / 11) * 0.1 / k; };
  // Pc spectrum
  const double mean = 4.45;
  const double sigma = .039/2.;
  const double ampl = 1.e4;
  auto fpc = [=](double W) { return ampl * TMath::Gaus(W, mean, sigma); };
  // Transformations and jacobians
  auto s2W = [](double s) { return sqrt(s); };
  auto s2k = [](double s) { return (s - M_PROTON) * ONE_OVER_2M_PROTON; };
  auto js2W = [](double s) { return 1 / (2 * sqrt(s)); };
  auto js2k = [](double) { return ONE_OVER_2M_PROTON; };
  // s range from photon energy
  TLorentzVector gmin{0, 0, 8, 8};
  TLorentzVector gmax{0, 0, 11, 11};
  TLorentzVector target{0, 0, 0, M_PROTON};
  const double smin = (gmin+target).M2();
  const double smax = (gmax+target).M2();
  // Functions and integrals
  auto fpc_s = [=](double* s, double* ) {
    return fpc(s2W(s[0])) * js2W(s[0]);
  };
  auto fbrems_s = [=](double* s, double* ) {
    return fbrems(s2k(s[0])) * js2k(s[0]);
  };
  auto fconv_s = [=](double* s, double* p) {
    return fpc_s(s, p) * fbrems_s(s, p);
  };

  TF1 tfpc("tfpc", fpc_s, smin, smax, 0);
  TF1 tfbrems("tfbrems", fbrems_s, smin, smax, 0);
  TF1 tfconv("tfconv", fconv_s, smin, smax, 0);
  std::cout << "Pc: " << tfpc.Integral(smin, smax) << std::endl;
  std::cout << "BS: " << tfbrems.Integral(smin, smax) << std::endl;
  std::cout << "conv: " << tfconv.Integral(smin, smax) << std::endl;
}

void testjpsi() {
  // transformations
  auto s2k = [](double s) { return (s - M2_PROTON) * ONE_OVER_2M_PROTON; };
  auto k2s = [](double k) { return (M2_PROTON + 2 * M_PROTON * k); };
  const double Ethres = (M2_JPSI + 2 * M_JPSI * M_PROTON) * ONE_OVER_2M_PROTON;
  // t limits
  auto tlim = [=](double s) {
    double Ebeam = s2k(s);
    TLorentzVector beam = {0, 0, Ebeam, Ebeam};
    TLorentzVector target = {0, 0, 0, M_PROTON};
    auto cm = beam + target;
    auto beta = cm.BoostVector();
    beam.Boost(-beta);
    target.Boost(-beta);
    cm.Boost(-beta);
    const double Etot = beam.E() + target.E();
    const double Ep = target.E();
    const double Pp = target.Vect().Mag();
    const double Epp = (Etot * Etot - M2_JPSI + M2_PROTON) / (2 * Etot);
    const double Ppp = std::sqrt(Epp * Epp - M2_PROTON);
    return std::pair<double, double>{
        2 * M2_PROTON - 2 * Ep * Epp - 2 * Pp * Ppp,
        2 * M2_PROTON - 2 * Ep * Epp + 2 * Pp * Ppp};
  };
  auto fjpsi = [=](double s, double t) {
    if (s2k(s) <= Ethres) {
      return 0.;
    }
    const auto trange = tlim(s);
    if (t < trange.first || t > trange.second) {
      return 0.;
    }
    const double x = (2.0 * M_PROTON * M_JPSI + M2_JPSI) / (s - M2_PROTON);
    const double v = 1 / (16. * TMath::Pi());
    const double ep = exp(t * 1.13);
    const double A_2g = 6.499e3 * v * (1 - x) * (1 - x) / M2_JPSI;
    return ep * A_2g;
  };
  // s range from photon energy
  const double Ebeam = 11;
  TLorentzVector gmin{0, 0, 8., 8.};
  TLorentzVector gthres{0, 0, Ethres, Ethres};
  TLorentzVector gmax{0, 0, Ebeam, Ebeam};
  TLorentzVector target{0, 0, 0, M_PROTON};
  const double smin = (gmin+target).M2();
  const double smax = (gmax+target).M2();
  const double sthres = (gthres + target).M2();
  std::cout << smin << " " << smax << std::endl;
  // get t range for integration
  const double tmin = tlim(smax).first;
  const double tmax = tlim(smax).second;
  // functions and integrals
  auto fjpsi2 = [=](double* x, double*) { return fjpsi(x[0], x[1]); };
  auto fjpsi_it = [=](double* ss, double*) {
    const double s = ss[0];
    if (s2k(s) <= Ethres) {
      return 0.;
    }
    const auto trange = tlim(s);
    auto fjpsi_t = [=](double* t, double*) { return fjpsi(s, t[0]); };
    TF1 tmp("tmp", fjpsi_t, trange.first, trange.second, 0);
    return tmp.Integral(trange.first, trange.second);
  };
  TF1* tfjpsi_it = new TF1("jpsi_it", fjpsi_it, sthres, smax, 0);
  TF2 tfjpsi("jpsi", fjpsi2, smin, smax, tmin, tmax, 0);
  std::cout << "jpsi: " << tfjpsi.Integral(sthres, smax, tmin, tmax, .001)
            << std::endl;
  tfjpsi_it->SetNpx(10000);
  // persistence/drawing
  f = new TFile{"dbg.run000003-100000.root", "read"};
  t = (TTree*)f->Get("jpsi_event");
  tf1 = tfjpsi_it;
  c = new TCanvas();
  //c->Divide(1,2);
  auto p = c->cd();
  p->SetGridx();
  p->SetGridy();
  p->SetLogy();
  tf1->Draw("");
  //p = c->cd(2);
  //kp->SetGridx();
  //p->SetGridy();
  //p->SetLogy();
  //t->Draw("W>>h(50,4.04498,4.63922)", "xsec_gen/100000*50/.59424", "");
  c->Print("tjpsi_it.pdf");
}

void tlimtest() {
  std::cout << "beam energy dependent t boundaries" << std::endl;
  //double Ebeam = 8.20787; // threshold
  double Ebeam = 11; // max
  TLorentzVector beam = {0, 0, Ebeam, Ebeam};
  TLorentzVector target = {0, 0, 0, M_PROTON};
  auto cm = beam + target;
  auto beta = cm.BoostVector();
  beam.Boost(-beta);
  target.Boost(-beta);
  cm.Boost(-beta);
  const double Etot = beam.E() + target.E();
  const double Ep = target.E();
  const double Pp = target.Vect().Mag();
  const double Epp = (Etot * Etot - M2_JPSI + M2_PROTON) / (2 * Etot);
  const double Ppp = std::sqrt(Epp * Epp - M2_PROTON);

  std::cout << "Threshold: "
            << (M2_JPSI + 2 * M_JPSI * M_PROTON) * ONE_OVER_2M_PROTON
            << std::endl;
  std::cout << "Ebeam: " << Ebeam << std::endl;

  std::cout << "tmin: " << 2 * M_PROTON - 2 * Ep * Epp - 2 * Pp * Ppp
            << std::endl;
  std::cout << "tmax: " << 2 * M_PROTON - 2 * Ep * Epp + 2 * Pp * Ppp
            << std::endl;
  std::cout << std::endl;
}

void test() {
   testbs();
  //tlim();
  //testjpsi();
}
