#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>

using std::cout;
using std::endl;

Double_t Mp = 0.9383;
Double_t Mj = 3.0969;

Double_t dSigmadt(Double_t s, Double_t t) {
  // return values are in nb/GeV^2

  Double_t v = 1.0 / (16.0 * TMath::Pi() * pow((s - Mp * Mp), 2));
  Double_t x = (2.0 * Mp * Mj + Mj * Mj) / (s - Mp * Mp);

  Double_t A_2g =
      6.499e3 * v * pow(1 - x, 2) * pow((s - Mp * Mp), 2) / (Mj * Mj);
  Double_t A_3g =
      2.894e3 * v * pow(1 - x, 0) * pow((s - Mp * Mp), 2) / (Mj * Mj * Mj * Mj);

  // Below are fits to Cornell data (without Brodsky's 2 or 3 gluon exchange
  // model)
  // return(0.94*exp(t*0.97));//Cornell, E = 9.3-10.4
  // return(1.10*exp(t*1.31));//Cornell, E = 10.4-11.1
  Double_t ep = exp(t * 1.13);

  return ((A_2g)*ep); // 2-gluon model only
}

Double_t xsec(Double_t W) {
  // This function returns the cross-section for a 5/2- isospin pentaquark of
  // mass 4450 when given a W.  The function is parameterized from Qian Wang's
  // paper: Phys.Rev. D92 (2015) 034022.  It ignores any angular considerations
  // for the cross-section and uses a digitization of the plot in Fig 3.
  // Returned result is in nb.

  Double_t x[28] = {4.2402e+0, 4.2694e+0, 4.2969e+0, 4.3227e+0, 4.3392e+0,
                    4.3546e+0, 4.3681e+0, 4.3838e+0, 4.3937e+0, 4.4073e+0,
                    4.4190e+0, 4.4284e+0, 4.4372e+0, 4.4444e+0, 4.46e+0,
                    4.4713e+0, 4.4822e+0, 4.4954e+0, 4.5136e+0, 4.5277e+0,
                    4.5453e+0, 4.5767e+0, 4.6126e+0, 4.6508e+0, 4.6956e+0,
                    4.7450e+0, 4.7943e+0, 4.8257e+0};

  Double_t y[28] = {1.4915e-3, 4.4072e-3, 1.4628e-2, 4.5656e-2, 1.3052e-1,
                    3.5261e-1, 1.1059e+0, 3.5976e+0, 1.0937e+1, 3.9491e+1,
                    1.2079e+2, 4.9569e+2, 1.8390e+3, 9.1724e+3, 1.2697e+3,
                    2.9682e+2, 7.8895e+1, 2.1468e+1, 5.8243e+0, 1.9329e+0,
                    7.2579e-1, 1.9087e-1, 5.9198e-2, 2.1300e-2, 8.1890e-3,
                    3.5764e-3, 1.7529e-3, 1.2854e-3};
  TGraph* gr = new TGraph(28, x, y);
  if (W < x[0] || W > x[27]) {
    return 0;
  }
  Double_t ret = gr->Eval(W);
  delete gr;
  return ret;
}

Bool_t specAcc(Float_t centp, Float_t centth, Float_t hAcc, Float_t vAcc,
               Float_t dpAccmin, Float_t dpAccmax, TLorentzVector particle) {

  Bool_t ret = 0;
  Float_t pmin = centp - centp * dpAccmin;
  Float_t pmax = centp + centp * dpAccmax;
  // cout<< particle.Theta()<<endl;
  particle.RotateY(-centth);

  Float_t horzA = asin(particle.X() / particle.Vect().Mag());
  Float_t vertA = asin(particle.Y() / particle.Vect().Mag());
  // cout<<horzA<<" "<<vertA<<" "<<particle.Vect().Mag()<<endl;
  if (particle.Vect().Mag() > pmin && particle.Vect().Mag() < pmax &&
      TMath::Abs(horzA) < hAcc && TMath::Abs(vertA) < vAcc) {

    ret = 1;
  }

  return ret;
}

void JPsi(Int_t n, Int_t seed, Double_t HMS_th, Double_t HMS_centp,
          Double_t SHMS_th, Double_t SHMS_centp) {

  TH1D* hW = new TH1D("hW", "hW", 100, 4.0, 4.7);
  TH1D* hWw = new TH1D("hWw", "hWw", 100, 4.0, 4.7);
  TH2D* he =
      new TH2D("he", "electron momentum and angle", 100, -1, 3, 100, 0, 9.2);
  TH2D* hp =
      new TH2D("hp", "positron momentum and angle", 100, -1, 3, 100, 0, 9.2);

  TFile f("JPsi.root", "recreate");

  TTree t2("t2", "JPsi");
  TTree INFO("INFO", "pentaquark");

  Int_t nelec = 0;
  Int_t npositron = 0;
  Int_t ncoin = 0;
  Double_t Wj = 92.9e-6;
  Double_t Me = 0.000511;
  Double_t ctheta_min = -1.0;
  Double_t ctheta_max = 1.0;

  Double_t ctheta_cm;
  Double_t theta_cm;
  Double_t phi;
  Double_t ctheta_cm1;
  Double_t theta_cm1;
  Double_t phi1;

  Double_t Ecm;
  Double_t Ep_cm;
  Double_t Ej_cm;
  Double_t Epp_cm;
  Double_t Pp_cm;
  Double_t Pj_cm;
  Double_t Ppp_cm;
  Double_t Pe_cm;
  Double_t Pps_cm;

  Double_t theta_j_lab;
  Double_t theta_pp_lab;
  Double_t theta_lab_electron;
  Double_t theta_lab_positron;
  Double_t phi_lab_electron;
  Double_t phi_lab_positron;
  Double_t phi_lab_pp;
  Double_t Pj_lab;
  Double_t Ppp_lab;
  Double_t Pelectron_lab;
  Double_t Ppositron_lab;

  Double_t En;
  Double_t t;
  Double_t t_min;
  Double_t t_max;
  Double_t t_w;
  Double_t s;
  Double_t sigma;
  Double_t beam;

  t2.Branch("sigma", &sigma, "sigma/D");
  t2.Branch("beam", &beam, "beam/D");
  t2.Branch("Ppp_lab", &Ppp_lab, "Ppp_lab/D");
  t2.Branch("Ppositron_lab", &Ppositron_lab, "Ppositron_lab/D");
  t2.Branch("Pelectron_lab", &Pelectron_lab, "Pelectron_lab/D");
  t2.Branch("theta_pp_lab", &theta_pp_lab, "theta_pp_lab/D");
  t2.Branch("theta_lab_electron", &theta_lab_electron, "theta_lab_electron/D");
  t2.Branch("theta_lab_positron", &theta_lab_positron, "theta_lab_positron/D");
  t2.Branch("phi_lab_pp", &phi_lab_pp, "phi_lab_pp/D");
  t2.Branch("phi_lab_positron", &phi_lab_positron, "phi_lab_positron/D");
  t2.Branch("phi_lab_electron", &phi_lab_electron, "phi_lab_electron/D");

  INFO.Branch("n", &n, "n/I");
  INFO.Branch("HMS_centp", &HMS_centp, "HMS_centp/D");
  INFO.Branch("HMS_th", &HMS_th, "HMS_th/D");
  INFO.Branch("SHMS_centp", &SHMS_centp, "SHMS_centp/D");
  INFO.Branch("SHMS_th", &SHMS_th, "SHMS_th/D");

  INFO.Fill();

  TRandom3* r3 = new TRandom3();

  r3->SetSeed(seed);

  TF1* jpsiBW = new TF1(
      "jpsiBW", "2.0*sqrt(2.0)*[0]*[1]*sqrt([0]*[0]*([0]*[0] + "
                "[1]*[1]))/(TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0]*([0]*[0] + "
                "[1]*[1]))))/(pow(x*x - [0]*[0],2) + [0]*[0]*[1]*[1])",
      Mj - 10 * Wj, Mj + 10 * Wj);

  jpsiBW->SetParameters(Mj, Wj);

  TLorentzVector photon, p, j, pp, electron, positron;

  for (Int_t i = 0; i < n; i++) {

    Double_t Mj2 = jpsiBW->GetRandom();

    En = r3->Uniform((Mj2 * Mj2 + 2 * Mj2 * Mp) / (2 * Mp), 11.0);

    photon.SetXYZM(0.0, 0.0, En, 0.0);

    p.SetXYZM(0.0, 0.0, 0.0, Mp);

    TVector3 a = (photon + p).BoostVector();

    photon.Boost(-a);

    p.Boost(-a);

    Ecm = photon.E() + p.E();

    Pp_cm = p.Vect().Mag();

    Ep_cm = p.E();

    Ej_cm = (Ecm * Ecm + Mj2 * Mj2 - Mp * Mp) / (2 * Ecm);

    Pj_cm = TMath::Sqrt(Ej_cm * Ej_cm - Mj2 * Mj2);

    Epp_cm = (Ecm * Ecm + Mp * Mp - Mj2 * Mj2) / (2 * Ecm);

    Ppp_cm = TMath::Sqrt(Epp_cm * Epp_cm - Mp * Mp);

    t_min = 2 * Mp * Mp - 2 * (Ep_cm * Epp_cm - Pp_cm * Ppp_cm * ctheta_min);

    t_max = 2 * Mp * Mp - 2 * (Ep_cm * Epp_cm - Pp_cm * Ppp_cm * ctheta_max);

    t = r3->Uniform(t_min, t_max);

    ctheta_cm = (t + 2 * Ep_cm * Epp_cm - 2 * Mp * Mp) / (2 * Pp_cm * Ppp_cm);

    theta_cm = TMath::ACos(ctheta_cm);

    Pe_cm = (TMath::Sqrt(0.5 * Mj2 * 0.5 * Mj2 - Me * Me));
    Pps_cm = (TMath::Sqrt(0.5 * Mj2 * 0.5 * Mj2 - Me * Me));

    s = Mp * Mp + 2 * (photon.E()) * Mp;

    phi = r3->Uniform(0.0, 2 * TMath::Pi());

    Double_t x = TMath::Sin(theta_cm) * TMath::Cos(phi);
    Double_t y = TMath::Sin(theta_cm) * TMath::Sin(phi);
    Double_t z = TMath::Cos(theta_cm);

    pp.SetXYZM(Pj_cm * (x), Pj_cm * (y), Pj_cm * (z), Mp);
    j.SetXYZM(Ppp_cm * (-1) * (x), Ppp_cm * (-1) * (y), Ppp_cm * (-1) * (z),
              Mj2);

    ctheta_cm1 = r3->Uniform(-1.0, 1.0);
    phi1 = r3->Uniform(0.0, 2 * TMath::Pi());
    theta_cm1 = TMath::ACos(ctheta_cm1);

    Double_t x1 = TMath::Sin(theta_cm1) * TMath::Cos(phi1);
    Double_t y1 = TMath::Sin(theta_cm1) * TMath::Sin(phi1);
    Double_t z1 = TMath::Cos(theta_cm1);

    electron.SetXYZM((x1)*Pe_cm, (y1)*Pe_cm, (z1)*Pe_cm, Me);
    positron.SetXYZM((-1) * (x1)*Pe_cm, (-1) * (y1)*Pe_cm, (-1) * (z1)*Pe_cm,
                     Me);

    // boosting electron and positron to Jpsi in the com
    TVector3 c = j.BoostVector();
    electron.Boost(c);
    positron.Boost(c);

    // boosting everything to lab frame
    TVector3 b = p.BoostVector();

    photon.Boost(-b);

    p.Boost(-b);

    j.Boost(-b);

    pp.Boost(-b);

    electron.Boost(-b);

    positron.Boost(-b);

    Double_t W = (photon + p).M();

    t_w = t_max - t_min;

    sigma = dSigmadt(s, t) * (t_w);

    theta_j_lab = j.Theta();
    theta_pp_lab = pp.Theta();
    theta_lab_electron = electron.Theta();
    theta_lab_positron = positron.Theta();
    phi_lab_electron = electron.Phi();
    phi_lab_positron = positron.Phi();
    phi_lab_pp = pp.Phi();

    Pj_lab = j.Vect().Mag();
    Ppp_lab = pp.Vect().Mag();
    Pelectron_lab = electron.Vect().Mag();
    Ppositron_lab = positron.Vect().Mag();

    beam = photon.E();

    // Bool_t
    // HMS=specAcc(4.2,21.0*3.14159/180.0,24e-3,70e-3,0.10,0.10,electron);

    // Bool_t
    // SHMS=specAcc(5.8,-15.0*3.14159/180.0,20e-3,50e-3,0.15,0.25,positron);

    Bool_t HMS = specAcc(HMS_centp, HMS_th, 24e-3, 70e-3, 0.10, 0.10, electron);

    Bool_t SHMS =
        specAcc(SHMS_centp, SHMS_th, 20e-3, 50e-3, 0.15, 0.25, positron);
    if (HMS) {
      nelec++;
      hp->Fill(positron.Theta(), positron.Vect().Mag());
    }

    if (SHMS) {
      npositron++;
      he->Fill(electron.Theta(), electron.Vect().Mag());
    }

    if (HMS && SHMS) {
      ncoin++;

      t2.Fill();
    }

    hW->Fill(W);
    hWw->Fill(W, sigma);
  }

  t2.Write();
  INFO.Write();
  cout << nelec << " " << npositron << " " << ncoin << endl;

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 1200);
  TCanvas* c2 = new TCanvas("c2", "c2", 800, 1200);

  c1->Divide(1, 2);
  c2->Divide(1, 2);

  c1->cd(1);
  he->Draw();

  c1->cd(2);
  hp->Draw();

  hWw->Divide(hW);

  c2->cd(1);
  hW->Draw();

  c2->cd(2);
  hWw->Draw();

  TGraphErrors* SLAC1 = new TGraphErrors(3);
  TGraphErrors* SLAC2 = new TGraphErrors(4);

  SLAC1->SetFillColor(1);
  SLAC1->SetMarkerStyle(4);
  SLAC1->SetPoint(0, 4.4, 0.5417519);
  SLAC1->SetPointError(0, 0.12, 0.1896132);
  SLAC1->SetPoint(1, 4.58, 0.4605594);
  SLAC1->SetPointError(1, 0.06, 0.1611958);
  SLAC1->SetPoint(2, 4.73, 0.4730425);
  SLAC1->SetPointError(2, 0.07, 0.1655649);

  SLAC2->SetFillColor(1);
  SLAC2->SetMarkerStyle(26);
  SLAC2->SetPoint(0, 4.32, 0.12);
  SLAC2->SetPointError(0, 0, 0.08);
  SLAC2->SetPoint(1, 4.76, 0.225);
  SLAC2->SetPointError(1, 0, 0.085);
  SLAC2->SetPoint(2, 5.12, 1.65);
  SLAC2->SetPointError(2, 0, 0.35);
  SLAC2->SetPoint(3, 5.4, 2.83);
  SLAC2->SetPointError(3, 0, 0.37);

  SLAC1->Draw("sameP");
  SLAC2->Draw("sameP");
}

