#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include <TStyle.h>

namespace pcsim {
namespace jleic_impl {
using namespace std;

class acceptance {
public:
  acceptance()
      : acceptance("det1", "v1", "ebeam_5GeV", "ibeam_100GeV", "solenoidon") {}
  acceptance(const char detector[], const char version[], const char ebeam[],
             const char ibeam[], const char solenoid[], bool debug = false);
  int get_acceptance(const int pid, const double* kin, int* region, double* acc,
                     bool debug = false);

protected:
  Int_t _status;

  TFile* _file_eledownstream;
  TFile* _file_central;
  TFile* _file_iondownstream;
  TFile* _file_ionforward;

  TFile* _file_eledownstream_neu;
  TFile* _file_central_neu;
  TFile* _file_iondownstream_neu;
  TFile* _file_ionforward_neu;

  //   const int mm=9;
  //   char * title[mm] =
  //   {"acc_eledipole1side","acc_eledipole3back","acc_barrel","acc_ionside","acc_eleside_outer","acc_eleside_inner","acc_iondownstream","acc_iondipole2back","acc_iondipole3back"};

  TH3F* hacc[9];
  TH3F* hacc_neu[9];
};

#ifdef JLEIC_IMPL

#ifndef EIC_ACCEPTANCE_PATH
#error EIC_ACCEPTANCE_PATH NOT DEFINED
#endif

acceptance::acceptance(const char detector[], const char version[],
                       const char ebeam[], const char ibeam[],
                       const char solenoid[], bool debug) {
  if (debug)
    std::cout << "construct acceptance class" << std::endl;
  _status = 0;
  bool control = true;
  // file names
  TString feledownstream, fcentral, fiondownstream, fionforward;
  TString feledownstream_neu, fcentral_neu, fiondownstream_neu, fionforward_neu;
  if (strcmp(detector, "det1") == 0 && strcmp(version, "v1") == 0) {
    _status = 1;
    if (debug)
      std::cout << detector << " " << version << " " << ebeam << " " << ibeam
                << " " << solenoid << std::endl;
    if (strcmp(ebeam, "ebeam_5GeV") == 0 &&
        strcmp(ibeam, "ibeam_100GeV") == 0) {
      if (strcmp(solenoid, "solenoidoff") == 0) {
        feledownstream = EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_"
                                             "solenoidoff_electron_5GeV_"
                                             "eledownstream_output_final.root";
        fcentral = EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_solenoidoff_"
                                       "proton_100GeV_central_output_final."
                                       "root";
        fiondownstream = EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_"
                                             "solenoidoff_proton_100GeV_"
                                             "iondownstream_output_final.root";
        fionforward =
            EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_solenoidoff_"
                                "proton_100GeV_ionforward_output_final."
                                "root";
      }
      if (strcmp(solenoid, "solenoidon") == 0) {
        feledownstream = EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_"
                                             "solenoidon_electron_5GeV_"
                                             "eledownstream_output_final.root";
        fcentral =
            EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_solenoidon_proton_"
                                "100GeV_central_output_final.root";
        fiondownstream =
            EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_solenoidon_"
                                "proton_100GeV_iondownstream_output_"
                                "final.root";
        fionforward =
            EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_solenoidon_"
                                "proton_100GeV_ionforward_output_final."
                                "root";
      }

      feledownstream_neu = EIC_ACCEPTANCE_PATH
          "acceptance_det1_hybrid_neutron_5GeV_eledownstream_output_final.root";
      fcentral_neu = EIC_ACCEPTANCE_PATH
          "acceptance_det1_hybrid_neutron_100GeV_central_output_final.root";
      fiondownstream_neu = EIC_ACCEPTANCE_PATH "acceptance_det1_hybrid_neutron_"
                                               "100GeV_iondownstream_output_"
                                               "final.root";
      fionforward_neu = EIC_ACCEPTANCE_PATH
          "acceptance_det1_hybrid_neutron_100GeV_ionforward_output_final.root";
    }
  } else {
    std::cout << "Error in acceptance::acceptance!" << std::endl;
    std::cout << "  No configuration is matched" << std::endl;
    control = false;
  }

  if (control) {
    _file_eledownstream = new TFile(feledownstream, "r");
    _file_central = new TFile(fcentral, "r");
    _file_iondownstream = new TFile(fiondownstream, "r");
    _file_ionforward = new TFile(fionforward, "r");

    hacc[0] = (TH3F*)_file_eledownstream->GetObjectChecked("acc_eledipole1side",
                                                           "TH3F");
    hacc[1] = (TH3F*)_file_eledownstream->GetObjectChecked("acc_eledipole3back",
                                                           "TH3F");
    hacc[2] = (TH3F*)_file_central->GetObjectChecked("acc_barrel", "TH3F");
    hacc[3] = (TH3F*)_file_central->GetObjectChecked("acc_ionside", "TH3F");
    hacc[4] =
        (TH3F*)_file_central->GetObjectChecked("acc_eleside_outer", "TH3F");
    hacc[5] =
        (TH3F*)_file_central->GetObjectChecked("acc_eleside_inner", "TH3F");
    hacc[6] = (TH3F*)_file_iondownstream->GetObjectChecked("acc_iondownstream",
                                                           "TH3F");
    hacc[7] =
        (TH3F*)_file_ionforward->GetObjectChecked("acc_iondipole2back", "TH3F");
    hacc[8] =
        (TH3F*)_file_ionforward->GetObjectChecked("acc_iondipole3back", "TH3F");

    _file_eledownstream_neu = new TFile(feledownstream_neu, "r");
    _file_central_neu = new TFile(fcentral_neu, "r");
    _file_iondownstream_neu = new TFile(fiondownstream_neu, "r");
    _file_ionforward_neu = new TFile(fionforward_neu, "r");

    hacc_neu[0] = (TH3F*)_file_eledownstream_neu->GetObjectChecked(
        "acc_eledipole1side", "TH3F");
    hacc_neu[1] = (TH3F*)_file_eledownstream_neu->GetObjectChecked(
        "acc_eledipole3back", "TH3F");
    hacc_neu[2] =
        (TH3F*)_file_central_neu->GetObjectChecked("acc_barrel", "TH3F");
    hacc_neu[3] =
        (TH3F*)_file_central_neu->GetObjectChecked("acc_ionside", "TH3F");
    hacc_neu[4] =
        (TH3F*)_file_central_neu->GetObjectChecked("acc_eleside_outer", "TH3F");
    hacc_neu[5] =
        (TH3F*)_file_central_neu->GetObjectChecked("acc_eleside_inner", "TH3F");
    hacc_neu[6] = (TH3F*)_file_iondownstream_neu->GetObjectChecked(
        "acc_iondownstream", "TH3F");
    hacc_neu[7] = (TH3F*)_file_ionforward_neu->GetObjectChecked(
        "acc_iondipole2back", "TH3F");
    hacc_neu[8] = (TH3F*)_file_ionforward_neu->GetObjectChecked(
        "acc_iondipole3back", "TH3F");

    if (debug) {
      TCanvas* c_acc = new TCanvas("acc", "acc", 1000, 900);
      c_acc->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc->cd(i + 1);
        hacc[i]->Draw();
      }

      TH2F* hacc_ThetaPhi[9];
      hacc_ThetaPhi[0] = (TH2F*)_file_eledownstream->GetObjectChecked(
          "acc_ThetaPhi_eledipole1side", "TH2F");
      hacc_ThetaPhi[1] = (TH2F*)_file_eledownstream->GetObjectChecked(
          "acc_ThetaPhi_eledipole3back", "TH2F");
      hacc_ThetaPhi[2] =
          (TH2F*)_file_central->GetObjectChecked("acc_ThetaPhi_barrel", "TH2F");
      hacc_ThetaPhi[3] = (TH2F*)_file_central->GetObjectChecked(
          "acc_ThetaPhi_ionside", "TH2F");
      hacc_ThetaPhi[4] = (TH2F*)_file_central->GetObjectChecked(
          "acc_ThetaPhi_eleside_outer", "TH2F");
      hacc_ThetaPhi[5] = (TH2F*)_file_central->GetObjectChecked(
          "acc_ThetaPhi_eleside_inner", "TH2F");
      hacc_ThetaPhi[6] = (TH2F*)_file_iondownstream->GetObjectChecked(
          "acc_ThetaPhi_iondownstream", "TH2F");
      hacc_ThetaPhi[7] = (TH2F*)_file_ionforward->GetObjectChecked(
          "acc_ThetaPhi_iondipole2back", "TH2F");
      hacc_ThetaPhi[8] = (TH2F*)_file_ionforward->GetObjectChecked(
          "acc_ThetaPhi_iondipole3back", "TH2F");
      TCanvas* c_acc_ThetaPhi =
          new TCanvas("acc_ThetaPhi", "acc_ThetaPhi", 1000, 900);
      c_acc_ThetaPhi->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc_ThetaPhi->cd(i + 1);
        hacc_ThetaPhi[i]->Draw("colz");
      }

      TH2F* hacc_ThetaP[9];
      hacc_ThetaP[0] = (TH2F*)_file_eledownstream->GetObjectChecked(
          "acc_ThetaP_eledipole1side", "TH2F");
      hacc_ThetaP[1] = (TH2F*)_file_eledownstream->GetObjectChecked(
          "acc_ThetaP_eledipole3back", "TH2F");
      hacc_ThetaP[2] =
          (TH2F*)_file_central->GetObjectChecked("acc_ThetaP_barrel", "TH2F");
      hacc_ThetaP[3] =
          (TH2F*)_file_central->GetObjectChecked("acc_ThetaP_ionside", "TH2F");
      hacc_ThetaP[4] = (TH2F*)_file_central->GetObjectChecked(
          "acc_ThetaP_eleside_outer", "TH2F");
      hacc_ThetaP[5] = (TH2F*)_file_central->GetObjectChecked(
          "acc_ThetaP_eleside_inner", "TH2F");
      hacc_ThetaP[6] = (TH2F*)_file_iondownstream->GetObjectChecked(
          "acc_ThetaP_iondownstream", "TH2F");
      hacc_ThetaP[7] = (TH2F*)_file_ionforward->GetObjectChecked(
          "acc_ThetaP_iondipole2back", "TH2F");
      hacc_ThetaP[8] = (TH2F*)_file_ionforward->GetObjectChecked(
          "acc_ThetaP_iondipole3back", "TH2F");
      TCanvas* c_acc_ThetaP =
          new TCanvas("acc_ThetaP", "acc_ThetaP", 1000, 900);
      c_acc_ThetaP->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc_ThetaP->cd(i + 1);
        hacc_ThetaP[i]->Draw("colz");
      }

      TH2F* hacc_PhiP[9];
      hacc_PhiP[0] = (TH2F*)_file_eledownstream->GetObjectChecked(
          "acc_PhiP_eledipole1side", "TH2F");
      hacc_PhiP[1] = (TH2F*)_file_eledownstream->GetObjectChecked(
          "acc_PhiP_eledipole3back", "TH2F");
      hacc_PhiP[2] =
          (TH2F*)_file_central->GetObjectChecked("acc_PhiP_barrel", "TH2F");
      hacc_PhiP[3] =
          (TH2F*)_file_central->GetObjectChecked("acc_PhiP_ionside", "TH2F");
      hacc_PhiP[4] = (TH2F*)_file_central->GetObjectChecked(
          "acc_PhiP_eleside_outer", "TH2F");
      hacc_PhiP[5] = (TH2F*)_file_central->GetObjectChecked(
          "acc_PhiP_eleside_inner", "TH2F");
      hacc_PhiP[6] = (TH2F*)_file_iondownstream->GetObjectChecked(
          "acc_PhiP_iondownstream", "TH2F");
      hacc_PhiP[7] = (TH2F*)_file_ionforward->GetObjectChecked(
          "acc_PhiP_iondipole2back", "TH2F");
      hacc_PhiP[8] = (TH2F*)_file_ionforward->GetObjectChecked(
          "acc_PhiP_iondipole3back", "TH2F");
      TCanvas* c_acc_PhiP = new TCanvas("acc_PhiP", "acc_PhiP", 1000, 900);
      c_acc_PhiP->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc_PhiP->cd(i + 1);
        hacc_PhiP[i]->Draw("colz");
      }

      TCanvas* c_acc_neu = new TCanvas("acc_neu", "acc_neu", 1000, 900);
      c_acc_neu->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc_neu->cd(i + 1);
        hacc_neu[i]->Draw();
      }

      TH2F* hacc_neu_ThetaPhi[9];
      hacc_neu_ThetaPhi[0] = (TH2F*)_file_eledownstream_neu->GetObjectChecked(
          "acc_ThetaPhi_eledipole1side", "TH2F");
      hacc_neu_ThetaPhi[1] = (TH2F*)_file_eledownstream_neu->GetObjectChecked(
          "acc_ThetaPhi_eledipole3back", "TH2F");
      hacc_neu_ThetaPhi[2] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaPhi_barrel", "TH2F");
      hacc_neu_ThetaPhi[3] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaPhi_ionside", "TH2F");
      hacc_neu_ThetaPhi[4] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaPhi_eleside_outer", "TH2F");
      hacc_neu_ThetaPhi[5] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaPhi_eleside_inner", "TH2F");
      hacc_neu_ThetaPhi[6] = (TH2F*)_file_iondownstream_neu->GetObjectChecked(
          "acc_ThetaPhi_iondownstream", "TH2F");
      hacc_neu_ThetaPhi[7] = (TH2F*)_file_ionforward_neu->GetObjectChecked(
          "acc_ThetaPhi_iondipole2back", "TH2F");
      hacc_neu_ThetaPhi[8] = (TH2F*)_file_ionforward_neu->GetObjectChecked(
          "acc_ThetaPhi_iondipole3back", "TH2F");
      TCanvas* c_acc_neu_ThetaPhi =
          new TCanvas("acc_neu_ThetaPhi", "acc_neu_ThetaPhi", 1000, 900);
      c_acc_neu_ThetaPhi->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc_neu_ThetaPhi->cd(i + 1);
        hacc_neu_ThetaPhi[i]->Draw("colz");
      }

      TH2F* hacc_neu_ThetaP[9];
      hacc_neu_ThetaP[0] = (TH2F*)_file_eledownstream_neu->GetObjectChecked(
          "acc_ThetaP_eledipole1side", "TH2F");
      hacc_neu_ThetaP[1] = (TH2F*)_file_eledownstream_neu->GetObjectChecked(
          "acc_ThetaP_eledipole3back", "TH2F");
      hacc_neu_ThetaP[2] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaP_barrel", "TH2F");
      hacc_neu_ThetaP[3] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaP_ionside", "TH2F");
      hacc_neu_ThetaP[4] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaP_eleside_outer", "TH2F");
      hacc_neu_ThetaP[5] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_ThetaP_eleside_inner", "TH2F");
      hacc_neu_ThetaP[6] = (TH2F*)_file_iondownstream_neu->GetObjectChecked(
          "acc_ThetaP_iondownstream", "TH2F");
      hacc_neu_ThetaP[7] = (TH2F*)_file_ionforward_neu->GetObjectChecked(
          "acc_ThetaP_iondipole2back", "TH2F");
      hacc_neu_ThetaP[8] = (TH2F*)_file_ionforward_neu->GetObjectChecked(
          "acc_ThetaP_iondipole3back", "TH2F");
      TCanvas* c_acc_neu_ThetaP =
          new TCanvas("acc_neu_ThetaP", "acc_neu_ThetaP", 1000, 900);
      c_acc_neu_ThetaP->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc_neu_ThetaP->cd(i + 1);
        hacc_neu_ThetaP[i]->Draw("colz");
      }

      TH2F* hacc_neu_PhiP[9];
      hacc_neu_PhiP[0] = (TH2F*)_file_eledownstream_neu->GetObjectChecked(
          "acc_PhiP_eledipole1side", "TH2F");
      hacc_neu_PhiP[1] = (TH2F*)_file_eledownstream_neu->GetObjectChecked(
          "acc_PhiP_eledipole3back", "TH2F");
      hacc_neu_PhiP[2] =
          (TH2F*)_file_central_neu->GetObjectChecked("acc_PhiP_barrel", "TH2F");
      hacc_neu_PhiP[3] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_PhiP_ionside", "TH2F");
      hacc_neu_PhiP[4] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_PhiP_eleside_outer", "TH2F");
      hacc_neu_PhiP[5] = (TH2F*)_file_central_neu->GetObjectChecked(
          "acc_PhiP_eleside_inner", "TH2F");
      hacc_neu_PhiP[6] = (TH2F*)_file_iondownstream_neu->GetObjectChecked(
          "acc_PhiP_iondownstream", "TH2F");
      hacc_neu_PhiP[7] = (TH2F*)_file_ionforward_neu->GetObjectChecked(
          "acc_PhiP_iondipole2back", "TH2F");
      hacc_neu_PhiP[8] = (TH2F*)_file_ionforward_neu->GetObjectChecked(
          "acc_PhiP_iondipole3back", "TH2F");
      TCanvas* c_acc_neu_PhiP =
          new TCanvas("acc_neu_PhiP", "acc_neu_PhiP", 1000, 900);
      c_acc_neu_PhiP->Divide(3, 3);
      for (Int_t i = 0; i < 9; i++) {
        c_acc_neu_PhiP->cd(i + 1);
        hacc_neu_PhiP[i]->Draw("colz");
      }
    }
  }
}

int acceptance::get_acceptance(const int pid, const double* kin, int* region,
                               double* acc, bool debug) {
  //   if (CheckConfiguration()) {return 1;}
  // kin: p/GeV, theta/Deg, phi/Deg
  // res: dp/GeV, dtheta/Deg, dphi/Deg, dz/cm
  /* pid: e+: 1      e-: -1
         mu+: 2     mu-: -2
         pi+: 3     pi-: -3
          K+: 4      K-: -4
           p: 5
  */
  //   if (pid == 5){
  //     Rp = proton_p;
  //     Rtheta = proton_theta;
  //     Rphi = proton_phi;
  //     Rz = proton_z;
  //   }
  //   else {
  //     std::cout << "Error in acceptance::get_acceptance!" << std::endl;
  //     std::cout << "  No pid is matched." << std::endl;
  //     return 1;
  //   }

  if (kin[0] <= 0) {
    std::cout << "Error in acceptance::get_acceptance!" << std::endl;
    std::cout << "  Nonphysical momentum." << std::endl;
    return 1;
  }
  if (kin[1] > 180.0 || kin[1] < 0.0) {
    std::cout << "Error in acceptance::get_acceptance!" << std::endl;
    std::cout << "  Theta angle " << kin[1] << " out of range [0, 180]."
              << std::endl;
    return 1;
  }
  if (kin[2] > 180.0 || kin[2] < -180.0) {
    std::cout << "Error in acceptance::get_acceptance!" << std::endl;
    std::cout << "  Phi angle " << kin[2] << " out of range [0, 180]."
              << std::endl;
    return 1;
  }

  //   if (CheckKinematic(kin, pid)) {return 2;}

  acc[0] = 0;
  region[0] = -1;

  double acc_region[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int counter = 0;
  for (Int_t i = 0; i < 9; i++) {
    //       if(i!=4) continue;
    acc_region[i] =
        hacc[i]->GetBinContent(hacc[i]->GetXaxis()->FindBin(kin[0]),
                               hacc[i]->GetYaxis()->FindBin(kin[1]),
                               hacc[i]->GetZaxis()->FindBin(kin[2]));
    if (acc_region[i] > 0) {
      counter++;
      acc[0] += acc_region[i];
      region[0] = i;
    }
  }

  if (debug) {
    if (counter > 1) {
      cout << "counter " << counter << endl;
      for (Int_t i = 0; i < 9; i++) {
        if (acc_region[i] > 0)
          cout << acc_region[i] << " " << i << endl;
      }
    }
  }

  return 0;
}
#endif
} // namespace jleic_impl
} // namespace pcsim
