#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include <TStyle.h>

namespace pcsim {
namespace jleic_impl {
using namespace std;

class pidentification {
public:
  pidentification()
      : pidentification("det1", "v1", "ebeam_5GeV", "ibeam_100GeV",
                        "solenoidon") {}
  pidentification(const char detector[], const char version[],
                  const char ebeam[], const char ibeam[], const char solenoid[],
                  bool debug = false);
  int get_pidentification(const int pid, const double* kin, const int region,
                          double* pident, bool debug = false);

protected:
  Int_t _status;

  //   TFile *  _file_eledownstream;
  //   TFile *  _file_central;
  //   TFile *  _file_iondownstream;
  //   TFile *  _file_ionforward;

  //   const int mm=9;
  //   char * title[9] =
  //   {"acc_eledipole1side","acc_eledipole3back","acc_barrel","acc_ionside","acc_eleside_outer","acc_eleside_inner","acc_iondownstream","acc_iondipole2back","acc_iondipole3back"};
};

pidentification::pidentification(const char detector[], const char version[],
                                 const char ebeam[], const char ibeam[],
                                 const char solenoid[], bool debug) {
  if (debug)
    std::cout << "construct pidentification class" << std::endl;
  _status = 0;
  bool control = true;
  // file names
  TString feledownstream, fcentral, fiondownstream, fionforward;
  if (strcmp(detector, "det1") == 0 && strcmp(version, "v1") == 0) {
    _status = 1;
    if (debug)
      std::cout << detector << " " << version << " " << ebeam << " " << ibeam
                << " " << solenoid << std::endl;
    if (strcmp(ebeam, "ebeam_5GeV") == 0 &&
        strcmp(ibeam, "ibeam_100GeV") == 0) {
      if (strcmp(solenoid, "solenoidoff") == 0) {
      }
      if (strcmp(solenoid, "solenoidon") == 0) {
      }
    }
  } else {
    std::cout << "Error in pid::pid!" << std::endl;
    std::cout << "  No configuration is matched" << std::endl;
    control = false;
  }

  if (control) {
    //     _file_eledownstream = new TFile(feledownstream, "r");
    //     _file_central = new TFile(fcentral, "r");
    //     _file_iondownstream = new TFile(fiondownstream, "r");
    //     _file_ionforward = new TFile(fionforward, "r");
  }
}

int pidentification::get_pidentification(const int pid, const double* kin,
                                         const int region, double* pident,
                                         bool debug) {
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
  //     std::cout << "Error in pid::GetPid!" << std::endl;
  //     std::cout << "  No pid is matched." << std::endl;
  //     return 1;
  //   }

  if (kin[0] <= 0) {
    std::cout << "Error in pid::GetPid!" << std::endl;
    std::cout << "  Nonphysical momentum." << std::endl;
    return 1;
  }
  if (kin[1] > 180.0 || kin[1] < 0.0) {
    std::cout << "Error in pid::GetPid!" << std::endl;
    std::cout << "  Theta angle " << kin[1] << " out of range [0, 180]."
              << std::endl;
    return 1;
  }
  if (kin[1] > 180.0 || kin[1] < -180.0) {
    std::cout << "Error in pid::GetPid!" << std::endl;
    std::cout << "  Phi angle " << kin[2] << " out of range [0, 180]."
              << std::endl;
    return 1;
  }

  if (region < 0 || region > 8) {
    std::cout << "Error in pid::GetPid!" << std::endl;
    std::cout << "secton " << region << std::endl;
    return 1;
  }

  pident[0] = 0;

  //   char * title[9] =
  //   {"acc_eledipole1side","acc_eledipole3back","acc_barrel","acc_ionside","acc_eleside_outer","acc_eleside_inner","acc_iondownstream","acc_iondipole2back","acc_iondipole3back"};

  // delta_p/p(%), theta(deg), phi(deg), vertexz(cm)
  double DEG = 180. / TMath::Pi();
  double pident_region[9] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9};

  pident[0] = pident_region[region];

  if (debug)
    cout << "region " << region << " pident " << pident[0] << endl;

  return 0;
}
} // namespace jleic_impl
} // namespace pcsim
