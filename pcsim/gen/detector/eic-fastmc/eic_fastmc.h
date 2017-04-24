#include "eic_fastmc/acceptance.h"
#include "eic_fastmc/pidentification.h"
#include "eic_fastmc/resolution.h"

using namespace std;

class eic_fastmc {
public:
  eic_fastmc()
      : eic_fastmc("det1", "v1", "ebeam_5GeV", "ibeam_100GeV", "solenoidon") {}

  eic_fastmc(const char detector[], const char version[], const char ebeam[],
             const char ibeam[], const char solenoid[], bool debug = false);
  int get_eic_fastmc(const int pid, const double* kin, int* region, double* acc,
                     double* res, double* kin_smeared, double* pident,
                     bool debug = false);

protected:
  Int_t _status;

  //   TFile *  _file_eledownstream;
  //   TFile *  _file_central;
  //   TFile *  _file_iondownstream;
  //   TFile *  _file_ionforward;

  //   const int mm=9;
  //   char * title[9] =
  //   {"acc_eledipole1side","acc_eledipole3back","acc_barrel","acc_ionside","acc_eleside_outer","acc_eleside_inner","acc_iondownstream","acc_iondipole2back","acc_iondipole3back"};

  acceptance* fun_acc;
  resolution* fun_res;
  pidentification* fun_pident;
};

eic_fastmc::eic_fastmc(const char detector[], const char version[],
                       const char ebeam[], const char ibeam[],
                       const char solenoid[], bool debug) {
  if (debug)
    std::cout << "construct eic_fastmc class" << std::endl;
  _status = 0;
  bool control = true;

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
    std::cout << "Error in eic_fastmc::eic_fastmc!" << std::endl;
    std::cout << "  No configuration is matched" << std::endl;
    control = false;
  }

  if (control) {
    fun_acc = new acceptance(detector, version, ebeam, ibeam, solenoid, debug);
    fun_res = new resolution(detector, version, ebeam, ibeam, solenoid, debug);
    fun_pident =
        new pidentification(detector, version, ebeam, ibeam, solenoid, debug);
  }
}

int eic_fastmc::get_eic_fastmc(const int pid, const double* kin, int* region,
                               double* acc, double* res, double* kin_smeared,
                               double* pident, bool debug) {
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
  //     std::cout << "Error in resolution::get_eic_fastmc!" << std::endl;
  //     std::cout << "  No pid is matched." << std::endl;
  //     return 1;
  //   }

  //   int get_acceptance(const int pid,const double * kin, int *region, double
  //   *acc,bool debug=false);
  fun_acc->get_acceptance(pid, kin, region, acc, debug);

  //   int get_resolution(const int pid,const double * kin,const int
  //   region,double *res,double *kin_smeared,bool debug=false);
  fun_res->get_resolution(pid, kin, region[0], res, kin_smeared, debug);

  //   int get_pidentification(const int pid,const double * kin,const int
  //   region,double *pident,bool debug=false);
  fun_pident->get_pidentification(pid, kin, region[0], pident, debug);

  //   std::cout << "  test " << std::endl;

  return 0;
}
