#include "eic_fastmc.h"

int test_eic_fastmc(){
  gStyle->SetOptStat(0);
  
  double kin[4];
  double kin_smeared[4];  
  double res[4];
  double acc[1];
  double pident[1];
  int region[1];  
  int pid=2212;
 
//   eic_fastmc(const char detector[],const char version[],const char ebeam[],const char ibeam[],const char solenoid[],bool debug=false);  
    eic_fastmc fun_eic_fastmc("det1","v1","ebeam_5GeV","ibeam_100GeV","solenoidoff",true);
  
       int n=1;	
      TH2F *hacc_ThetaPhi=new TH2F("acc_ThetaPhi","acc_ThetaPhi;vertex #phi (deg);vertex #theta (deg)",360*n,-180,180,180*n,0,180);

      for (Int_t j=-180*n;j<180*n;j++) {
      for (Int_t i=0;i<180*n;i++) {
	
	kin[0]=95; kin[1]=i/float(n); kin[2]=j/float(n); kin[3]=sqrt(95*95+0.938*0.938); 

// int eic_fastmc::get_eic_fastmc(const int pid,const double * kin,int *region,double *acc,double *res,double *kin_smeared,double *pident,bool debug)
	fun_eic_fastmc.get_eic_fastmc(pid,kin,region,acc,res,kin_smeared,pident,true);
	
      cout << "acc " << acc[0] << " region " << region[0] << "pident " << pident[0] << endl;	
      cout << "kin " << kin[0] << " " << kin[1] << " " << kin[2] << " " << kin[3] << endl;       
      cout << "kin_smeared " << kin_smeared[0] << " " << kin_smeared[1] << " " << kin_smeared[2]  << " " << kin_smeared[3] << endl;       
	
	hacc_ThetaPhi->SetBinContent(j+180*n+1,i+1,acc[0]);
      }
      }
      
      TCanvas *c_acc_total = new TCanvas("acc_total","acc_total",1000,900);
      hacc_ThetaPhi->SetMaximum(2);
      hacc_ThetaPhi->Draw("colz");
    
  return 0;
}
