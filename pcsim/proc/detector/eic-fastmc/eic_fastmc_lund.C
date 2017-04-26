#include "eic_fastmc.h"
#include "readLund.h"
#include <unistd.h>

int eic_fastmc_lund(){
  gStyle->SetOptStat(0);
  
  double kin[4];
  double kin_smeared[4];  
  double res[4];
  double acc[1];
  double pident[1];
  int region[1];  
  int pid=2212;

  // JF begin of one event  (Read lund file) .... 
  vector<Particle> LundEvent = readLund();

  printf("Event size=%d\n",LundEvent.size());
  // loop over particles in single event 
 
  for (ip=0; ip<LundEvent.size(); ip++) {
    
    printf("id =%d px=%f py=%f\n",LundEvent[ip].tmp0,LundEvent[ip].px,LundEvent[ip].py);

  }

  sleep(10);
  // JF  end of one event  (Read lund file) .... 

  
//   eic_fastmc(const char detector[],const char version[],const char ebeam[],const char ibeam[],const char solenoid[],bool debug=false);  
    eic_fastmc fun_eic_fastmc("det1","v1","ebeam_5GeV","ibeam_100GeV","solenoidoff",false);
  
       int n=1;	
      TH2F *hacc_ThetaPhi=new TH2F("acc_ThetaPhi","acc_ThetaPhi;vertex #phi (deg);vertex #theta (deg)",360*n,-180,180,180*n,0,180);

      for (Int_t j=-180*n;j<180*n;j++) {
      for (Int_t i=0;i<180*n;i++) {
	
	kin[0]=95; kin[1]=i/float(n); kin[2]=j/float(n); kin[3]=sqrt(95*95+0.938*0.938); 
	
// int eic_fastmc::get_eic_fastmc(const int pid,const double * kin,int *region,double *acc,double *res,double *kin_smeared,double *pident,bool debug)
	fun_eic_fastmc.get_eic_fastmc(pid,kin,region,acc,res,kin_smeared,pident,false);
	
      cout << "acc " << acc[0] << " region " << region[0] << endl;	
      cout << "kin_smeared " << kin_smeared[0] << " " << kin_smeared[1] << " " << kin_smeared[2]  << " " << kin_smeared[3] << endl;    
      cout << "pident " << pident[0] << endl; 
      
	
	hacc_ThetaPhi->SetBinContent(j+180*n+1,i+1,acc[0]);
      }
      }
      
      TCanvas *c_acc_total = new TCanvas("acc_total","acc_total",1000,900);
      hacc_ThetaPhi->SetMaximum(2);
      hacc_ThetaPhi->Draw("colz");
    
  return 0;
}
