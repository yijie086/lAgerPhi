#include "acceptance.h"
#include "resolution.h"
#include "pidentification.h"

int test_all(){
  gStyle->SetOptStat(0);
  
  double kin[4];
  double kin_smeared[4];  
  double res[4];
  double acc[1];
  double pident[1];
  int region[1];  
  int pid=2212;
  

    acceptance fun_acc("det1","v1","ebeam_5GeV","ibeam_100GeV","solenoidoff");
//     acceptance fun_acc("det1","v1","ebeam_5GeV","ibeam_100GeV","solenoidoff",true);  
  
       int n=1;	
      TH2F *hacc_ThetaPhi=new TH2F("acc_ThetaPhi","acc_ThetaPhi;vertex #phi (deg);vertex #theta (deg)",360*n,-180,180,180*n,0,180);

      for (Int_t j=-180*n;j<180*n;j++) {
      for (Int_t i=0;i<180*n;i++) {
	
	kin[0]=95; kin[1]=i/float(n); kin[2]=j/float(n); kin[3]=sqrt(95*95+0.938*0.938); 
// 	fun_acc.get_acceptance(pid,kin,region,acc,true);	
	fun_acc.get_acceptance(pid,kin,region,acc);		
// 	cout << "acc " << acc[0] << " region " << region[0] << endl;	
	hacc_ThetaPhi->SetBinContent(j+180*n+1,i+1,acc[0]);
      }
      }
      
      TCanvas *c_acc_total = new TCanvas("acc_total","acc_total",1000,900);
      hacc_ThetaPhi->SetMaximum(2);
      hacc_ThetaPhi->Draw("colz");

  resolution fun_res("det1","v1","ebeam_5GeV","ibeam_100GeV","solenoidoff",true);
  
//   fun_res.get_resolution(pid,kin,region[0],res,kin_smeared,true);		  //   
    fun_res.get_resolution(pid,kin,0,res,kin_smeared,true);		  
  
   pidentification  fun_pident("det1","v1","ebeam_5GeV","ibeam_100GeV","solenoidoff",true);

    fun_pident.get_pidentification(pid,kin,0,pident,true);		  
    
  return 0;
}
