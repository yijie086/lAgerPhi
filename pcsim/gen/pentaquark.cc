#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <iostream>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TF1.h"
#include "TVector3.h"

Double_t xsec(Double_t W){
  //This function returns the cross-section for a 5/2- isospin pentaquark of mass 4450 when given a W.  The function is parameterized from Qian Wang's paper: Phys.Rev. D92 (2015) 034022.  It ignores any angular considerations for the cross-section and uses a digitization of the plot in Fig 3.  Returned result is in nb.

  Double_t x[28] = {4.2402e+0, 4.2694e+0, 4.2969e+0, 4.3227e+0, 4.3392e+0, 4.3546e+0, 4.3681e+0, 4.3838e+0, 4.3937e+0, 4.4073e+0, 4.4190e+0, 4.4284e+0, 4.4372e+0, 4.4444e+0, 4.46e+0, 4.4713e+0, 4.4822e+0, 4.4954e+0, 4.5136e+0, 4.5277e+0, 4.5453e+0, 4.5767e+0, 4.6126e+0, 4.6508e+0, 4.6956e+0, 4.7450e+0, 4.7943e+0, 4.8257e+0};

  Double_t y[28] = {1.4915e-3, 4.4072e-3, 1.4628e-2, 4.5656e-2, 1.3052e-1, 3.5261e-1, 1.1059e+0, 3.5976e+0, 1.0937e+1, 3.9491e+1, 1.2079e+2, 4.9569e+2, 1.8390e+3, 9.1724e+3, 1.2697e+3, 2.9682e+2, 7.8895e+1, 2.1468e+1, 5.8243e+0, 1.9329e+0, 7.2579e-1, 1.9087e-1, 5.9198e-2, 2.1300e-2, 8.1890e-3, 3.5764e-3, 1.7529e-3, 1.2854e-3};
  TGraph *gr = new TGraph(28,x,y);
  if(W < x[0] || W > x[27]){
    return 0;
  }
  Double_t ret = gr->Eval(W);
  delete gr;
  return ret;
}


Bool_t specAcc(Float_t centp, Float_t centth, Float_t hAcc, Float_t vAcc, Float_t dpAccmin, Float_t dpAccmax, TLorentzVector particle)
{

Bool_t ret=0;
Float_t pmin=centp-centp*dpAccmin;
Float_t pmax=centp+centp*dpAccmax;
//cout<< particle.Theta()<<endl;
particle.RotateY(-centth);

Float_t horzA=asin(particle.X()/particle.Vect().Mag());
Float_t vertA=asin(particle.Y()/particle.Vect().Mag());
//cout<<horzA<<" "<<vertA<<" "<<particle.Vect().Mag()<<endl;
if (particle.Vect().Mag()>pmin && particle.Vect().Mag() < pmax && TMath::Abs(horzA)<hAcc && TMath::Abs(vertA)<vAcc)
{

  ret=1;
}

  return ret;
}

 void pentaquark(Int_t n, Int_t seed, Double_t HMS_th, Double_t HMS_centp, Double_t SHMS_th, Double_t SHMS_centp) {

  TH1D*h1= new TH1D("h1","Jpsi Breit Wigner Distribution",100,3.094,3.099);
  TH1D*h2= new TH1D("h2","Pentaquark Breit Wigner Distribution",100,4.2,4.7);
  TH2D*he= new TH2D("he","electron momentum and angle",100,-1,3,100,0,9.2);
  TH2D*hp= new TH2D("hp","positron momentum and angle",100,-1,3,100,0,9.2);
  TH1D*hW= new TH1D("hW","hW",100, 4.0,5.0);
  TH1D*hWw= new TH1D("hWw","hWw",100, 4.0,5.0);



  

  TFile f("pentaquark.root","recreate");
  TTree t1("t1","pentaquark");
  TTree INFO("INFO","pentaquark");

  Double_t Mp=0.9383;
  Double_t Mpc=4.4498;
  Double_t Mj=3.0969;
  Double_t Wj=92.9e-6;
  Double_t Wpc=39e-3;

  TVector3 a,b;
  TLorentzVector photon, p, pc, pp, pj;
  


  Int_t nelec=0;
  Int_t npositron=0;
  Int_t ncoin=0;
  Double_t theta_cm;
  Double_t theta_cm1;
  Double_t ctheta_cm;
  Double_t ctheta_cm1;
  Double_t phi;
  Double_t phi1;
  Double_t theta_lab_j;
  Double_t theta_lab_pp;
  Double_t theta_lab_electron;
  Double_t theta_lab_positron;
  Double_t phi_lab_electron;
  Double_t phi_lab_positron;
  Double_t phi_lab_pp;
  Double_t Pp_lab;
  Double_t Pj_lab;
  Double_t Ppp_lab;
  Double_t Ppositron_lab;
  Double_t Pelectron_lab;
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t x1;
  Double_t y1;
  Double_t z1;
  Float_t xA;
  Float_t yA;
  Double_t cs;
  Double_t W;
  Double_t beam;
  
  

  Double_t Me=0.000511;

  TLorentzVector electron, positron;


 


  
  t1.Branch("cs",&cs,"cs/D");
  t1.Branch("beam",&beam,"beam/D");
  t1.Branch("Ppp_lab",&Ppp_lab,"Ppp_lab/D");
  t1.Branch("Ppositron_lab",&Ppositron_lab,"Ppositron_lab/D");
  t1.Branch("Pelectron_lab",&Pelectron_lab,"Pelectron_lab/D");
  t1.Branch("theta_lab_pp",&theta_lab_pp,"theta_lab_pp/D");
  t1.Branch("theta_lab_electron",&theta_lab_electron,"theta_lab_electron/D");
  t1.Branch("theta_lab_positron",&theta_lab_positron,"theta_lab_positron/D");
  t1.Branch("phi_lab_electron",&phi_lab_electron,"phi_lab_electron/D");
  t1.Branch("phi_lab_positron",&phi_lab_positron,"phi_lab_positron/D");
  t1.Branch("phi_lab_pp",&phi_lab_pp,"phi_lab_pp/D");

  INFO.Branch("n",&n,"n/I");
  INFO.Branch("HMS_centp",&HMS_centp,"HMS_centp/D");
  INFO.Branch("HMS_th",&HMS_th,"HMS_th/D");
  INFO.Branch("SHMS_centp",&SHMS_centp,"SHMS_centp/D");
  INFO.Branch("SHMS_th",&SHMS_th,"SHMS_th/D");
 
  INFO.Fill();


  
  
  TRandom3*r3=new TRandom3 ();



  TF1 *jpsiBW= new TF1("jpsiBW","2.0*sqrt(2.0)*[0]*[1]*sqrt([0]*[0]*([0]*[0] + [1]*[1]))/(TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0]*([0]*[0] + [1]*[1]))))/(pow(x*x - [0]*[0],2) + [0]*[0]*[1]*[1])",Mj - 10*Wj, Mj + 10*Wj);
  jpsiBW->SetParameters(Mj,Wj);
  TF1 *pcBW= new TF1("pc","2.0*sqrt(2.0)*[0]*[1]*sqrt([0]*[0]*([0]*[0] + [1]*[1]))/(TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0]*([0]*[0] + [1]*[1]))))/(pow(x*x - [0]*[0],2) + [0]*[0]*[1]*[1])",Mpc - 10*Wpc, Mpc + 10*Wpc);
  pcBW->SetParameters(Mpc,Wpc);

  

  for (Int_t i=0;i<n;i++) {

     Double_t Mj2=jpsiBW->GetRandom();
     Double_t Mpc2=pcBW->GetRandom();
     Double_t Ep_cm=((Mpc2*Mpc2 - Mj2*Mj2+Mp*Mp)/(2*Mpc2));
     Double_t Ej_cm=((Mpc2*Mpc2+Mj2*Mj2-Mp*Mp)/(2*Mpc2));
     Double_t Pp_cm=(TMath::Sqrt(Ep_cm*Ep_cm - Mp*Mp));
     Double_t Pj_cm=(TMath::Sqrt(Ej_cm*Ej_cm - Mj2*Mj2));
     Double_t Pe_cm=(TMath::Sqrt(0.5*Mj2*0.5*Mj2-Me*Me));
     Double_t Pps_cm=(TMath::Sqrt(0.5*Mj2*0.5*Mj2-Me*Me));


    //In pentaquark rest/CoM frame
    r3->SetSeed(i+10456);
    ctheta_cm=r3->Uniform(-1.0,1.0);
    phi=r3->Uniform(0.1,2*TMath::Pi());
    theta_cm=TMath::ACos(ctheta_cm);
  
    ctheta_cm1=r3->Uniform(-1.0,1.0);
    phi1=r3->Uniform(0.0,2*TMath::Pi());
    theta_cm1=TMath::ACos(ctheta_cm1);

    x=TMath::Sin(theta_cm)*TMath::Cos(phi);
    y=TMath::Sin(theta_cm)*TMath::Sin(phi);
    z=TMath::Cos(theta_cm);

    x1=TMath::Sin(theta_cm1)*TMath::Cos(phi1);
    y1=TMath::Sin(theta_cm1)*TMath::Sin(phi1);
    z1=TMath::Cos(theta_cm1);

    photon.SetXYZM(0.0,0.0,((Mpc2*Mpc2-Mp*Mp)/(2*Mpc2)),0.0);
    p.SetXYZM(0.0,0.0,(-1)*((Mpc2*Mpc2-Mp*Mp)/(2*Mpc2)),Mp);
    pc.SetPxPyPzE(0.0,0.0,0.0,Mpc2);

    a=p.BoostVector();
  
    photon.Boost(-a);
    p.Boost(-a);

    //cout<<p.E()<< "" <<photon.E()<<endl; 

    //cout<<photon.E()<<endl;

    W=(photon+p).M(); 
    cs=xsec(W);


    pc.SetPxPyPzE(0.0,0.0,0.0,Mpc2);
    pp.SetXYZM(Pp_cm*(-1)*(x),Pp_cm*(-1)*(y),Pp_cm*(-1)*(z),Mp);
    pj.SetXYZM(Pj_cm*(x),Pj_cm*(y),Pj_cm*(z),Mj2);



    //making electron and positron in J/psi rest frame.
    electron.SetXYZM((x1)*Pe_cm,(y1)*Pe_cm,(z1)*Pe_cm,Me);
    positron.SetXYZM((-1)*(x1)*Pe_cm,(-1)*(y1)*Pe_cm,(-1)*(z1)*Pe_cm,Me);


    //boosting electron and positron to Pentaquark/CoM frame
    b=pj.BoostVector();

    electron.Boost(b);
    positron.Boost(b);


    //Boost everything back to lab frame

    electron.Boost(-a);
    positron.Boost(-a);
    pj.Boost(-a);
    pc.Boost(-a);
    pp.Boost(-a);

    theta_lab_j = pj.Theta();
    theta_lab_pp = pp.Theta();
    theta_lab_electron=electron.Theta();
    theta_lab_positron=positron.Theta();
    phi_lab_electron=electron.Phi();
    phi_lab_positron=positron.Phi();
    phi_lab_pp=pp.Phi();
    Pj_lab = pj.Vect().Mag();
    Ppp_lab = pp.Vect().Mag();
    Pelectron_lab = electron.Vect().Mag();
    Ppositron_lab = positron.Vect().Mag();
    beam=photon.E();
    //cout<<(photon + p).E()<< "  " <<(pp + electron + positron).E()<<endl;
  

    
  
    
    //Bool_t HMS=specAcc(3.65,0.35,24e-3,70e-3,0.10,0.10,electron);

    //Bool_t SHMS=specAcc(3.95,-0.45,20e-3,50e-3,0.15,0.25,positron);

    Bool_t HMS=specAcc(HMS_centp,HMS_th,24e-3,70e-3,0.10,0.10,electron);

    Bool_t SHMS=specAcc(SHMS_centp,SHMS_th,20e-3,50e-3,0.15,0.25,positron);

   
    if(HMS){
      nelec++;
    hp->Fill(positron.Theta(),positron.Vect().Mag());
}


    if(SHMS){
      npositron++;
    he->Fill(electron.Theta(),electron.Vect().Mag());
}

    if(HMS&&SHMS){
    ncoin++;

    t1.Fill();
   
}
    
    
    h1->Fill((electron+positron).M());
    h2->Fill((pp+pj).M());
    hW->Fill(W);
    hWw->Fill(W,cs);
    
  }


    
  t1.Write();
  INFO.Write();
  cout<<nelec<<" "<<npositron<<" "<<ncoin<<endl;

  TCanvas*c1=new TCanvas("c1","c1",800,1200);
  TCanvas*c2=new TCanvas("c2","c2",800,1200);
  
  c1->Divide(1,2);
  c2->Divide(1,2);
 

  c1->cd(1);
  h1->Draw();
  h1->GetXaxis()->SetTitle("(electron+positron) mass (GeV)");
  
  c1->cd(2);
  h2->Draw();
  h2->GetXaxis()->SetTitle("(jpsi+pp) mass (GeV)");

  c1->cd(1);
  he->Draw();
  
  c1->cd(2);
  hp->Draw();

  hWw->Divide(hW);
 
  c2->cd(1);
  hW->Draw();

  c2->cd(2);
  hWw->Draw();

  
  



}
