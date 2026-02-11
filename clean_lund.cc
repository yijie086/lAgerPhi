// =====================================================
// clean_lund.cc
//
// Clean LUND file + physics rejection sampling
// No external dependency (no exprtk)
//
// Compile:
//   g++ -std=c++17 -O2 clean_lund.cc -o clean_lund
//
// Usage:
//   ./clean_lund in.lund out.lund FMAX summary.txt tag
// =====================================================

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <random>

static constexpr double Mp = 0.9382720813;

// ---------- particle struct ----------
struct Particle {
  int idx;
  double d2;
  int ist, pid, parent, daughter;
  double px,py,pz,E,m,vx,vy,vz;
};

// ---------- bad number ----------
inline bool bad(double x){
  return std::isnan(x)||!std::isfinite(x);
}

// ---------- physics function ----------
// >>>>>>>>> MODIFY THIS ONLY <<<<<<<<<<
//
// Example DVCS-friendly function:
//
//  - suppress large Q2
//  - prefer moderate xB
//  - flat in others
//
double f_phys(double Q2,double xB,double y,double W,double t)
{
  double a = 0.431061 - 0.569157/xB + 0.216588 /(xB*xB) - 0.017522 /(xB*xB*xB);
  return 0.1/a;
}

// ---------- kinematics ----------
bool computeKin(
  const std::vector<Particle>& ps,
  double beamE,
  double& Q2,double& xB,double& y,double& W,double& t)
{
  const Particle* e=nullptr;
  const Particle* p=nullptr;

  for(auto& x:ps){
    if(x.pid==11 && x.ist==1 && !e) e=&x;
    if(x.pid==2212 && x.ist==1 && !p) p=&x;
  }

  if(!e) return false;

  double Ee=e->E;
  double pe=std::sqrt(e->px*e->px+e->py*e->py+e->pz*e->pz);
  double cosT=e->pz/pe;

  Q2=2*beamE*Ee*(1-cosT);
  y=(beamE-Ee)/beamE;
  xB=Q2/(2*Mp*(beamE-Ee));
  W=std::sqrt(Mp*Mp+2*Mp*(beamE-Ee)-Q2);

  if(p){
    double dqE=beamE-Ee;
    double dqx=-e->px;
    double dqy=-e->py;
    double dqz=beamE-e->pz;

    double dx=p->px-dqx;
    double dy=p->py-dqy;
    double dz=p->pz-dqz;
    double dE=p->E-dqE;

    t=dE*dE-(dx*dx+dy*dy+dz*dz);
  }
  else t=NAN;

  return !(bad(Q2)||bad(xB)||bad(y)||bad(W));
}

// =====================================================
// main
// =====================================================
int main(int argc,char* argv[])
{
  if(argc<6){
    std::cerr<<"Usage:\n"
             <<"clean_lund in.lund out.lund FMAX summary.txt tag\n";
    return 1;
  }

  std::string inF=argv[1];
  std::string outF=argv[2];
  double FMAX=std::stod(argv[3]);
  std::string summary=argv[4];
  std::string tag=argv[5];

  std::ifstream fin(inF);
  std::ofstream fout(outF);

  if(!fin||!fout){
    std::cerr<<"File error\n";
    return 2;
  }

  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> U(0,1);

  long long accepted=0, generated=0;

  std::string line;

  while(std::getline(fin,line))
  {
    if(line.empty()||line[0]=='#'){
      fout<<line<<"\n";
      continue;
    }

    std::istringstream hs(line);

    int npart,i2,i3,bpid,i8,i9;
    double d4,d5,beamE,w;

    hs>>npart>>i2>>i3>>d4>>d5>>bpid>>beamE>>i8>>i9>>w;
    if(!hs) continue;

    std::vector<Particle> ps;
    bool badEv=bad(beamE);

    for(int i=0;i<npart;i++)
    {
      std::getline(fin,line);
      std::istringstream psx(line);

      Particle p;

      psx>>p.idx>>p.d2>>p.ist>>p.pid>>p.parent>>p.daughter
         >>p.px>>p.py>>p.pz>>p.E>>p.m>>p.vx>>p.vy>>p.vz;

      if(!psx||bad(p.px)||bad(p.py)||bad(p.pz)||bad(p.E))
        badEv=true;

      ps.push_back(p);
    }

    generated++;
    if(badEv) continue;

    double Q2,xB,y,W,t;
    if(!computeKin(ps,beamE,Q2,xB,y,W,t)) continue;

    double fv=f_phys(Q2,xB,y,W,t);
    double P=fv/FMAX;
    if(P<0) P=0;
    if(P>1) P=1;
    if(U(rng)>P) continue;

    // ---------- write event ----------
    fout<<npart<<" "<<i2<<" "<<i3<<" "<<d4<<" "<<d5<<" "
        <<bpid<<" "<<beamE<<" "<<i8<<" "<<i9<<" "<<w<<"\n";

    int idx=1;
    for(auto& p:ps){
      fout<<idx++<<" "<<p.d2<<" "<<p.ist<<" "<<p.pid<<" "
          <<p.parent<<" "<<p.daughter<<" "
          <<p.px<<" "<<p.py<<" "<<p.pz<<" "<<p.E<<" "
          <<p.m<<" "<<p.vx<<" "<<p.vy<<" "<<p.vz<<"\n";
    }

    accepted++;
  }

  std::cout<<"ACCEPTED_EVENTS="<<accepted<<"\n";

  std::ofstream rep(summary,std::ios::app);
  rep<<"tag="<<tag<<" generated="<<generated
     <<" accepted="<<accepted<<"\n";
}
