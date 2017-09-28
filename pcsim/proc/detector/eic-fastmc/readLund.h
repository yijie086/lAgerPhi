#include<iostream>
#include<fstream>
#include <vector>;

static int FIRST=0;
double tmp0, tmp, px, py, pz, En, tmp1,tmp2, Vx, Vy,Vz;
int pdef, type, parent, daughter;
int nparticles,bemaPol;


struct Particle {  
  double tmp0;
  int type;
  int parent;
  int pdef;
  int daughter;
  
  double px;
  double py;
  double pz;
  double En;
  double tmp2;
  double Vx;
  double Vy;
  double Vz;

   
 //.....
  
}; 
  


vector<Particle>  readLund(){
  ifstream myfile;
  vector<Particle> LundEvent ;
  Particle particle;
  
  if(FIRST==0) { FIRST=1;
    myfile.open("pythia.txt");
    if(!myfile)
      {
	std::cout << " Can't open input file " << "pythia.txt" << ". Exiting. " << std::endl;
	exit(1);
      }
    std::cout << " open input file " << "pythia.txt" << ". OK. " << std::endl;
  }
  
  //printf("=====> LUND 0 !!!!  \n");
  myfile >> nparticles ;
  std::cout << " nparticles" << nparticles << std::endl;
  
  for(unsigned i=0; i<9; i++)
    {        
      myfile >> tmp1;                             
    }
  
  for(int p=0; p<nparticles; p++)
    {
      myfile >> tmp0 >> tmp >> type >> pdef >> parent >> daughter >> px >> py >> pz >> En >> tmp2 >> Vx >> Vy >> Vz;
      std::cout << " LUND:LINE:: ip=" << tmp0 << " " <<  type << " "  << pdef << " "  << parent << " "  << daughter << " "  << px << " "  << py << " "  << pz << " "  << En << " "  << tmp2 << " "  << Vx  << " " << Vy << " "  << Vz  << " " << std::endl;

      particle.px=px;
      particle.py=py;
      particle.pz=pz;
      particle.En=En;
      particle.parent=parent;
      particle.daughter=daughter;
      particle.pdef=pdef;
      particle.type=type;
      particle.tmp0=tmp0;
      particle.tmp2=tmp2;
      particle.Vx=Vx;
      particle.Vy=Vy;
      particle.Vz=Vz;
 

      LundEvent.push_back(particle);
    }

  
  return LundEvent ;
  
}
