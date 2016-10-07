// compile with 
// g++ -o ../build/simRead_andrea simRead_andrea.cpp `root-config --cflags --glibs`
// syntax
// simRead_base `ls out*`

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TObjArray.h"
#include "TObject.h"
#include <algorithm>    // std::sort
#include <numeric>      // std::accumulate
#include <vector>  
#include "TH1F.h"
#include "TCanvas.h"

	struct fotone_ottico
	{
	    float energia;
	    float tempo;
		float x;
		float y;
		
		bool operator<(const fotone_ottico& a) const 
		{
    		return tempo < a.tempo;
		}

	};
	


int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>"); //needed by ROOT to deal with standard vectors
  
  //-------------------
  // Input Files
  //-------------------
  TChain *tree =  new TChain("tree"); // read input files
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  // find the number of channels directly from the tchain file
  // before creating the variables
  // first, get the list of leaves
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
  for(int i = 0 ; i < nLeaves ; i++)
  {
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  int numOfCry = 0;
  std::string det_prefix("detector");
  std::string cry_prefix("cry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    //     leavesName.push_back(leavescopy->At(i)->GetName());
    if (!leavesName[i].compare(0, det_prefix.size(), det_prefix))
      numOfCh++;
    if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
      numOfCry++;
  }
  //the string "cry" appears 4 times per crystal..
  numOfCry = numOfCry / 4;
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;
  std::cout << "Number of Crystals \t= "<< numOfCry << std::endl;
  
  
  //------------------
  // Input TTree
  //------------------
  
  //create the branches in the input ttree and connect to the variables
  
  // global variables
  // these are 1 number per TTree entry - so 1 number per gamma shot
  Long64_t Seed;                      // seed of the simulation (read every time, but always the same) 
  int Run;                            // run id (usually just 1)(read every time, but always the same) 
  int Event;                          // event id
  float totalEnergyDeposited;         // total energy deposited in this event, in all the matrix
  int NumOptPhotons;                  // number of optical photons generated in this event, in the entire matrix
  int NumCherenkovPhotons;            // number of Cherenkov photons generated in this event, in the entire matrix
  
  //energy deposition events
  // these are 1 array of std::vector<float> per TTree entry. Explanation:
  // one TTree entry is 1 gamma shot, but
  // since the energy deposition events per gamma shot can be whatever, from 0 to N
  // we saved the data in std::vectors since they have variable size
  // if in a particular sim event the gamma deposited energy in N different points
  // the std::vectors are effectively arrays of length N
  // holding info of each energy deposition event, in this way
  // pEdep[N] = E0,E1,E2,E3,...,EN
  // px[N]    = x0,x1,x2,x3,...,xN
  // py[N]    = y0,y1,y2,y3,...,yN
  // pz[N]    = z0,z1,z2,z3,...,zN
  // where En is the energy deposited by the gamma in the n-th interaction, xn the x position of this interaction, and so on.
  // Just to complicate life (not really, in fact to make it simpler) the std::vectors are already sorted per crystal here, so
  // connecting them to the TBranches is easier. Remember that in the output from the simulation, the vectors above are TTree entries called
  // cryN
  // cryNPosXEnDep
  // cryNPosYEnDep
  // cryNPosZEnDep
  // here we need to store these values in some variables once we read them while looping on the TTree, and we don't know 
  // a priori the number of crystals. So instead of a simple pEdep std::vector, we need an array of pEdep std::vectors, with 
  // dynamic length (and same holds for px, py and pz). Hence, We don't declare a std::vector, but a pointer to std::vector, which
  // allows to make an array of std::vectors run time.
  // Finally, ROOT comes in the way and forces us not to pass a std::vector to the SetBranchAddress function, but a f...ing **ptr because
  // of template classes. All very complicated but in the end this is how it work:
  // in the event loop, after GetEvent(i) the variable pEdep[i] is a std::vector for the i-th crystal, that has N elements, one per each 
  // energy deposition occurred in that crystal.
  std::vector<float> **pEdep;         // for each energy deposition event, the amount of energy deposited
  std::vector<float> **px;            // for each energy deposition event, the x position 
  std::vector<float> **py;            // for each energy deposition event, the y position 
  std::vector<float> **pz;            // for each energy deposition event, the z position
  // create the arrays 
  pEdep = new std::vector<float>* [numOfCry];
  px    = new std::vector<float>* [numOfCry];
  py    = new std::vector<float>* [numOfCry];
  pz    = new std::vector<float>* [numOfCry];
  // inizialize to 0... or you'll have bad surprises
  for (int i = 0 ; i < numOfCry ; i++)
  {
    pEdep[i] = 0; 
    px[i] = 0;
    py[i] = 0;
    pz[i] = 0;
  }
  
  //Total number of photons detected in this event
  // for each TTree entry, a simple number saying how many optical photons entered that
  // specific detector, passed the PDE check and where "detected" (i.e. saved)
  Short_t  *detector;
  detector = new Short_t [numOfCh];
  
  //hits on detectors events 
  //these are 1 std::vector<float> per TTree entry
  //For each gamma shot, a great number of opticals is generated in the matrix
  //some of them enter the detectors, pass the PDE check and are "detected" (i.e. saved). So
  //each one of these std::vector has length M, where M is the number of optical photons
  //saved for that gamma event. Notice that this M is DIFFERENT from N of the energy deposition part.
  //data is then arranged as above, so for example
  // pOpticalX[M]               = x0,x1,x2,x3,...,xN
  // pOpticalPreMomentumX[M]    = m0,m1,m2,m3,...,mN
  // pOpticalPostMomentumX[M]   = v0,v1,v2,v3,...,vN
  // pGlobalTime[M]             = t0,t1,t2,t3,...,tN
  // and so on.
  
  //Position of optical photon when entering detector
  std::vector<float> *pOpticalX;            // x position of optical photon when entering detector
  std::vector<float> *pOpticalY;            // y position of optical photon when entering detector
  std::vector<float> *pOpticalZ;            // z position of optical photon when entering detector
  //inizialize...
  pOpticalX = 0;
  pOpticalY = 0;
  pOpticalZ = 0;
  //Momentum before and after entering the detector
  //before
  std::vector<float> *pOpticalPreMomentumX;            // x component of momentum unitary vector of optical photon before entering detector
  std::vector<float> *pOpticalPreMomentumY;            // y component of momentum unitary vector of optical photon before entering detector
  std::vector<float> *pOpticalPreMomentumZ;            // z component of momentum unitary vector of optical photon before entering detector
  //inizialize...
  pOpticalPreMomentumX = 0;
  pOpticalPreMomentumY = 0;
  pOpticalPreMomentumZ = 0;
  //after
  std::vector<float> *pOpticalPostMomentumX;            // x component of momentum unitary vector of optical photon after entering detector
  std::vector<float> *pOpticalPostMomentumY;            // y component of momentum unitary vector of optical photon after entering detector
  std::vector<float> *pOpticalPostMomentumZ;            // z component of momentum unitary vector of optical photon after entering detector
  //inizialize...
  pOpticalPostMomentumX = 0;
  pOpticalPostMomentumY = 0;
  pOpticalPostMomentumZ = 0;
  //Global time (i.e. from emission of primary gamma for this event) of detection for the optical
  std::vector<float> *pGlobalTime;
  pGlobalTime = 0;//inizialize...
  //Type of optical photon (0 = scintillation, 1 = Cherenkov, 2 = other (should not exist))
  std::vector<int> *pPhotonType;
  pPhotonType = 0;//inizialize...
  //Energy of optical photon (in eV)
  std::vector<float> *pPhotonEnergy;
  pPhotonEnergy = 0;//inizialize...
  
  
  //------------------------
  // Set Branch Addresses
  //------------------------
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  for (int i = 0 ; i < numOfCry ; i++)
  {
    std::stringstream snames;
    snames << "cry" << i;
    tree->SetBranchAddress(snames.str().c_str(),&pEdep[i]);
    snames.str("");
    snames<< "cry" << i << "PosXEnDep";    
    tree->SetBranchAddress(snames.str().c_str(),&px[i]);
    snames.str("");
    snames<< "cry" << i << "PosYEnDep";
    tree->SetBranchAddress(snames.str().c_str(),&py[i]);
    snames.str("");
    snames<< "cry" << i << "PosZEnDep";
    tree->SetBranchAddress(snames.str().c_str(),&pz[i]);
  }
  for (int i = 0 ; i < numOfCh ; i++) 
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }
  tree->SetBranchAddress("PositionX",&pOpticalX);
  tree->SetBranchAddress("PositionY",&pOpticalY);
  tree->SetBranchAddress("PositionZ",&pOpticalZ);
  
  tree->SetBranchAddress("PreMomentumX",&pOpticalPreMomentumX);
  tree->SetBranchAddress("PreMomentumY",&pOpticalPreMomentumY);
  tree->SetBranchAddress("PreMomentumZ",&pOpticalPreMomentumZ);
  tree->SetBranchAddress("PostMomentumX",&pOpticalPostMomentumX);
  tree->SetBranchAddress("PostMomentumY",&pOpticalPostMomentumY);
  tree->SetBranchAddress("PostMomentumZ",&pOpticalPostMomentumZ);
  
  tree->SetBranchAddress("GlobalTime",&pGlobalTime);
  tree->SetBranchAddress("PhotonType",&pPhotonType);
  tree->SetBranchAddress("PhotonEnergy",&pPhotonEnergy);
  
  
  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  
  TH1F *timeHisto[numOfCry];
  TCanvas *timeCanvas[numOfCry];
  
  for(int i = 0; i < numOfCry ; i++)
  {
		std::ostringstream s;
		s << "cristallo_" << i;
		std::string name(s.str());
 		timeHisto[i]= new TH1F(name.c_str(),name.c_str(), 100, 0,2);
		name+="_c";
		timeCanvas[i]=new TCanvas(name.c_str(),name.c_str(),1200,800);
  }
		

  for(int iEvent = 256; iEvent < 257 ; iEvent++)
  {  
		
		tree->GetEvent(iEvent);
		int CrystalsHit = 0;      // counter of how many crystals were hit by this gamma
		
		float* energyPerCrystal;  // total energy deposited in each crystal
		energyPerCrystal = new float[numOfCry];
		for(int i = 0; i < numOfCry ; i++)
		{
		energyPerCrystal[i] = 0;
		}
		
		// some example of global variable
		std::cout << "Event Numb \t\t\t = "              << Event << std::endl;
		std::cout << "Total energy deposited \t\t = "  << totalEnergyDeposited << std::endl;
		std::cout << "Total Numb of Opticals produced\t = "  << NumOptPhotons+NumCherenkovPhotons << std::endl;
		
		int primi_fotoni = 10;

		//creo un array di struct, tante quante sono i fotoni ottici; ogni elemento dell'array contiene energia, tempo di arrivo, x e y di rivelazione
		fotone_ottico fotone[pGlobalTime->size()];
		if(pGlobalTime->size())
		{
			for (int i=0; i<pGlobalTime->size(); i++) 
			{
				fotone[i].energia = pPhotonEnergy->at(i);
				fotone[i].tempo = pGlobalTime->at(i);
				fotone[i].x = pOpticalX->at(i);
				fotone[i].y = pOpticalY->at(i);
			}
			std::sort(fotone, fotone+pGlobalTime->size());
		}

		// one example of energy deposition related std::vectors
		for(int i = 0; i < numOfCry ; i++) //crystal by crystal
		{
		if(px[i]->size()) //if the i-th crystal had energy deposition (--> if(0) => false) 
		{
		
			std:: cout << "Le posizioni in cui il gamma ha deposto energia sono: "<<std::endl;    /////
			CrystalsHit++; //count the crystal as hit
			for(int j = 0; j < px[i]->size(); j++) //run on the hits in this crystal
			{
				energyPerCrystal[i] += pEdep[i]->at(j); //add the energy deposited
		
				std::cout << "(" << px[i]->at(j) << "," << py[i]->at(j) << "," << pz[i]->at(j) << ")" << std::endl;   ///////
			}
			//output to user
			std::cout << "Crystal[" << i <<"]\t\t\t = "  << energyPerCrystal[i] << std::endl;

		}

		}
		std::cout << "Crystals Hit\t\t\t = " << CrystalsHit << std::endl;

		// cerco di capire in che cristallo ha deposto en 


		int index_i;
		int index_j;

		int numero_cristallo; //////

			
		float cry = 3.13;
		float gap = 0.07;

		float cryst_center_x;
		float cryst_center_y;

		float en_depositata=0.;
		float tempo=0;
		float primi_fotoni_vec[10];
		float sigma=0;
		float somma =0;

		// qua ho messo mano
		if(CrystalsHit==1)
		{

			
			for(int i = 0; i < numOfCry ; i++) //crystal by crystal (serve solo per cercare il cristallo in cui ha interagito il gamma)
			{
					if(px[i]->size()) //if the i-th crystal had energy deposition (--> if(0) => false) 
					{
						std::cout << "Il cristallo che è stato colpito è: "<< i << std::endl;
						// converto il numero del cristallo nelle sue coordinate (i,j) col quale viene indicizzato nella matrice
						index_i = i/8;
						index_j = i%8;
						std::cout << "Le coordinate del cristallo sono:\t" << "(" << index_i <<","<<index_j<<")"<<std::endl;
						// converto le coordinate del cristallo nelle (x,y) del suo centro
						index_i = index_i-sqrt(numOfCry)/2;
						index_j = index_j-sqrt(numOfCry)/2;
						std::cout << "Le coordinate shiftate del cristallo sono:\t" << "(" << index_i <<","<<index_j<<")"<<std::endl;
						cryst_center_x = gap/2 + cry/2 + index_i*(gap + cry); 
						cryst_center_y = gap/2 + cry/2 + index_j*(gap + cry); 
						std::cout << "Le coordinate (x,y) del cristallo sono:\t" << "(" << cryst_center_x <<","<< cryst_center_y<<")"<<std::endl;
						std::cout << "I range in cui acquisire fotoni per questo cristallo:\t" << "(" << cryst_center_x-(cry/2) <<","<< cryst_center_x+(cry/2) << ")" <<" x "<< "("  << cryst_center_y-(cry/2) <<","<< cryst_center_y+(cry/2) << ")" << std::endl; 

						for (int a = -1; a<2; a++)
						{
								for (int b = -1; b<2; b++)
								{
										int k =0;
                                                                                std::cout << a << " " <<  b << std::endl;
										for (int m=0; m < pOpticalX->size(); m++)
										{
//                                                                                         std::cout << "m="<< m <<std::endl;
											if ( ((cryst_center_y + b*(gap + cry))-(cry/2)) < (fotone[m].y)  &&  (fotone[m].y) < (cryst_center_y+b*(gap+cry)+(cry/2))  &&   ((cryst_center_x + a*(gap + cry)) - (cry/2) ) < (fotone[m].x)  &&  (fotone[m].x < (cryst_center_x+a*(gap+cry)+(cry/2))) )	
											{
												en_depositata += pPhotonEnergy->at(m);
												if (k<primi_fotoni)
												{
                                                                                                        std::cout << "k=" << k <<std::endl;
													tempo += fotone[m].tempo;
													k++;
													primi_fotoni_vec[k]=fotone[m].tempo;
												}
											}
										}
										double mean = tempo/10;
										if(a==0 && b==0)
										{
											std::cout << "L'mmpc in corrispondenza del cristallo colpito ha raccolto:" << en_depositata <<std::endl;
										}
										
										std::cout << "l'energia depositata nel cristallo di centro (" << cryst_center_x+a*(gap+cry) <<","<< cryst_center_y + b*(gap + cry)<<") è: \t" << en_depositata <<std::endl;
										
										std::cout << "Il tempo medio di arrivo dei primi 10 fotoni su questo mppc è: \t"<< mean <<std::endl;
										
                                                                                
										for (int l=0;l<primi_fotoni;l++){somma += pow(primi_fotoni_vec[l]-mean,2);}
										sigma = sqrt(somma/primi_fotoni);
										std::cout << "con sigma di:\t\t\t\t\t\t\t"<< sigma <<std::endl;
										std::cout << "----------------------------------------------------" << std::endl;
										if (a==0 && b==0)
										{
											timeHisto[i]->Fill(mean);
										}
										en_depositata=0;
										somma=0;
										sigma=0;
										tempo=0;
								}
						}
						
						for (int i=0; i < pOpticalX->size(); i++){en_depositata += pPhotonEnergy->at(i);}
						std::cout << "l'energia depositata dai fotoni ottici dopo l'interazione in questo cristallo è: " << en_depositata <<std::endl;
						

					}
			}
		}

	

		// one example of optical photon related std::vectors
		// get the total number of optical detected 
		// the size of any optical photon related std::vector is M (see above)
		std::cout << "Number of Opticals detected\t = " << pGlobalTime->size() << std::endl;
		
		std::cout << "-----------------------------" << std::endl;
		std::cout << std::endl;

		//progress feedback to the user...
		counter++;
		int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
		if( (perc % 10) == 0 )
		{
		std::cout << "\r";
		std::cout << perc << "% done... ";
		}
  }
  
  
  for(int i = 0; i < numOfCry ; i++)
  {
		
			std::cout << i<<std::endl;
			std::ostringstream p;
			p << "cristallo" << i;
			std::string name(p.str());
			std::cout << name <<std::endl;
			timeCanvas[i]->cd();
			timeHisto[i]->Draw();
			std::cout << name <<std::endl;
			//timeCanvas[i]->Print("name.gif");
			std::cout << name <<std::endl;
			
			std::cout << "prova" << std::endl;
		
  }
  std::cout << "prova1" << std::endl;
//   TFile* provaHistogram = new TFile("prova.root", "new");
  std::cout << "prova2" << std::endl;
  return 0;
}
