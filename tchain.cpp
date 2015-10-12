// little program to merge the ttrees and add the new branches needed for "pet" analisys
// compile with
// g++ -o tchain tchain.cpp `root-config --cflags --glibs`

#include <iostream>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include <sstream>
#include <string>



int main(int argc, char** argv)
{
  
  TChain *chain =  new TChain(argv[1]);
  if(argc < 4){ 
    std::cout << std::endl;
    std::cout << "ERROR! Too few parameters" << std::endl;
    std::cout << "USAGE:" << std::endl;
    std::cout << "tchain [name] [mode] [file1.root ... fileN.root] " << std::endl;
    std::cout << "[name] of the ttree in input files" << std::endl;
    std::cout << "[mode] can be:" << std::endl;
    std::cout << "normal = just merge the root files" << std::endl;
    std::cout << "pet    = merge root files and add the relevant branches " << std::endl;
    std::cout << std::endl;
    return 1;
  }
  else {
    std::string mode = argv[2];
    if(mode == "normal"){
      for (int i = 3 ; i < argc ; i++)
      {
	std::cout << "Adding file " << argv[i] << std::endl;
	chain->Add(argv[i]);
      }
      chain->Merge("output_normal.root");
    }
    else{
      if(mode == "pet"){
	for (int i = 3 ; i < argc ; i++)
        {
	  std::cout << "Adding file " << argv[i] << std::endl;
	  chain->Add(argv[i]);
        }
	
	//TFile* fOutput = new TFile("output.root","recreate");
	//fOutput->cd();
	//chain->Merge("output.root");
	
	Float_t xmppc[32]={-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0};
	Float_t ymppc[32]={-4.8,0,-4.8,0,-4.8,0,-4.8,0,-1.6,0,-1.6,0,-1.6,0,-1.6,0,1.6,0,1.6,0,1.6,0,1.6,0,4.8,0,4.8,0,4.8,0,4.8,0};
	//variables for the flood histogram computation
	Float_t columnsum;
	Float_t rowsum;
	Float_t total;
	Float_t posx;
	Float_t posy;
	
	//create the new ttree
	std::stringstream snames,stypes;
	std::string names,types;
	
	TTree *t1 = new TTree("adc","adc"); 
	
	//variables
	long long int ExtendedTimeTag,t1_ExtendedTimeTag;
	long long int DeltaTimeTag,t1_DeltaTimeTag;
	int TriggerChannel;
	int charge[32],t1_charge[32];
	long long int counter = 0;
	float floodx,floody,firstonsecond;
	
	//branches
	TBranch *b_ExtendedTimeTag;
	TBranch *b_DeltaTimeTag;
	TBranch *b_charge[32];
	TBranch *b_TriggerChannel;
	TBranch *b_floodx;
	TBranch *b_floody;
	TBranch *b_firstonsecond;
	//std::stringstream snames[32];
	//std::string names[32];
	chain->SetBranchAddress("ExtendedTimeTag", &ExtendedTimeTag, &b_ExtendedTimeTag);
	chain->SetBranchAddress("DeltaTimeTag", &DeltaTimeTag, &b_DeltaTimeTag);
	for(int i=0; i<32; i++)
	{
	  snames.str(std::string());
	  stypes.str(std::string());
	  snames << "ch" << i;
	  names = snames.str();
	  chain->SetBranchAddress(names.c_str(), &charge[i], &b_charge[i]);
	}
	
	
	
	
	
	//create the second ttree
	//first 2 branches of the ttree
	t1->Branch("ExtendedTimeTag",&t1_ExtendedTimeTag,"ExtendedTimeTag/L"); 	//absolute time tag of the event
	t1->Branch("DeltaTimeTag",&t1_DeltaTimeTag,"DeltaTimeTag/L"); 			//delta time from previous event
	//branches of the 32 channels data
	for (int i = 0 ; i < 32 ; i++){
	  //empty the stringstreams
	  snames.str(std::string());
	  stypes.str(std::string());
	  t1_charge[i] = 0;
	  snames << "ch" << i;
	  stypes << "ch" << i << "/I";
	  names = snames.str();
	  types = stypes.str();
	  t1->Branch(names.c_str(),&t1_charge[i],types.c_str());
	}
	t1->Branch("TriggerChannel",&TriggerChannel,"TriggerChannel/I");		//TriggerChannel = channel where the signal is the highest
	t1->Branch("FloodX",&floodx,"FloodX/F");					//x position in the complete flood histogram
	t1->Branch("FloodY",&floody,"FloodY/F");					//y position in the complete flood histogram
	t1->Branch("FirstOnSecond",&firstonsecond,"FirstOnSecond/F");		//ratio highest signal to second highest 
	
	
	Int_t nevent = chain->GetEntries();
	for (Int_t i=0;i<nevent;i++) {
	  chain->GetEvent(i);              //read complete accepted event in memory
	  double maxCharge = 0;
	  double secondCharge = 0;
	  columnsum=0;
	  rowsum=0;
	  total=0;
	  
	  
	  t1_ExtendedTimeTag = ExtendedTimeTag;
	  t1_DeltaTimeTag = DeltaTimeTag;
	  for (int i = 0 ; i < 32 ; i++){
	    t1_charge[i] = charge[i];
	    if (charge[i] > maxCharge){
	      maxCharge = charge[i];
	      TriggerChannel = i;
	    }
	    total+=charge[i];
	    rowsum += charge[i]*xmppc[i];
	    columnsum += charge[i]*ymppc[i];
	  }
	  //find second highest charge
	  for (int i = 0 ; i < 32 ; i++){
	    if (charge[i] != maxCharge){//TODO here actually I'm not considering the case (unlikely) there are two channels with same value = maxcharge
	      if (charge[i] > secondCharge){
		secondCharge = charge[i];
	      }
	    }
	  }
	  //compute flood x and y
	  floodx=rowsum/total;
	  floody=columnsum/total;
	  //compute first on second ratio
	  firstonsecond = maxCharge / secondCharge;
	  
	  t1->Fill();//fills the tree with the data
	  //counter to give a feedback to the user
	  counter++;
	  if( (counter % 1000) == 0 )
	    std::cout << counter << std::endl;
	  
	}
	
	TFile* fOutput = new TFile("output_pet.root","recreate");
	fOutput->cd();
	t1->Write();
	fOutput->Close();
	
      }
      else{
	std::cout << std::endl;
        std::cout << "ERROR! Invalid choice of [mode]" << std::endl;
        std::cout << "USAGE:" << std::endl;
	std::cout << "tchain [mode] [file1.root ... fileN.root] " << std::endl;
	std::cout << "[mode] can be:" << std::endl;
	std::cout << "normal = just merge the root files" << std::endl;
	std::cout << "pet    = merge root files and add the relevant branches " << std::endl;
	std::cout << std::endl;
      }
    }   
    return 0;  
  }
  
}