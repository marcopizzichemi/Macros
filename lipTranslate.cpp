#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>

// ****** Compilation command *************************************************
// g++ -o lipTranslate lipTranslate.cpp `root-config --cflags --glibs`
// ****************************************************************************


int main (int argc, char *argv[])
{
  TString InputfileName = argv[1];
  
  TString outputFile = "out_" +  InputfileName;
  //open the root file
  TFile *fOutput = new TFile(outputFile,"recreate");
  //declare the for output TTree
  TTree *t1;
  t1 = new TTree("adc","adc"); 
  std::stringstream snames[32],stypes[32];
  std::string names[32],types[32];
  
  //variables
  long long int ExtendedTimeTag;
  long long int DeltaTimeTag;
  long long int nEvents = 0;
  Short_t charge[32];
  
  //create the second ttree, t1
  //first 2 branches of the ttree
  t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/L"); 	//absolute time tag of the event
  t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/L"); 			//delta time from previous event
  //branches of the 32 channels data
  for (int i = 0 ; i < 32 ; i++){
    //empty the stringstreams
    snames[i].str(std::string());
    stypes[i].str(std::string());
    charge[i] = 0;
    snames[i] << "ch" << i;
    stypes[i] << "ch" << i << "/S";
    names[i] = snames[i].str();
    types[i] = stypes[i].str();
    t1->Branch(names[i].c_str(),&charge[i],types[i].c_str());
  }
  
  
  TFile *f = new TFile(InputfileName);
  
  
  f->cd();
  
  //declare the ttree
  TTree *tree;
  //get the ttree from tfile
  f->GetObject("lmData",tree);
  
  //variables
  Float_t		step1;
  Float_t		step2;
  UShort_t		mh_n1;
  UShort_t		mh_j1;
  UInt_t		mt_dt1;
  Long64_t		time1;
  UShort_t		channel1;
  Float_t		tot1;
  UShort_t		tac1;
  Double_t		channelIdleTime1;
  Double_t		tacIdleTime1;
  Float_t		tqT1;
  Float_t		tqE1;
  Int_t 		xi1;
  Int_t			yi1;
  Float_t		x1;
  Float_t		y1;
  Float_t		z1;
  UShort_t		mh_n2;
  UShort_t		mh_j2;
  UInt_t		mt_dt2;
  Long64_t		time2;
  UShort_t		channel2;
  Float_t		tot2;
  UShort_t		tac2;
  Double_t		channelIdleTime2;
  Double_t		tacIdleTime2;
  Float_t		tqT2;
  Float_t		tqE2;
  Int_t 		xi2;
  Int_t			yi2;
  Float_t		x2;
  Float_t		y2;
  Float_t		z2;
  
  
  
  //branches
  TBranch		*b_step1;
  TBranch		*b_step2;
  TBranch		*b_mh_n1;
  TBranch		*b_mh_j1;
  TBranch		*b_mt_dt1;
  TBranch		*b_time1;
  TBranch		*b_channel1;
  TBranch		*b_tot1;
  TBranch		*b_tac1;
  TBranch		*b_channelIdleTime1;
  TBranch		*b_tacIdleTime1;
  TBranch		*b_tqT1;
  TBranch		*b_tqE1;
  TBranch		*b_xi1;
  TBranch		*b_yi1;
  TBranch		*b_x1;
  TBranch		*b_y1;
  TBranch		*b_z1;
  TBranch		*b_mh_n2;
  TBranch		*b_mh_j2;
  TBranch		*b_mt_dt2;
  TBranch		*b_time2;
  TBranch		*b_channel2;
  TBranch		*b_tot2;
  TBranch		*b_tac2;
  TBranch		*b_channelIdleTime2;
  TBranch		*b_tacIdleTime2;
  TBranch		*b_tqT2;
  TBranch		*b_tqE2;
  TBranch		*b_xi2;
  TBranch		*b_yi2;
  TBranch		*b_x2;
  TBranch		*b_y2;
  TBranch		*b_z2;
  
  
  //set branch...
  tree->SetBranchAddress("step1",&step1,&b_step1);
  tree->SetBranchAddress("step2",&step2,&b_step2);
  tree->SetBranchAddress("mh_n1",&mh_n1,&b_mh_n1);
  tree->SetBranchAddress("mh_j1",&mh_j1,&b_mh_j1);
  tree->SetBranchAddress("mt_dt1",&mt_dt1,&b_mt_dt1);
  tree->SetBranchAddress("time1",&time1,&b_time1);
  tree->SetBranchAddress("channel1",&channel1,&b_channel1);
  tree->SetBranchAddress("tot1",&tot1,&b_tot1);
  tree->SetBranchAddress("tac1",&tac1,&b_tac1);
  tree->SetBranchAddress("channelIdleTime1",&channelIdleTime1,&b_channelIdleTime1);
  tree->SetBranchAddress("tacIdleTime1",&tacIdleTime1,&b_tacIdleTime1);
  tree->SetBranchAddress("tqT1",&tqT1,&b_tqT1);
  tree->SetBranchAddress("tqE1",&tqE1,&b_tqE1);
  tree->SetBranchAddress("xi1",&xi1,&b_xi1);
  tree->SetBranchAddress("yi1",&yi1,&b_yi1);
  tree->SetBranchAddress("x1",&x1,&b_x1);
  tree->SetBranchAddress("y1",&y1,&b_y1);
  tree->SetBranchAddress("z1",&z1,&b_z1);
  tree->SetBranchAddress("mh_n2",&mh_n2,&b_mh_n2);
  tree->SetBranchAddress("mh_j2",&mh_j2,&b_mh_j2);
  tree->SetBranchAddress("mt_dt2",&mt_dt2,&b_mt_dt2);
  tree->SetBranchAddress("time2",&time2,&b_time2);
  tree->SetBranchAddress("channel2",&channel2,&b_channel2);
  tree->SetBranchAddress("tot2",&tot2,&b_tot2);
  tree->SetBranchAddress("tac2",&tac2,&b_tac2);
  tree->SetBranchAddress("channelIdleTime2",&channelIdleTime2,&b_channelIdleTime2);
  tree->SetBranchAddress("tacIdleTime2",&tacIdleTime2,&b_tacIdleTime2);
  tree->SetBranchAddress("tqT2",&tqT2,&b_tqT2);
  tree->SetBranchAddress("tqE2",&tqE2,&b_tqE2);
  tree->SetBranchAddress("xi2",&xi2,&b_xi2);
  tree->SetBranchAddress("yi2",&yi2,&b_yi2);
  tree->SetBranchAddress("x2",&x2,&b_x2);
  tree->SetBranchAddress("y2",&y2,&b_y2);
  tree->SetBranchAddress("z2",&z2,&b_z2);
  
  Long64_t nentries = tree->GetEntries();
  
  std::cout << "nentries = " << nentries << std::endl;
  
  UShort_t mh_n1_check = 0;
  UShort_t mh_n2_check = 0;
  UShort_t mh_j1_check = 0;
  UShort_t mh_j2_check = 0;
  Long64_t time1_check = 0;
  
  //temporary output for sanity checks
  TFile *fTemp = new TFile("temporary.root","recreate");
  fTemp->cd();
  //declare the ttree
  TTree *tempTree;
  tempTree = new TTree("lmData","lmData");
  
  //variables
  Float_t		temp_step1;
  Float_t		temp_step2;
  UShort_t		temp_mh_n1;
  UShort_t		temp_mh_j1;
  UInt_t		temp_mt_dt1;
  Long64_t		temp_time1;
  UShort_t		temp_channel1;
  Float_t		temp_tot1;
  UShort_t		temp_tac1;
  Double_t		temp_channelIdleTime1;
  Double_t		temp_tacIdleTime1;
  Float_t		temp_tqT1;
  Float_t		temp_tqE1;
  Int_t 		temp_xi1;
  Int_t			temp_yi1;
  Float_t		temp_x1;
  Float_t		temp_y1;
  Float_t		temp_z1;
  UShort_t		temp_mh_n2;
  UShort_t		temp_mh_j2;
  UInt_t		temp_mt_dt2;
  Long64_t		temp_time2;
  UShort_t		temp_channel2;
  Float_t		temp_tot2;
  UShort_t		temp_tac2;
  Double_t		temp_channelIdleTime2;
  Double_t		temp_tacIdleTime2;
  Float_t		temp_tqT2;
  Float_t		temp_tqE2;
  Int_t 		temp_xi2;
  Int_t			temp_yi2;
  Float_t		temp_x2;
  Float_t		temp_y2;
  Float_t		temp_z2;
  
  tempTree->Branch("step1"		,&temp_step1,                         "step1/F"		);
  tempTree->Branch("step2"		,&temp_step2,                         "step2/F"		);
  tempTree->Branch("mh_n1"		,&temp_mh_n1,                         "mh_n1/s"		);
  tempTree->Branch("mh_j1"		,&temp_mh_j1,                         "mh_j1/s"		);
  tempTree->Branch("mt_dt1"		,&temp_mt_dt1,                         "mt_dt1/i"		);
  tempTree->Branch("time1"		,&temp_time1,                         "time1/L"		);
  tempTree->Branch("channel1"		,&temp_channel1,                         "channel1/s"		);
  tempTree->Branch("tot1"		,&temp_tot1,                         "tot1/F"		);
  tempTree->Branch("tac1"		,&temp_tac1,                         "tac1/s"		);
  tempTree->Branch("channelIdleTime1"	,&temp_channelIdleTime1,                         "channelIdleTime1/D"	);
  tempTree->Branch("tacIdleTime1"	,&temp_tacIdleTime1,                         "tacIdleTime1/D"	);
  tempTree->Branch("tqT1"		,&temp_tqT1,                         "tqT1/F"		);
  tempTree->Branch("tqE1"		,&temp_tqE1,                         "tqE1/F"		);
  tempTree->Branch("xi1"		,&temp_xi1,                         "xi1/I"		);
  tempTree->Branch("yi1"		,&temp_yi1,                         "yi1/I"		);
  tempTree->Branch("x1"			,&temp_x1,                         "x1/F"			);
  tempTree->Branch("y1"			,&temp_y1,                         "y1/F"			);
  tempTree->Branch("z1"			,&temp_z1,                         "z1/F"			);
  tempTree->Branch("mh_n2"		,&temp_mh_n2,                         "mh_n2/s"		);
  tempTree->Branch("mh_j2"		,&temp_mh_j2,                         "mh_j2/s"		);
  tempTree->Branch("mt_dt2"		,&temp_mt_dt2,                         "mt_dt2/i"		);
  tempTree->Branch("time2"		,&temp_time2,                         "time2/L"		);
  tempTree->Branch("channel2"		,&temp_channel2,                         "channel2/s"		);
  tempTree->Branch("tot2"		,&temp_tot2,                         "tot2/F"		);
  tempTree->Branch("tac2"		,&temp_tac2,                         "tac2/s"		);
  tempTree->Branch("channelIdleTime2"	,&temp_channelIdleTime2,                         "channelIdleTime2/D"	);
  tempTree->Branch("tacIdleTime2"	,&temp_tacIdleTime2,                         "tacIdleTime2/D"	);
  tempTree->Branch("tqT2"		,&temp_tqT2,                         "tqT2/F"		);
  tempTree->Branch("tqE2"		,&temp_tqE2,                         "tqE2/F"		);
  tempTree->Branch("xi2"		,&temp_xi2,                         "xi2/I"		);
  tempTree->Branch("yi2"		,&temp_yi2,                         "yi2/I"		);
  tempTree->Branch("x2"			,&temp_x2,                         "x2/F"			);
  tempTree->Branch("y2"			,&temp_y2,                         "y2/F"			);
  tempTree->Branch("z2"			,&temp_z2,                         "z2/F"			);
  
  
  //define output variables
  // 1  2  3  4 
  // 5  6  7  8
  // 9  10 11 12
  // 13 14 15 16
  
  //new lip channels
  //57 61 45 41
  //53 49 33 37
  //59 63 47 43
  //55 51 35 39
  
  
  Short_t output_array[32];
  for (int i = 0 ; i < 32 ; i++)
  {
    output_array[i] = 0;
  }
  
  //map from lip channels
  std::map<int,int> lipToDigitizerMap;
  lipToDigitizerMap[55] = 0;
  lipToDigitizerMap[51] = 2;
  lipToDigitizerMap[35] = 4;
  lipToDigitizerMap[39] = 6;
  lipToDigitizerMap[59] = 8;
  lipToDigitizerMap[63] = 10;
  lipToDigitizerMap[47] = 12;
  lipToDigitizerMap[43] = 14;
  lipToDigitizerMap[53] = 16;
  lipToDigitizerMap[49] = 18;
  lipToDigitizerMap[33] = 20;
  lipToDigitizerMap[37] = 22;
  lipToDigitizerMap[57] = 24;
  lipToDigitizerMap[61] = 26;
  lipToDigitizerMap[45] = 28;
  lipToDigitizerMap[41] = 30;
  
  int dummycount = 0;
  
  
  
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    nEvents++;
    if(nEvents % 10000 == 0) 
    { 
      std::cout <<  nEvents << " events\r"; 
    }
    tree->GetEntry(jentry);
//     std::cout 		<< 
//     step1 		<< "\t" << 
//     step2               << "\t" <<
//     mh_n1               << "\t" <<
//     mh_j1               << "\t" <<
//     mt_dt1              << "\t" <<
//     time1               << "\t" <<
//     channel1            << "\t" <<
//     tot1                << "\t" <<
//     tac1                << "\t" <<
//     channelIdleTime1    << "\t" <<
//     tacIdleTime1        << "\t" <<
//     tqT1                << "\t" <<
//     tqE1                << "\t" <<
//     xi1                 << "\t" <<
//     yi1                 << "\t" <<
//     x1                  << "\t" <<
//     y1                  << "\t" <<
//     z1                  << "\t" <<
//     mh_n2               << "\t" <<
//     mh_j2               << "\t" <<
//     mt_dt2              << "\t" <<
//     time2               << "\t" <<
//     channel2            << "\t" <<
//     tot2                << "\t" <<
//     tac2                << "\t" <<
//     channelIdleTime2    << "\t" <<
//     tacIdleTime2        << "\t" <<
//     tqT2                << "\t" <<
//     tqE2                << "\t" <<
//     xi2                 << "\t" <<
//     yi2                 << "\t" <<
//     x2                  << "\t" <<
//     y2                  << "\t" <<
//     z2                  << std::endl;
    
    
    if(time1 == time1_check)
    {
      if(mh_j1 == 0)
      {
	//std::cout << tot2 << "\t" << channel2 << "\t" << lipToDigitizerMap[channel2] << "\t";  
	int index = lipToDigitizerMap[channel2];
	output_array[index] = (Short_t) round(tot2);
      }
      else 
      {
	//we skip events with mh_j1 != 0 because for what concerns module 2 they are just a repetition
      }
    }
    else 
    {
      if(mh_j1 == 0)
      {
	if (time1_check != 0)
	{
	  //write old array
	  //std::cout << std::endl;
	  for (int i = 0 ; i < 32 ; i++)
	  {
	    //std::cout << output_array[i] << "\t";
	    charge[i] = output_array[i];
	    ExtendedTimeTag = time1;
	    DeltaTimeTag = time2;
	  }
	  //std::cout << std::endl;
	  t1->Fill();
	  
	  //reset to 0 output_array
	  for (int i = 0 ; i < 32 ; i++)
	  {
	    output_array[i] = 0;
	    charge[i] = 0;
	  }
	}
	
	//begin new array
	//std::cout << std::endl;
	time1_check = time1; 
	//std::cout << time1_check << "\t" << tot2 << "\t" << channel2 << "\t" << lipToDigitizerMap[channel2] << "\t";
	int index = (int) lipToDigitizerMap[channel2];
	output_array[index] = (int) round(tot2);
	
	//and fill the sanity check ttree
	if(mh_j1 == 0 && mh_j2 == 0)//simple sanity check
	{
	  temp_step1=                   step1;
	  temp_step2=                   step2;
	  temp_mh_n1=                   mh_n1;
	  temp_mh_j1=                   mh_j1;
	  temp_mt_dt1=                  mt_dt1;
	  temp_time1=                   time1;
	  temp_channel1=                channel1;
	  temp_tot1=                    tot1;
	  temp_tac1=                    tac1;
	  temp_channelIdleTime1=        channelIdleTime1;
	  temp_tacIdleTime1=            tacIdleTime1;
	  temp_tqT1=                    tqT1;
	  temp_tqE1=                    tqE1;
	  temp_xi1=                     xi1;
	  temp_yi1=                     yi1;
	  temp_x1=                      x1;
	  temp_y1=                      y1;
	  temp_z1=                      z1;
	  temp_mh_n2=                   mh_n2;
	  temp_mh_j2=                   mh_j2;
	  temp_mt_dt2=                  mt_dt2;
	  temp_time2=                   time2;
	  temp_channel2=                channel2;
	  temp_tot2=                    tot2;
	  temp_tac2=                    tac2;
	  temp_channelIdleTime2=        channelIdleTime2;
	  temp_tacIdleTime2=            tacIdleTime2;
	  temp_tqT2=                    tqT2;
	  temp_tqE2=                    tqE2;
	  temp_xi2=                     xi2;
	  temp_yi2=                     yi2;
	  temp_x2=                      x2;
	  temp_y2=                      y2;
	  temp_z2=                      z2;
	  tempTree->Fill();
	  dummycount++;
	}
      }
    }
  }
  
  
  
  fOutput->cd();
  t1->Write();
  
  
  fTemp->cd();
  tempTree->Write();
  
  
  f->Close();
  fOutput->Close();
  fTemp->Close();
  
  std::cout << std::endl << dummycount << std::endl;
  
  return 0;
}