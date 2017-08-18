// compile with
// g++ -o ../build/simRead_doi simRead_doi.cpp `root-config --cflags --glibs` && cp structDictionary.C ../build/
// syntax
// simRead_doi `ls out*`

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
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
#include "TGraph.h"
#include "TNtuple.h"
#include <cmath>        // std::abs

#include "../../Macros/code/struct.hh"


bool compare_by_GlobalTime(const optPhot a, const optPhot b)
{
  return a.GlobalTime < b.GlobalTime;
}

struct evento
{
  int id;
  float z;
  float t[9];
};

int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>"); //needed by ROOT to deal with standard vectors

  //HACK to use the dictionary easily
  std::string fullFileName = "";
  // Code taken from: http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
  std::string path = "";
  pid_t pid = getpid();
  char buf[20] = {0};
  sprintf(buf,"%d",pid);
  std::string _link = "/proc/";
  _link.append( buf );
  _link.append( "/exe");
  char proc[512];
  int ch = readlink(_link.c_str(),proc,512);
  if (ch != -1) {
    proc[ch] = 0;
    path = proc;
    std::string::size_type t = path.find_last_of("/");
    path = path.substr(0,t);
  }
  fullFileName = path + std::string("/");
  //now even worse, assuming the executable is in build and the .C macro in code
  // std::string command = ".L " + fullFileName.substr(0,fullFileName.size()-7) + "/code/structDictionary.C+";
  std::string command = ".L " + fullFileName + "structDictionary.C+";
  // std::cout << fullFileName << std::endl;
  // std::cout << "command " << command << std::endl;
  gROOT->ProcessLine(command.c_str());

  //-------------------
  // speed up IO
  //-------------------
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(NULL);



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
  std::string det_prefix("detector");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    //     leavesName.push_back(leavescopy->At(i)->GetName());
    if (!leavesName[i].compare(0, det_prefix.size(), det_prefix))
    numOfCh++;
    // if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
    // numOfCry++;
  }
  //the string "cry" appears 4 times per crystal..
  // numOfCry = numOfCry / 4;
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;

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
  int  interazioni_avvenute;


  //NEW energy deposition, each gamma 511 event has a std::vector of struct (type enDep) with all the data of each energy deposition
  std::vector<enDep> *energyDeposition = 0;

  //Total number of photons detected in this event
  // for each TTree entry, a simple number saying how many optical photons entered that
  // specific detector, passed the PDE check and where "detected" (i.e. saved)
  Short_t  *detector;
  detector = new Short_t [numOfCh];

  //NEW optical photons. for each gamma 511 event, every optical photon detected is a struct of type optPhot. a std::vector<optPhot> is saved for each gamma 511
  std::vector<optPhot> *photons = 0;
  std::vector<enDep> *EventDepositions = 0;

  std::vector<optPhot> fotoni;


  //------------------------
  // Set Branch Addresses
  //------------------------
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);

  tree->SetBranchAddress("optical",&photons);
  tree->SetBranchAddress("energyDeposition",&energyDeposition);

  for (int i = 0 ; i < numOfCh ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }

  long long int DeltaTimeTag,ExtendedTimeTag;
  Short_t charge[16];


  // correct one
  Double_t xmppc[16]={-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};
  Double_t ymppc[16]={-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};

  double columsum = 0;
  double rowsum = 0;
  double total = 0;
  double floodx = 0;
  double floody = 0;
  double floodz = 0;

  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "***************************************************" << nEntries << std::endl;
  std::cout << "nEntries = " << nEntries << std::endl;
  std::cout << "***************************************************" << nEntries << std::endl;

  long int numero_ottici=0;
  long int foundCandidate = 0;
  long int singleCounter = 0;
  long int doubleCounter = 0;
  long int tripleCounter = 0;
  long int multipleCounter = 0;

  //--------------------------------------------
  //Set here some variables
  //--------------------------------------------
  float z_medio_minimo=20; //imposto uno z minimo alto
  float t_medio_minimo=10; //imposto un t minimo alto
  int cristallo_centrale=28;
  int numero_primi_fotoni=5;
  //---------------------------------------------

  TCanvas *timeCanvas[9];
  TH1F *timeHisto[9];
  TNtuple *ntupla[9];
  TH2F *scatterplot[9];
  TCanvas *scatterCanvas[9];
  int progressivo=0;

  float z_medio;
  float t_medio;
  float delta_time;
  float en_parziale;

  //---------------------------------------------------------//
  // ROBE PRELIMINARI SUI CANVAS E I TEMPI MEDI??
  //---------------------------------------------------------//
  for(int a = -1; a < 2; a++)
  {
    // std::cout<< "prova3" <<std::endl;
    for(int b = -8; b < 9; b=b+8)
    {
      std::ostringstream s,u,v;
      s << "cristallo_" << cristallo_centrale+a+b; //nome per gli istogrammi coi tempi
      u << "ntupla" << cristallo_centrale+a+b; //nome per le ntuple
      v << "t_medio_" << cristallo_centrale+a+b; //nome per gli scatterplot

      std::string histog(s.str());
      std::string ntuple(u.str());
      std::string scatter(v.str());

      timeHisto[progressivo]= new TH1F(histog.c_str(),histog.c_str(), 150, 0,12);
      timeCanvas[progressivo]=new TCanvas(histog.c_str(),histog.c_str(),1200,800);
      ntupla[progressivo]=new TNtuple(ntuple.c_str(),ntuple.c_str(),"centrale+a+b:z_medio:t_medio");

      //Gli assi degli scatterplot variano in base a che mppc stiamo considerando
      if (a==0&&b==0)
      {
        scatterplot[progressivo]= new TH2F(scatter.c_str(),scatter.c_str(), 50, 0,20,20,0.8,1.1);
      }
      else {scatterplot[progressivo]= new TH2F(scatter.c_str(),scatter.c_str(), 50, 0,20,20+abs(a*b),1,2.3+abs(a*b)-abs(a*b/1.5));}

      scatterCanvas[progressivo]=new TCanvas(scatter.c_str(),scatter.c_str(),1200,800);
      progressivo++;

    }
  }


  //---------------------------------------------------------//
  // MAIN LOOP ON EVENTS
  //---------------------------------------------------------//

  std::vector<evento> eventi;


  for(int iEvent = 0; iEvent < nEntries; iEvent++)
  {
    tree->GetEvent(iEvent);
    std::vector<int> crystals;

    // riordino i fotoni
    std::sort(photons->begin(), photons->end(), compare_by_GlobalTime );

    int cry;

    //---------------------------------------------------------//
    // LOOP ON EVENTS IN ORDER TO CONSIDER ONLY SINGL CRYSTAL DEPOSITIONS
    //---------------------------------------------------------//
    for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++) //run on energy depositions and find in how many crystals energy was deposited
    {
      //read the crystal where energy was deposited
      cry = energyDeposition->at(eEvent).CrystalID;
      //loop in the crystals found
      //look for the same id
      bool sameID = false;
      for(int j = 0 ; j < crystals.size(); j++)
      {
        if(crystals[j] == cry) sameID = true;
      }
      if(!sameID) crystals.push_back(cry); //se la deposizione è stata singola tengo l'evento (o il cristallo??)

    }

    int progressivo=0;

    //---------------------------------------------------------//
    // PERFORM ANALYSIS ONLY IN EVENTS IN THE CRYSTAL WE ARE CONSIDERING
    //---------------------------------------------------------//
    if(crystals[0] == cristallo_centrale)
    {
      if(crystals.size() == 1)
      {
        evento evento_in_analisi;
        evento_in_analisi.id=iEvent;

        //calcolo z e t medio per questo evento
        for (int depEvent = 0; depEvent < energyDeposition->size(); depEvent++)
        {
          z_medio+=((energyDeposition->at(depEvent).DepositionZ)+10)*(energyDeposition->at(depEvent).EnergyDeposited);
          en_parziale+=energyDeposition->at(depEvent).EnergyDeposited;
        }

        if(en_parziale > 0.5)
        {
          z_medio=(z_medio/en_parziale);
          if(z_medio<z_medio_minimo){z_medio_minimo=z_medio;}

          evento_in_analisi.z=z_medio;

          for(int a = -1; a < 2; a++)
          {
            for(int b = -8; b < 9; b=b+8)
            {
              int k=0;
              float tempo_parziale=0;
              float energia_parziale=0;
              float zEnergy_parziale=0;
              for(long int oEvent = 0; (oEvent < photons->size())&&(k<numero_primi_fotoni) ; oEvent++) //run on energy depositions and find in how many crystals energy was deposited
              {
                int iDetector = (int) ((photons->at(oEvent).PositionX + 12.7) /  3.2);
                int jDetector = (int) ((photons->at(oEvent).PositionY + 12.7) /  3.2);
                float energy = photons->at(oEvent).PhotonEnergy;

                //per ciascuno dei 9 mppc, calcolo il tempo medio di arrivo dei fotoni
                if (iDetector*8+jDetector == cristallo_centrale+a+b)
                {
                  tempo_parziale += photons->at(oEvent).GlobalTime;
                  k++;
                }
              }

              t_medio=tempo_parziale/numero_primi_fotoni;
              evento_in_analisi.t[progressivo]=t_medio;


              //se sto guardando l'mppc centrale, faccio la roba sul t medio minimo
              if(a==0&&b==0)
              {
                if((t_medio<t_medio_minimo)&&(t_medio!=0))
                {
                  t_medio_minimo=t_medio;
                }
              }
              // delta_time=t_medio-t_medio_minimo;


              // timeHisto[progressivo]->Fill(delta_time);
              // timeHisto[progressivo]->Fill(t_medio-0.749989);
              timeHisto[progressivo]->Fill(t_medio);

              // std::cout<< "provo a riempire l'histo 3"<<std::endl;
              // ntupla[progressivo]->Fill(cristallo_centrale+a+b,z_medio,t_medio-0.749989);
              ntupla[progressivo]->Fill(cristallo_centrale+a+b,z_medio,t_medio);

              // std::cout<< "il tempo totale dei primi 5 è: "<<tempo_parziale <<std::endl;
              // std::cout<< "diviso 5 è: "<<tempo_parziale/numero_primi_fotoni <<std::endl;
              // scatterplot[progressivo]->Fill(z_medio,t_medio-0.749989);
              scatterplot[progressivo]->Fill(z_medio,t_medio);

              progressivo++;

            }
          }

          eventi.push_back(evento_in_analisi);
        }
      }
    }

    numero_ottici=numero_ottici+photons->size();
    counter++;

    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... "<<std::endl;
      //std::cout << counter << std::endl;
    }


  }

  //---------------------------------------------------------//
  // ANALYSIS & OUTPUT
  //---------------------------------------------------------//
  char name[10], title[20];
  TObjArray Histogram_list(0);
  TObjArray NTuple_list(0);
  TObjArray ScatterPlot_list(0);
  progressivo=0;
  for(int a = -1; a < 2; a++)
  {
    // std::cout<< "prova3" <<std::endl;
    for(int b = -8; b < 9; b=b+8)
    {
      sprintf(name,"h%d",cristallo_centrale+b+a);
      sprintf(title,"histo nr:%d",cristallo_centrale+b+a);
      Histogram_list.Add(timeHisto[progressivo]);
      NTuple_list.Add(ntupla[progressivo]);

      // scatterCanvas[progressivo]->cd();
      // scatterplot[progressivo]->Draw("COLZ");

      // Qlist.Add(scatterCanvas[progressivo]);
      ScatterPlot_list.Add(scatterplot[progressivo]);
      progressivo++;
    }
  }


  std::ofstream myfile;


  myfile.open ("output_simRead_doi.txt");

  for(unsigned int i=0; i<eventi.size(); i++)
  myfile << eventi[i].id <<"\t"<< eventi[i].z <<"\t"<< eventi[i].t[0] << "\t"<< eventi[i].t[1]<< "\t"<< eventi[i].t[2]<< "\t"<< eventi[i].t[3]<< "\t"<< eventi[i].t[4]<< "\t"<< eventi[i].t[5]<< "\t"<< eventi[i].t[7]<< "\t"<< eventi[i].t[7]<<
  "\t"<< eventi[i].t[8] <<"\n";

  myfile.close();


  std::string outFile = "FileOut.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  //  numero_interazioni->Write();


  Histogram_list.Write();
  NTuple_list.Write();
  ScatterPlot_list.Write();


  // TCanvas *canvas0 = new TCanvas("canvas0","A Simple Graph Example",200,10,700,500); //canvas per i singoli istogrammi
  // canvas0->cd();
  // lunghezza->Draw();
  // fOut->cd();
  // canvas0->Write();
  //
  // TCanvas *canvas1 = new TCanvas("canvas1","A Simple Graph Example",200,10,700,500); //canvas per i singoli istogrammi
  // canvas1->cd();
  // numero->Draw();
  // fOut->cd();
  // canvas1->Write();
  //


  // averageZvsWi->Write();
  // wivsRatio->Write();
  // wMeasuredvsRatioW->Write();
  // t1->Write();

  //   f1->Close();

  // std::string outFile = "analysis_OUT.root";
  // TFile* fOut = new TFile(outFile.c_str(),"recreate");

  // flood->Write();
  // positions->Write();
  //   f1->Close();

  fOut->Close();

  // std::cout << t_medio_minimo<<std::endl;

  return 0;
}
