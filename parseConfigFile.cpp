#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

struct Params_t {
  std::string saturationFileName;
  bool coincidence;
  bool deepAnalysis;
  bool ChannelOn[32];
  int ChannelModule[32];
  std::string ChannelLabel[32]; 
  int ChannelPosition[32];
} Params;

int parseConfigFile(std::istream& FileHeader, Params_t &Params);
int printConfigFile(Params_t &Params);



int main()
{
  std::fstream configFile("analysis_config.txt");
//   std::string saturationFileName,analysisChoice;
  //std::string ChannelLabel[32];
//   bool deepAnalysis = false;
//   saturationFileName = "";
//   analysisChoice ="no";

  //parse the config file
  parseConfigFile(configFile,Params);
  printConfigFile(Params);
  
  
  
  return 0;
}


int parseConfigFile(std::istream& FileHeader,Params_t &Params)
{
  // defines the strings to look for
//   std::string FileName = "name of data file := ";       
//   std::string NumberFormat = "!number format := ";             
//   std::string SizeX = "!matrix size [1] := ";         
//   std::string SizeY = "!matrix size [2] := ";
//   std::string SizeZ = "!matrix size [3] := "; 
//   std::string SpacingX = "scaling factor (mm/pixel) [1] := ";	 
//   std::string SpacingY = "scaling factor (mm/pixel) [2] := ";	 
//   std::string SpacingZ = "scaling factor (mm/pixel) [3] := ";
  std::stringstream ChIdentifierStream[32];
  std::string ChIdentifierString[32];
  std::string saturation = "!Saturation File Name =";
  std::string analysisString = "!Deep Analysis =";
  std::string coincidenceString = "!Coincidence =";
  //std::string analysisChoice;
  for (int i = 0; i < 32 ; i++)
  {
    ChIdentifierStream[i] << "!Channel [" << i << "] =";
    ChIdentifierString[i] = ChIdentifierStream[i].str(); //set the search channel strings
    Params.ChannelLabel[i] = "VOID"; //Initialize to VOID, then only the specified will be changed
    Params.ChannelOn[i] = false; //Initialize to off all channels 
    Params.deepAnalysis = false;
    Params.coincidence = false;
    Params.saturationFileName = "";
  }
  std::string s1; 
  //looks for the strings defined above, assing the rest of the line to the correct variable      
  while(getline(FileHeader,s1))	
  {
    unsigned int i;
    //look for global parameters
    i = s1.find(saturation); 
    if (i==0)  
    {
      Params.saturationFileName.assign(s1, saturation.size() , s1.size());	
      Params.saturationFileName.erase(remove_if(Params.saturationFileName.begin(), Params.saturationFileName.end(), isspace), Params.saturationFileName.end());
    }
    i = s1.find(analysisString);
    if (i==0)  
    {
      std::string analysis;
      analysis.assign(s1, analysisString.size() , s1.size());
      analysis.erase(remove_if(analysis.begin(), analysis.end(), isspace), analysis.end());
      if(analysis == "yes" |analysis == "Yes" |  analysis == "YES" | analysis == "on" | analysis == "On" | analysis == "ON" | analysis == "1" ){
	Params.deepAnalysis = true;
      } 
    }
    i = s1.find(coincidenceString);
    if (i==0)  
    {
      std::string coincidence;
      coincidence.assign(s1, coincidenceString.size() , s1.size());	
      coincidence.erase(remove_if(coincidence.begin(), coincidence.end(), isspace), coincidence.end());
      if(coincidence == "yes" |coincidence == "Yes" |  coincidence == "YES" | coincidence == "on" | coincidence == "On" | coincidence == "ON" | coincidence == "1" ){
	Params.coincidence = true;
      }
    }
    //channel parameters
    for (int j = 0 ; j < 32 ; j++)
    {
      i = s1.find(ChIdentifierString[j]); //find the info about the channel
      if (i==0) //if the channel is in the config file  
      {
	//first, turn on the channel
	Params.ChannelOn[j] = true;
	//then, tokenize the string
	std::string parametersString,noLeadingSpace;
        parametersString.assign(s1, ChIdentifierString[j].size() , s1.size());
	noLeadingSpace = parametersString.substr(parametersString.find_first_not_of(" "),parametersString.length()-parametersString.find_first_not_of(" "));
	
	//std::cout << "noLeadingSpace " << noLeadingSpace << std::endl;
	std::string::size_type pos = noLeadingSpace.find_first_of(" ,.-\t");
	std::string token = noLeadingSpace.substr(0, pos);
	std::string::size_type pos1 = noLeadingSpace.find_first_of(" ,.-\t",pos+1);
	std::string token1 = noLeadingSpace.substr(pos+1, pos1-1);
	std::string::size_type pos2 = noLeadingSpace.find_first_of(" ,.-\t",pos1+1);
	std::string token2 = noLeadingSpace.substr(pos1+1, pos2-1);
	//std::cout << token << " " << token1 << " " << token2 << std::endl;
	token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
	token1.erase(remove_if(token1.begin(), token1.end(), isspace), token1.end());
	token2.erase(remove_if(token2.begin(), token2.end(), isspace), token2.end());
	//std::cout << token << " " << token1 << " " << token2 << std::endl;
	//finally assign to the Params
	Params.ChannelModule[j] = atoi(token.c_str());
	Params.ChannelLabel[j] = token1;
	Params.ChannelPosition[j] = atoi(token2.c_str());
      }
    }
    
  }
  return 0;
}

int printConfigFile(Params_t &Params)
{
  std::cout << "/***********************************************/" << std::endl;
  std::cout << "|                                               |" << std::endl;
  std::cout << "|          Printing config parameters           |" << std::endl;
  std::cout << "|                                               |" << std::endl;
  std::cout << "/***********************************************/" << std::endl;
  std::cout << "Saturation File =\t" 			<< Params.saturationFileName 		<< std::endl;
  std::cout << "Coincidence \t=\t" 			<< Params.coincidence 		<< std::endl;
  std::cout << "Deep analysis \t=\t" 			<< Params.deepAnalysis 		<< std::endl;
  for (int i = 0 ; i < 32 ; i++){
    std::cout << "Channel [" << i << "] \t=\t"	  	<< Params.ChannelOn[i] 		<< "\t"
                                                        << Params.ChannelModule[i] 		<< "\t"
                                                        << Params.ChannelLabel[i] 		<< "\t"
                                                        << Params.ChannelPosition[i] 		<< std::endl;
    
  }
    
  
  
  
  
  
  return 0;
}