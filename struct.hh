#include <vector>

struct optPhot
{
  // float energy;
  // float time;
  float PositionX;
  float PositionY;
  float PositionZ;
  float PreMomentumX;
  float PreMomentumY;
  float PreMomentumZ;
  float PostMomentumX;
  float PostMomentumY;
  float PostMomentumZ;
  float GlobalTime;
  int   PhotonType;
  float PhotonEnergy;
  int TrackID;
  int ParentID;
  int   OriginCrystalI;
  int   OriginCrystalJ;
  int   ExitFace;
  //exit face legend:
  // 0 - no exit was found among the next possibilities
  // 1 - exit from front face (face coupled to MPPC)
  // 2 - exit from back face (face opposite to MPPC)
};

struct enDep
{
  int   CrystalI;
  int   CrystalJ;
  int   CrystalID;
  float EnergyDeposited;
  float DepositionTime;
  float DepositionX;
  float DepositionY;
  float DepositionZ;
  int TrackID;
  int ParentID;
  std::string ProcessName;
  std::string ParticleName;
  int PDGEncoding;
};
