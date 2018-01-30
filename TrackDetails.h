#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TVector3.h"

// - Bayeux
#include "bayeux/dpp/base_module.h"
#include "bayeux/mctools/simulated_data.h"
#include "bayeux/genbb_help/primary_particle.h"
#include "bayeux/genbb_help/primary_event.h"
#include "bayeux/datatools/service_manager.h"
#include "bayeux/geomtools/manager.h"
#include "bayeux/geomtools/geometry_service.h"
#include "bayeux/geomtools/line_3d.h"
#include "bayeux/geomtools/helix_3d.h"
#include "bayeux/geomtools/geomtools.h"

// - Falaise
#include "falaise/snemo/datamodels/calibrated_data.h"
#include "falaise/snemo/datamodels/tracker_clustering_data.h"
#include "falaise/snemo/datamodels/tracker_clustering_solution.h"
#include "falaise/snemo/datamodels/particle_track_data.h"



class TrackDetails{
  TVector3 foilmostVertex_;
  TVector3 direction_;
  bool vertexOnFoil_;
  TVector3 projectedVertex_;
  enum Particle { ELECTRON, GAMMA, ALPHA, UNKNOWN };
  Particle particleType_=UNKNOWN;
  double mainwallFraction_=0;
  double xwallFraction_=0;
  double vetoFraction_=0;
  int firstHitType_=-1;
  double energy_=0;
  int charge_=1;
  bool makesTrack_=false;
  snemo::datamodel::particle_track track_;
  bool hasTrack_=false;
  bool CheckFoilmostVertex(snemo::datamodel::particle_track track);
public:

  // SuperNEMO constants
  int MAINWALL=1302;
  int XWALL=1232;
  int GVETO=1252;
  
  TrackDetails();
  TrackDetails(snemo::datamodel::particle_track track);
  void Initialize(snemo::datamodel::particle_track track);
  bool Initialize();

  bool IsGamma();
  bool IsAlpha();
  bool IsElectron();
  bool IsNegativeElectron();
  bool IsPositron();
  bool MakesTrack();
  
  // Charge
  int GetCharge();
  
  // Where did it hit?
  bool HitMainwall();
  bool HitXwall();
  bool HitGammaVeto();
  int GetFirstHitType();

  // Energies
  double GetEnergy();
  double GetMainwallFraction();
  double GetXwallFraction();
  double GetVetoFraction();

  // Foilmost vertex
  double getFoilmostVertexX();
  double getFoilmostVertexY();
  double getFoilmostVertexZ();
  TVector3 getFoilmostVertex();
  
  // Direction at foilmost end
  double getDirectionX();
  double getDirectionY();
  double getDirectionZ();
  TVector3 getDirection();
  
  // Foil-projected vertex
  double getProjectedVertexX();
  double getProjectedVertexY();
  double getProjectedVertexZ();
  TVector3 getProjectedVertex();
  
  
  
};

