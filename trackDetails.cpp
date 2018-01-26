#include "TrackDetails.h"

using namespace std;

TrackDetails::TrackDetails()
{};

TrackDetails::TrackDetails(snemo::datamodel::particle_track track)
{
  this->Initialize(track);
}

void TrackDetails::Initialize(snemo::datamodel::particle_track track)
{
  
  // Populate everything you can
  switch (track.get_charge())
  {
    case snemo::datamodel::particle_track::NEUTRAL:
    {
      particleType_=GAMMA;
      double thisEnergy=0;
      double thisXwallEnergy=0;
      double thisVetoEnergy=0;
      double thisMainWallEnergy=0;
      
      // Store the gamma candidate energies
      double firstHitTime=-1.;
      int firstHitType=0;
      // Store the gamma candidate energies
      for (unsigned int hit=0; hit<track.get_associated_calorimeter_hits().size();++hit)
      {
        
        const snemo::datamodel::calibrated_calorimeter_hit & calo_hit = track.get_associated_calorimeter_hits().at(hit).get();
        double thisHitEnergy=calo_hit.get_energy();
        thisEnergy +=  thisHitEnergy;
        int hitType=calo_hit.get_geom_id().get_type();
        if (hitType==MAINWALL)
        thisMainWallEnergy+= thisHitEnergy;
        else if (hitType==XWALL)
        thisXwallEnergy+= thisHitEnergy;
        else if (hitType==GVETO)
        thisVetoEnergy+= thisHitEnergy;
        else cout<<"WARNING: Unknown calorimeter type "<<hitType<<endl;
        
        // Get the coordinates of the hit with the earliest time
        if (firstHitTime==-1 || calo_hit.get_time()<firstHitTime)
        {
          firstHitTime=calo_hit.get_time();
          // Find out which calo wall
          firstHitType=hitType;
        }
      }

      energy_=thisEnergy;
      firstHitType_=firstHitType;
      // And the fraction of the energy deposited in each wall
      mainwallFraction_=thisMainWallEnergy/thisEnergy;
      xwallFraction_=thisXwallEnergy/thisEnergy;
      vetoFraction_=thisVetoEnergy/thisEnergy;
    } // end case neutral (gammas)
    
    default:
    {
      particleType_=UNKNOWN;
    }
  }//end switch
} // end Initialize

bool TrackDetails::IsGamma()
{
  return (particleType_==GAMMA);
}


double TrackDetails::GetEnergy()
{
  return (energy_);
}

// Fraction of particle's calo energy that is deposited in the main calo wall (France and Italy sides)
double TrackDetails::GetMainwallFraction()
{
  return (mainwallFraction_);
}
// Fraction of particle's calo energy that is deposited in the X-wall (tunnel & mountain ends)
double TrackDetails::GetXwallFraction()
{
  return (xwallFraction_);
}

// Fraction of particle's calo energy that is deposited in the gamma veto (top / bottom)
double TrackDetails::GetVetoFraction()
{
  return (vetoFraction_);
}

// Where did it hit first?
int TrackDetails::GetFirstHitType()
{
  return (firstHitType_);
}
bool TrackDetails::HitMainwall()
{
  return (firstHitType_ == MAINWALL);
}
bool TrackDetails::HitXwall()
{
  return (firstHitType_ == XWALL);
}
bool TrackDetails::HitGammaVeto()
{
  return (firstHitType_ == GVETO);
}

