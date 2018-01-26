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
      // There could be multiple hits for a gamma so we need to add them up
      for (unsigned int hit=0; hit<track.get_associated_calorimeter_hits().size();++hit)
      {
        
        const snemo::datamodel::calibrated_calorimeter_hit & calo_hit = track.get_associated_calorimeter_hits().at(hit).get();
        double thisHitEnergy=calo_hit.get_energy();
        // Sum the energies
        thisEnergy +=  thisHitEnergy;
        
        // We want to know what fraction of the energy was deposited in each calo wall
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
          // Find out which calo wall it hit first
          firstHitType=hitType;
        }
      }

      energy_=thisEnergy;
      firstHitType_=firstHitType;
      // And the fraction of the energy deposited in each wall
      mainwallFraction_=thisMainWallEnergy/thisEnergy;
      xwallFraction_=thisXwallEnergy/thisEnergy;
      vetoFraction_=thisVetoEnergy/thisEnergy;
      return;
    } // end case neutral (gammas)
    
    // Any of these will make a track
    case snemo::datamodel::particle_track::POSITIVE:
    case snemo::datamodel::particle_track::NEGATIVE:
    case snemo::datamodel::particle_track::UNDEFINED: // Used for straight tracks
    makesTrack_=true;
    break;
    default:
    {
      particleType_=UNKNOWN; // Nothing we can do here so we are done
      return;
    }
  }//end switch
  
  
  
  // Now we have only charged particles remaining there are a few things we can do:
  // Identify electron candidates
  // Identify alpha candidates
  // Get edgemost inner vertex, regardless of whether they have associated calorimeters etc

  if (!track.has_trajectory())
  {
    particleType_=UNKNOWN;
    return;
  }
  
  const snemo::datamodel::tracker_trajectory & the_trajectory = track.get_trajectory();
  const snemo::datamodel::tracker_cluster & the_cluster = the_trajectory.get_cluster();
    
  // ALPHA candidates are undefined charge particles associated with a delayed hit and no associated hit
  if (track.get_charge()==snemo::datamodel::particle_track::UNDEFINED && !track.has_associated_calorimeter_hits() && the_cluster.is_delayed()>0)
  {
    particleType_=ALPHA;
  }
  // ELECTRON candidates are prompt and have an associated calorimeter hit. No charge requirement as yet
  else if (the_cluster.is_delayed()<=0 && track.has_associated_calorimeter_hits())
  {
    particleType_=ELECTRON;
    charge_=track.get_charge();
  }
  
} // end Initialize


// What particle is it?
bool TrackDetails::IsGamma()
{
  return (particleType_== GAMMA);
}
bool TrackDetails::IsElectron()
{
  return (particleType_== ELECTRON);
}
bool TrackDetails::IsAlpha()
{
  return (particleType_== ALPHA);
}
bool TrackDetails::IsNegativeElectron()
{
  return (particleType_== ELECTRON && charge_==snemo::datamodel::particle_track::POSITIVE);
}
bool TrackDetails::IsPositron()
{
  return (particleType_== ELECTRON && charge_==snemo::datamodel::particle_track::NEGATIVE);
}
int TrackDetails::GetCharge()
{
  return charge_;
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

// Does it make a track? (charged particle)
bool TrackDetails::MakesTrack()
{
  return makesTrack_;
}
