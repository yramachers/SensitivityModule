#include "TrackDetails.h"

using namespace std;

TrackDetails::TrackDetails()
{};

TrackDetails::TrackDetails(snemo::datamodel::particle_track track)
{
  
  foilmostVertex_.SetXYZ(-9999,-9999,-9999);
  direction_.SetXYZ(-9999,-9999,-9999);
  projectedVertex_.SetXYZ(-9999,-9999,-9999);
  this->Initialize(track);
}

void TrackDetails::Initialize(snemo::datamodel::particle_track track)
{
  track_=track;
  hasTrack_=true;
  this->Initialize();
}

// Populates all relevant branche based on type of particle (gamma, alpha, electron)
// Returns true if it has identified a particle type and initialized
// Returns false if it can't work out what sort of particle it is

bool TrackDetails::Initialize()
{
  if (!hasTrack_) return false; // You can't get the track details unless there is a track
  charge_=(int)track_.get_charge();
  // Populate everything you can about the track
  switch (charge_)
  {
    case snemo::datamodel::particle_track::NEUTRAL:
    {
      particleType_=GAMMA;
      PopulateCaloHits();
      return true;
    } // end case neutral (gammas)
    
    // Any of these will make a track
    case snemo::datamodel::particle_track::POSITIVE:
    case snemo::datamodel::particle_track::NEGATIVE:
    case snemo::datamodel::particle_track::UNDEFINED: // Used for straight tracks
    makesTrack_=true;
    break;
    default:
    {
      particleType_= UNKNOWN; // Nothing we can do here so we are done
      return false;
    }
  }//end switch
  

  if (!track_.has_trajectory())
  {
    particleType_=UNKNOWN;
    return false;
  }
  
  // Now we have only charged particles remaining there are a few things we can do:
  // Identify electron candidates
  // Identify alpha candidates
  // Get edgemost inner vertex, regardless of whether they have associated calorimeters etc

  const snemo::datamodel::tracker_trajectory & the_trajectory = track_.get_trajectory();
  const snemo::datamodel::tracker_cluster & the_cluster = the_trajectory.get_cluster();
  
  // Number of hits and lengths of track
  trackerHitCount_ = the_cluster.get_number_of_hits(); // Currently a track only contains 1 cluster
  trackLength_ = the_trajectory.get_pattern().get_shape().get_length();
  
  // Get details about the vertex position
  vertexOnFoil_ = SetFoilmostVertex();
  if (SetDirection()) SetProjectedVertex(); // Can't project if no direction!
  
  // ALPHA candidates are undefined charge particles associated with a delayed hit and no associated hit
  if (track_.get_charge()==snemo::datamodel::particle_track::UNDEFINED && !track_.has_associated_calorimeter_hits() && the_cluster.is_delayed()>0)
  {
    particleType_=ALPHA;
    delayTime_ = (the_cluster.get_hit(0).get_delayed_time());
    return true;
  }
  // ELECTRON candidates are prompt and have an associated calorimeter hit. No charge requirement as yet
  else if (the_cluster.is_delayed()<=0 && track_.has_associated_calorimeter_hits())
  {
    particleType_=ELECTRON;
    charge_= track_.get_charge();
    PopulateCaloHits(); // As it has an associated hit, we can calculate the hit fractions
    return true;
  }
  
  
  return false; // Not an alpha or an electron, what could it be?
} // end Initialize

bool TrackDetails::PopulateCaloHits()
{
  if ( !hasTrack_) return false;
  double thisEnergy=0;
  double thisXwallEnergy=0;
  double thisVetoEnergy=0;
  double thisMainWallEnergy=0;
  double firstHitTime=-1.;
  int firstHitType=0;
  // Store the gamma candidate energies
  // There could be multiple hits for a gamma so we need to add them up
  for (unsigned int hit=0; hit<track_.get_associated_calorimeter_hits().size();++hit)
  {
    
    const snemo::datamodel::calibrated_calorimeter_hit & calo_hit = track_.get_associated_calorimeter_hits().at(hit).get();
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
  return true;
}



// Return true if vertex is on the foil
// Populate the inner vertex
bool TrackDetails::SetFoilmostVertex()
{
  if ( !hasTrack_) return false;
  double thisInnerVertex=0; // stores the y coordinate of the vertex closest to the source foil
  double closestX=9999;
  bool hasVertexOnFoil=false;
  if (track_.has_vertices()) // There isn't any time ordering to the vertices so check them all
  {
    for (unsigned int iVertex=0; iVertex<track_.get_vertices().size();++iVertex)
    {
      const geomtools::blur_spot & vertex = track_.get_vertices().at(iVertex).get();
      if (snemo::datamodel::particle_track::vertex_is_on_source_foil(vertex))
      {
        hasVertexOnFoil = true;
      }
      const geomtools::vector_3d & vertexTranslation = vertex.get_placement().get_translation();
      // Get details for the vertex nearest the source foil, which is at x = 0
      if (TMath::Abs(vertexTranslation.x()) < closestX) // this is nearer the foil
      {
        closestX=vertexTranslation.x();
        foilmostVertex_.SetXYZ(vertexTranslation.x(),vertexTranslation.y(),vertexTranslation.z());
      } // end for each vertex
    }
  }
  return hasVertexOnFoil;
}
// Populates the direction_ vector with the direction of the track at the foilmost end
// Returns true if you managed to set it, false if not
bool TrackDetails::SetDirection()
{
  if ( !hasTrack_) return false;
  if (!track_.has_trajectory()) return false; // Can't get the direction without a trajectory!
  const snemo::datamodel::base_trajectory_pattern & the_base_pattern = track_.get_trajectory().get_pattern();
  
  if (the_base_pattern.get_pattern_id()=="line") {
    const geomtools::line_3d & the_shape = (const geomtools::line_3d&)the_base_pattern.get_shape();
    // Find the two ends of the track
    geomtools::vector_3d one_end=the_shape.get_first();
    geomtools::vector_3d the_other_end=the_shape.get_last();
    // which is which?
    geomtools::vector_3d foilmost_end = ((TMath::Abs(one_end.x()) < TMath::Abs(the_other_end.x())) ? one_end: the_other_end);
    geomtools::vector_3d outermost_end = ((TMath::Abs(one_end.x()) >= TMath::Abs(the_other_end.x())) ? one_end: the_other_end);
    geomtools::vector_3d direction = the_shape.get_direction_on_curve(the_shape.get_first()); // Only the first stores the direction for a line track
    int multiplier = (direction.x() * outermost_end.x() > 0)? 1: -1; // If the direction points the wrong way, reverse it
    // This will always point inwards towards the foil
    direction_.SetXYZ(direction.x() * multiplier, direction.y() * multiplier, direction.z() * multiplier);
  } //end line track
  else {
    const geomtools::helix_3d & the_shape = (const geomtools::helix_3d&)the_base_pattern.get_shape();
    // Find the two ends of the track
    geomtools::vector_3d one_end=the_shape.get_first();
    geomtools::vector_3d the_other_end=the_shape.get_last();
    // which is which?
    geomtools::vector_3d foilmost_end = ((TMath::Abs(one_end.x()) < TMath::Abs(the_other_end.x())) ? one_end: the_other_end);
    geomtools::vector_3d outermost_end = ((TMath::Abs(one_end.x()) >= TMath::Abs(the_other_end.x())) ? one_end: the_other_end);
    
    geomtools::vector_3d direction = the_shape.get_direction_on_curve(foilmost_end); // Not the same on a curve
    int multiplier = (direction.x() * outermost_end.x() > 0)? 1: -1; // If the direction points the wrong way, reverse it
    // This will also point in towards the foil. Is that misleading in the case of a track that curves towards the foil and then out again? Not a problem when looking for bb events, but would it be misleading in cases of tracks from the wires?
    direction_.SetXYZ(direction.x() * multiplier, direction.y() * multiplier, direction.z() * multiplier);
  }// end helix track
  return true;
}

// Populate the projectedVertex_ vector with where the vertex would be if it were projected back to the foil
// At the moment this uses a simple linear projection; would be better to project the helix
// Return false if the vertex does not project back to the foil (track would not intersect foil or we don't have enough info)
bool TrackDetails::SetProjectedVertex()
{
  // Check that we have the necessary to do this calculation
  if (GetFoilmostVertexX()==-9999 || GetDirectionX()==-9999 || !hasTrack_) return false;
  
  double scale=foilmostVertex_.X()/direction_.X();
  projectedVertex_=foilmostVertex_ - scale*direction_; // The second term is the extension to the track to project it back with a straight line
  // The direction has been chosen so it will always point outwards from the foil.
  // The calculation should always give a projected X coordinate of 0
  // But if it projects in such a way that the y or z values are outside the detector, we should return false
  if (TMath::Abs(projectedVertex_.Y()) > MAXY || TMath::Abs(projectedVertex_.Z()) > MAXZ)
  {
    return false;
  }
  return true;
}
  
// Getters for the vertex information

// Foilmost vertex
double TrackDetails::GetFoilmostVertexX()
{
  return foilmostVertex_.X();
}
double TrackDetails::GetFoilmostVertexY()
{
  return foilmostVertex_.Y();
}
double TrackDetails::GetFoilmostVertexZ()
{
  return foilmostVertex_.Z();
}
TVector3 TrackDetails::GetFoilmostVertex()
{
  return foilmostVertex_;
}
bool TrackDetails::HasFoilVertex()
{
  return vertexOnFoil_;
}
// Foil-projected vertex
double TrackDetails::GetProjectedVertexX()
{
  return projectedVertex_.X();
}
double TrackDetails::GetProjectedVertexY()
{
  return projectedVertex_.Y();
}
double TrackDetails::GetProjectedVertexZ()
{
  return projectedVertex_.Z();
}
TVector3 TrackDetails::GetProjectedVertex()
{
  return projectedVertex_;
}
// Track direction at the inner vertex
double TrackDetails::GetDirectionX()
{
  return direction_.X();
}
double TrackDetails::GetDirectionY()
{
  return direction_.Y();
}
double TrackDetails::GetDirectionZ()
{
  return direction_.Z();
}
TVector3 TrackDetails::GetDirection()
{
  return direction_;
}



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


// For anything that hits the calo wall
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

// Details of the track
double TrackDetails::GetTrackLength()
{
  return (trackLength_);
}

double TrackDetails::GetDelayTime()
{
  return (delayTime_);
}
int TrackDetails::GetTrackerHitCount()
{
  return trackerHitCount_;
}


// Does it make a track? (charged particle)
bool TrackDetails::MakesTrack()
{
  return makesTrack_;
}
