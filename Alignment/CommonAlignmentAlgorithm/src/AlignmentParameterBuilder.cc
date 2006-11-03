/** \file AlignableParameterBuilder.cc
 *
 *  $Date: 2006/11/03 11:00:55 $
 *  $Revision: 1.6 $
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CommonDetAlgo/interface/AlgebraicObjects.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignableDet.h"
#include "Alignment/CommonAlignment/interface/AlignableComposite.h"

#include "Alignment/CommonAlignmentParametrization/interface/RigidBodyAlignmentParameters.h"
#include "Alignment/CommonAlignmentParametrization/interface/CompositeRigidBodyAlignmentParameters.h"

// This class's header

#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterBuilder.h"


//__________________________________________________________________________________________________
AlignmentParameterBuilder::AlignmentParameterBuilder( AlignableTracker* alignableTracker )
{

  theAlignableTracker = alignableTracker;

  theTrackerAlignableId = new TrackerAlignableId();
  theOnlyDS=false;
  theOnlySS=false;
  theSelLayers=false;
  theMinLayer=-1;
  theMaxLayer=999;

}

//__________________________________________________________________________________________________
int AlignmentParameterBuilder::addSelections(const edm::ParameterSet &pset)
{
   const std::vector<std::string> selections = 
     pset.getParameter<std::vector<std::string> >("alignableParamSelector");
   bool allOk = true;

   int addedSets = 0;
   for (std::vector<std::string>::const_iterator itSel = selections.begin();
	itSel != selections.end(); ++itSel) {
     std::vector<std::string> decompSel(this->decompose(*itSel, ','));
     if (decompSel.size() < 2) {
       edm::LogError("Alignment") << "@SUB=AlignmentParameterBuilder::addSelections"
				  << "ignoring " << *itSel << "from alignableParamSelector: "
				  << "should have at least 2 ','-separated parts";
       allOk = false;
       continue;
     }

     const std::vector<bool> paramSel(this->decodeParamSel(decompSel[1]));
     if (paramSel.size() != RigidBodyAlignmentParameters::N_PARAM) {
       allOk = false; // Error already from decodeParamSel
       continue;
     }
     if (decompSel.size() > 2) {
       edm::LogWarning("Alignment") << "@SUB=AlignmentParameterBuilder::addSelections"
				    << "r/phi/eta/z-range selection not yet implemented...";
     }
     this->addSelection(decompSel[0], paramSel);
     ++addedSets;
   }

   if (allOk) return addedSets;
   else { // @SUB-syntax is not supported by exception, but anyway useful information...
     throw cms::Exception("BadConfig") <<"@SUB=AlignmentParameterBuilder::addSelections"
				       << ": Problems decoding 'alignableParamSelector'.";
     return -1;
   }
}

//__________________________________________________________________________________________________
std::vector<bool> AlignmentParameterBuilder::decodeParamSel(const std::string &selString) const
{

  if (selString.length() != RigidBodyAlignmentParameters::N_PARAM) {
    edm::LogError("Alignment") <<"@SUB=AlignmentParameterBuilder::decodeSelections"
			       << "selectionString has wrong size " << selString.length()
			       << " instead of " << RigidBodyAlignmentParameters::N_PARAM;
    return std::vector<bool>();
  } else {
    std::vector<bool> result(RigidBodyAlignmentParameters::N_PARAM, false);
    // shifts
    if (selString.substr(0,1)=="1") result[RigidBodyAlignmentParameters::dx] = true;
    if (selString.substr(1,1)=="1") result[RigidBodyAlignmentParameters::dy] = true;
    if (selString.substr(2,1)=="1") result[RigidBodyAlignmentParameters::dz] = true;
    // rotations
    if (selString.substr(3,1)=="1") result[RigidBodyAlignmentParameters::dalpha] = true;
    if (selString.substr(4,1)=="1") result[RigidBodyAlignmentParameters::dbeta] = true;
    if (selString.substr(5,1)=="1") result[RigidBodyAlignmentParameters::dgamma] = true;

    return result;
  }

}


//__________________________________________________________________________________________________
std::vector<std::string> 
AlignmentParameterBuilder::decompose(const std::string &s, std::string::value_type delimiter) const
{

  std::vector<std::string> result;

  std::string::size_type previousPos = 0;
  while (true) {
    const std::string::size_type delimiterPos = s.find(delimiter, previousPos);
    if (delimiterPos == std::string::npos) {
      result.push_back(s.substr(previousPos)); // until end
      break;
    }
    result.push_back(s.substr(previousPos, delimiterPos - previousPos));
    previousPos = delimiterPos + 1;
  }

  return result;
}

//__________________________________________________________________________________________________
void AlignmentParameterBuilder::addSelection(const std::string &name, const std::vector<bool> &sel)
{

  edm::LogWarning("Alignment") << "[AlignmentParameterBuilder] Called for selection >" << name<<"<";

  if      (name == "AllDets")       addAllDets(sel);
  else if (name == "AllRods")       addAllRods(sel);
  else if (name == "AllLayers")     addAllLayers(sel);
  else if (name == "AllComponents") addAllComponents(sel);
  else if (name == "AllAlignables") addAllAlignables(sel);

  // TIB+TOB
  else if (name == "BarrelRods")    add(theAlignableTracker->barrelRods(),sel);
  else if (name == "BarrelDets")    add(theAlignableTracker->barrelGeomDets(),sel);
  else if (name == "BarrelLayers")  add(theAlignableTracker->barrelLayers(),sel);

  else if (name == "BarrelDSRods") {
    theOnlyDS = true;
    add(theAlignableTracker->barrelRods(), sel);
    theOnlyDS = false;
  }
  else if (name == "BarrelSSRods") {
    theOnlySS = true;
    add(theAlignableTracker->barrelRods(), sel);
    theOnlySS = false;
  }

  // PXBarrel
  else if (name == "PixelHalfBarrelDets")
	add(theAlignableTracker->pixelHalfBarrelGeomDets(),sel);
  else if (name == "PixelHalfBarrelLadders") 
	add(theAlignableTracker->pixelHalfBarrelLadders(),sel);
  else if (name == "PixelHalfBarrelLayers")  
	add(theAlignableTracker->pixelHalfBarrelLayers(),sel);

  else if (name == "PixelHalfBarrelLaddersLayers12") {
    theSelLayers=true; theMinLayer=1; theMaxLayer=2;
    add(theAlignableTracker->pixelHalfBarrelLadders(),sel);
    theSelLayers = false;
  }


  // PXEndcap
  else if (name == "PXECDets")      add(theAlignableTracker->pixelEndcapGeomDets(),sel);
  else if (name == "PXECPetals")    add(theAlignableTracker->pixelEndcapPetals(),sel);
  else if (name == "PXECLayers")    add(theAlignableTracker->pixelEndcapLayers(),sel);

  // Pixel Barrel+endcap
  else if (name == "PixelDets") {
    add(theAlignableTracker->pixelHalfBarrelGeomDets(),sel);
    add(theAlignableTracker->pixelEndcapGeomDets(),sel);
  }
  else if (name == "PixelRods") {
    add(theAlignableTracker->pixelHalfBarrelLadders(),sel);
    add(theAlignableTracker->pixelEndcapPetals(),sel);
  }
  else if (name == "PixelLayers") {
    add(theAlignableTracker->pixelHalfBarrelLayers(),sel);
    add(theAlignableTracker->pixelEndcapLayers(),sel);
  }

  // TID
  else if (name == "TIDLayers")     add(theAlignableTracker->TIDLayers(),sel);
  else if (name == "TIDRings")      add(theAlignableTracker->TIDRings(),sel);
  else if (name == "TIDDets")       add(theAlignableTracker->TIDGeomDets(),sel);

  // TEC
  else if (name == "TECDets")       add(theAlignableTracker->endcapGeomDets(),sel); 
  else if (name == "TECPetals")     add(theAlignableTracker->endcapPetals(),sel);
  else if (name == "TECLayers")     add(theAlignableTracker->endcapLayers(),sel);

  // StripEndcap (TID+TEC)
  else if (name == "EndcapDets") {
    add(theAlignableTracker->TIDGeomDets(),sel);
    add(theAlignableTracker->endcapGeomDets(),sel); 
  }
  else if (name == "EndcapPetals") {
    add(theAlignableTracker->TIDRings(),sel);
    add(theAlignableTracker->endcapPetals(),sel);
  }
  else if (name == "EndcapLayers") {
    add(theAlignableTracker->TIDLayers(),sel);
    add(theAlignableTracker->endcapLayers(),sel);
  }

  // Strip Barrel+endcap
  else if (name == "StripDets") {
    add(theAlignableTracker->barrelGeomDets(),sel);
    add(theAlignableTracker->TIDGeomDets(),sel);
    add(theAlignableTracker->endcapGeomDets(),sel); 
  }
  else if (name == "StripRods") {
    add(theAlignableTracker->barrelRods(),sel);
    add(theAlignableTracker->TIDRings(),sel);
    add(theAlignableTracker->endcapPetals(),sel);
  }
  else if (name == "StripLayers") {
    add(theAlignableTracker->barrelLayers(),sel);
    add(theAlignableTracker->TIDLayers(),sel);
    add(theAlignableTracker->endcapLayers(),sel);
  }


  // Custom scenarios

  else if (name == "ScenarioA") {
    std::vector<bool> mysel(6,false);
    // pixel barrel dets x,y,z
    mysel[RigidBodyAlignmentParameters::dx]=true;
    mysel[RigidBodyAlignmentParameters::dy]=true;
    mysel[RigidBodyAlignmentParameters::dz]=true;
    add(theAlignableTracker->pixelHalfBarrelGeomDets(),mysel);
    // strip barrel double sided
    theOnlyDS = true;
    add(theAlignableTracker->barrelRods(),mysel);
    theOnlyDS = false;
    // strip barrel single sided
    mysel[RigidBodyAlignmentParameters::dy]=false;
    theOnlySS = true;
    add(theAlignableTracker->barrelRods(),mysel);
    theOnlySS = false;
  }

  else if (name == "ScenarioB") {
    std::vector<bool> mysel(6,false);
    // pixel barrel ladders x,y,z
    mysel[RigidBodyAlignmentParameters::dx]=true;
    mysel[RigidBodyAlignmentParameters::dy]=true;
    mysel[RigidBodyAlignmentParameters::dz]=true;
    add(theAlignableTracker->pixelHalfBarrelLadders(),mysel);
    // strip barrel layers double sided
    theOnlyDS = true;
    add(theAlignableTracker->barrelLayers(),mysel);
    theOnlyDS = false;
    // strip barrel layers single sided
    mysel[RigidBodyAlignmentParameters::dy]=false;
    theOnlySS = true;
    add(theAlignableTracker->barrelLayers(),mysel);
    theOnlySS = false;
  }


  else if (name == "CustomStripLayers") {
    std::vector<bool> mysel(6,false);
    mysel[RigidBodyAlignmentParameters::dx]=true;
    mysel[RigidBodyAlignmentParameters::dy]=true;
    mysel[RigidBodyAlignmentParameters::dz]=true;
    // strip barrel layers double sided
    theOnlyDS = true;
    add(theAlignableTracker->barrelLayers(),mysel);
    theOnlyDS = false;
    // strip barrel layers single sided
    mysel[RigidBodyAlignmentParameters::dz]=false;
    theOnlySS = true;
    add(theAlignableTracker->barrelLayers(),mysel);
    theOnlySS = false;
    // TID
    mysel[RigidBodyAlignmentParameters::dz]=true;
    add(theAlignableTracker->TIDLayers(),mysel);
    // TEC
    mysel[RigidBodyAlignmentParameters::dz]=false;
    add(theAlignableTracker->endcapLayers(),mysel);
  }

  else if (name == "CustomStripRods") {
    std::vector<bool> mysel(6,false);
    mysel[RigidBodyAlignmentParameters::dx]=true;
    mysel[RigidBodyAlignmentParameters::dy]=true;
    mysel[RigidBodyAlignmentParameters::dz]=true;
    // strip barrel layers double sided
    theOnlyDS = true;
    add(theAlignableTracker->barrelRods(),mysel);
    theOnlyDS = false;
    // strip barrel layers single sided
    mysel[RigidBodyAlignmentParameters::dy]=false;
    theOnlySS = true;
    add(theAlignableTracker->barrelRods(),mysel);
    theOnlySS = false;
    // TID
    mysel[RigidBodyAlignmentParameters::dy]=true;
    add(theAlignableTracker->TIDRings(),mysel);
    // TEC
    mysel[RigidBodyAlignmentParameters::dz]=false;
    add(theAlignableTracker->endcapPetals(),mysel);
  }

  else if (name == "CSA06Selection") {
    std::vector<bool> mysel(6,false);
    mysel[RigidBodyAlignmentParameters::dx]=true;
    mysel[RigidBodyAlignmentParameters::dy]=true;
    mysel[RigidBodyAlignmentParameters::dz]=true;
    mysel[RigidBodyAlignmentParameters::dalpha]=true;
    mysel[RigidBodyAlignmentParameters::dbeta]=true;
    mysel[RigidBodyAlignmentParameters::dgamma]=true;
//  TOB outermost layer (5) kept fixed
    theSelLayers=true; theMinLayer=1; theMaxLayer=5;
//  TOB rods double sided   
    theOnlyDS=true;
    add(theAlignableTracker->outerBarrelRods(),mysel);
    theOnlyDS=false;
// TOB rods single sided   
    mysel[RigidBodyAlignmentParameters::dy]=false;
    mysel[RigidBodyAlignmentParameters::dz]=false;
    theOnlySS=true;
    add(theAlignableTracker->outerBarrelRods(),mysel);
    theOnlySS=false;
    mysel[RigidBodyAlignmentParameters::dy]=true;
    mysel[RigidBodyAlignmentParameters::dz]=true;
//
    theSelLayers=false; 
 // TIB dets double sided   
    theOnlyDS=true;
    add(theAlignableTracker->innerBarrelGeomDets(),mysel);
    theOnlyDS=false;
 // TIB dets single sided   
    mysel[RigidBodyAlignmentParameters::dy]=false;
    mysel[RigidBodyAlignmentParameters::dz]=false;
    theOnlySS=true;
    add(theAlignableTracker->innerBarrelGeomDets(),mysel);
    theOnlySS=false;
  }

  else { // @SUB-syntax is not supported by exception, but anyway useful information... 
    throw cms::Exception("BadConfig") <<"@SUB=AlignmentParameterBuilder::addSelection"
				      << ": Selection '" << name << "' invalid!";
  }
  edm::LogInfo("Warning") << "[AlignmentParameterBuilder] Added " 
			  << theAlignables.size()<< " alignables in total";

}


//__________________________________________________________________________________________________
void AlignmentParameterBuilder::addAllDets(const std::vector<bool> &sel)
{

  add(theAlignableTracker->barrelGeomDets(),sel);          // TIB+TOB
  add(theAlignableTracker->endcapGeomDets(),sel);          // TEC
  add(theAlignableTracker->TIDGeomDets(),sel);             // TID
  add(theAlignableTracker->pixelHalfBarrelGeomDets(),sel); // PixelBarrel
  add(theAlignableTracker->pixelEndcapGeomDets(),sel);     // PixelEndcap

  edm::LogInfo("Alignment") << "Initialized for "
			    << theAlignables.size() << " dets";
}


//__________________________________________________________________________________________________
void AlignmentParameterBuilder::addAllRods(const std::vector<bool> &sel)
{
  add(theAlignableTracker->barrelRods(),sel);
  add(theAlignableTracker->pixelHalfBarrelLadders(),sel);
  add(theAlignableTracker->endcapPetals(),sel);
  add(theAlignableTracker->TIDRings(),sel);
  add(theAlignableTracker->pixelEndcapPetals(),sel);

  edm::LogInfo("Alignment") << "Initialized for "
			    << theAlignables.size() << " rods";
}


//__________________________________________________________________________________________________
void AlignmentParameterBuilder::addAllLayers(const std::vector<bool> &sel)
{
  add(theAlignableTracker->barrelLayers(),sel);
  add(theAlignableTracker->pixelHalfBarrelLayers(),sel);
  add(theAlignableTracker->endcapLayers(),sel);
  add(theAlignableTracker->TIDLayers(),sel);
  add(theAlignableTracker->pixelEndcapLayers(),sel);

  edm::LogInfo("Alignment") << "Initialized for "
			    << theAlignables.size() << " layers";

}


//__________________________________________________________________________________________________
void AlignmentParameterBuilder::addAllComponents(const std::vector<bool> &sel)
{
  add(theAlignableTracker->components(),sel);
  edm::LogInfo("Alignment") << "Initialized for "
			    << theAlignables.size() 
			    << " Components (HalfBarrel/Endcap)";
}


//__________________________________________________________________________________________________
void AlignmentParameterBuilder::addAllAlignables(const std::vector<bool> &sel)
{

  add(theAlignableTracker->barrelGeomDets(),sel);          
  add(theAlignableTracker->endcapGeomDets(),sel);          
  add(theAlignableTracker->TIDGeomDets(),sel);             
  add(theAlignableTracker->pixelHalfBarrelGeomDets(),sel); 
  add(theAlignableTracker->pixelEndcapGeomDets(),sel);     

  add(theAlignableTracker->barrelRods(),sel);
  add(theAlignableTracker->pixelHalfBarrelLadders(),sel);
  add(theAlignableTracker->endcapPetals(),sel);
  add(theAlignableTracker->TIDRings(),sel);
  add(theAlignableTracker->pixelEndcapPetals(),sel);

  add(theAlignableTracker->barrelLayers(),sel);
  add(theAlignableTracker->pixelHalfBarrelLayers(),sel);
  add(theAlignableTracker->endcapLayers(),sel);
  add(theAlignableTracker->TIDLayers(),sel);
  add(theAlignableTracker->pixelEndcapLayers(),sel);

  add(theAlignableTracker->components(),sel);


  edm::LogInfo("Alignment") << "Initialized for "
			    << theAlignables.size() 
			    << " Components (HalfBarrel/Endcap)";

}


//__________________________________________________________________________________________________
void AlignmentParameterBuilder::add(const std::vector<Alignable*> &alignables,
				    const std::vector<bool> &sel)
{

  int num_adu = 0;
  int num_det = 0;
  int num_hlo = 0;

  // loop on Alignable objects
  for ( std::vector<Alignable*>::const_iterator ia=alignables.begin();
        ia!=alignables.end();  ia++ ) {
    Alignable* ali=(*ia);

    // select on single/double sided barrel layers
	std::pair<int,int> tl=theTrackerAlignableId->typeAndLayerFromAlignable( ali );
    int type = tl.first;
    int layer = tl.second;

    bool keep=true;
    if (theOnlySS) // only single sided
      if ( (abs(type)==3 || abs(type)==5) && layer<=2 ) 
		keep=false;

    if (theOnlyDS) // only double sided
      if ( (abs(type)==3 || abs(type)==5) && layer>2 )
		keep=false;

    // reject layers
    if ( theSelLayers && (layer<theMinLayer || layer>theMaxLayer) )  
	  keep=false;


    if (keep) {

	  AlgebraicVector par(6,0);
	  AlgebraicSymMatrix cov(6,0);

	  AlignableDet* alidet = dynamic_cast<AlignableDet*>(ali);
	  if (alidet !=0) { // alignable Det
		RigidBodyAlignmentParameters* dap = 
		  new RigidBodyAlignmentParameters(ali,par,cov,sel);
		ali->setAlignmentParameters(dap);
		num_det++;
	  } else { // higher level object
		CompositeRigidBodyAlignmentParameters* dap = 
		  new CompositeRigidBodyAlignmentParameters(ali,par,cov,sel);
		ali->setAlignmentParameters(dap);
		num_hlo++;
	  }

	  theAlignables.push_back(ali);
	  num_adu++;

    }
  }

  edm::LogWarning("Alignment") << "Added " << num_adu 
			       << " Alignables, of which " << num_det << " are Dets and "
			       << num_hlo << " are higher level.";

}


//__________________________________________________________________________________________________
void AlignmentParameterBuilder::fixAlignables(int n)
{

  if (n<1 || n>3) {
    edm::LogError("BadArgument") << " n = " << n << " is not in [1,3]";
    return;
  }

  std::vector<Alignable*> theNewAlignables;
  int i=0;
  int imax = theAlignables.size();
  for ( std::vector<Alignable*>::const_iterator ia=theAlignables.begin();
        ia!=theAlignables.end();  ia++ ) 
	{
	  i++;
	  if ( n==1 && i>1 ) 
		theNewAlignables.push_back(*ia);
	  else if ( n==2 && i>1 && i<imax ) 
		theNewAlignables.push_back(*ia);
	  else if ( n==3 && i>2 && i<imax) 
		theNewAlignables.push_back(*ia);
	}

  theAlignables = theNewAlignables;

  edm::LogWarning("Alignment") << "removing " << n 
			       << " alignables, so that " << theAlignables.size() 
			       << " alignables left";
  
}

