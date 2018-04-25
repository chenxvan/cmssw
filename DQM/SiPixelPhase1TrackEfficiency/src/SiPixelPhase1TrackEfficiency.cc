// -*- C++ -*-
//
// Package:     SiPixelPhase1TrackEfficiency
// Class:       SiPixelPhase1TrackEfficiency
//

// Original Author: Marcel Schneider

#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

/*
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCabling.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingMap.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "DataFormats/SiPixelDetId/interface/PixelFEDChannel.h"
#include "CondFormats/SiPixelObjects/interface/PixelROC.h"
*/


///commnet

namespace {

class SiPixelPhase1TrackEfficiency final : public SiPixelPhase1Base {
  enum {
    VALID,
    MISSING,
    INACTIVE,
    EFFICIENCY,
    VERTICES
  };

  public:
  explicit SiPixelPhase1TrackEfficiency(const edm::ParameterSet& conf);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > clustersToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
  edm::EDGetTokenT<TrajTrackAssociationCollection>  trajTrackCollectionToken_;  
  edm::EDGetTokenT<MeasurementTrackerEvent> tracker_; //new
  //  edm::EDGetTokenT<PixelFEDChannelCollection> pixelFEDChannelCollectionToken_;
  bool applyVertexCut_;

  const TrackerTopology*                trackerTopology_;
  const Propagator*                     trackerPropagator_;
  const MeasurementEstimator*           chi2MeasurementEstimator_;
  //  bool firstEvent_;
  //  const TrackerGeometry* trackerGeometry_ = nullptr;
  //  const SiPixelFedCabling* cablingMap = nullptr;

};

SiPixelPhase1TrackEfficiency::SiPixelPhase1TrackEfficiency(const edm::ParameterSet& iConfig) :
  SiPixelPhase1Base(iConfig)//,
 {
  tracker_ = consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker")); 
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  vtxToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryvertices"));
  applyVertexCut_=iConfig.getUntrackedParameter<bool>("VertexCut",true);
  trajTrackCollectionToken_ = consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("trajectoryInput"));
  clustersToken_=consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("clusters")); 
  //  pixelFEDChannelCollectionToken_ = consumes<PixelFEDChannelCollection>(edm::InputTag("siPixelDigis"));
  //  firstEvent_=true;
}

void SiPixelPhase1TrackEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if( !checktrigger(iEvent,iSetup,DCS) ) return;

  // get geometry
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  assert(tracker.isValid());

  // get primary vertex
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken( vtxToken_, vertices);

  // TrackerTopology for module informations
  edm::ESHandle<TrackerTopology> trackerTopologyHandle;
  iSetup.get<TrackerTopologyRcd>().get(trackerTopologyHandle);
  trackerTopology_ = trackerTopologyHandle.product();

  // Tracker propagator for propagating tracks to other layers
  edm::ESHandle<Propagator> propagatorHandle;
  iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", propagatorHandle);
  std::unique_ptr<Propagator> propagatorUniquePtr(propagatorHandle.product() -> clone());
  trackerPropagator_ = propagatorUniquePtr.get();
  const_cast<Propagator*>(trackerPropagator_) -> setPropagationDirection(oppositeToMomentum);

  // Measurement estimator
  edm::ESHandle<Chi2MeasurementEstimatorBase> chi2MeasurementEstimatorHandle;
  iSetup.get<TrackingComponentsRecord>().get("Chi2", chi2MeasurementEstimatorHandle);
  chi2MeasurementEstimator_ = chi2MeasurementEstimatorHandle.product();
    
  //Tracker
  edm::Handle<MeasurementTrackerEvent> trackerMeas;
  iEvent.getByToken(tracker_, trackerMeas);

  edm::ESHandle<MeasurementTracker> measurementTrackerHandle;
  iSetup.get<CkfComponentsRecord>().get(measurementTrackerHandle);

  //vertices
  if (!vertices.isValid()) return;
  histo[VERTICES].fill(vertices->size(),DetId(0),&iEvent);
  if (applyVertexCut_ &&  vertices->empty()) return;

  // should be used for weird cuts
  //const auto primaryVertex = vertices->at(0); 

  // get the map
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken( tracksToken_, tracks);
  if (!tracks.isValid()) return;
  
  //new
  edm::Handle<TrajTrackAssociationCollection>  trajTrackCollectionHandle;
  iEvent.getByToken(trajTrackCollectionToken_, trajTrackCollectionHandle);
  if (!trajTrackCollectionHandle.isValid()) return;
  
   //Access Pixel Clusters
  edm::Handle< edmNew::DetSetVector<SiPixelCluster> > siPixelClusters;
  iEvent.getByToken(clustersToken_, siPixelClusters);
  if(!siPixelClusters.isValid()) return;  
//   

   edm::ESHandle<PixelClusterParameterEstimator> cpEstimator;
   iSetup.get<TkPixelCPERecord>().get("PixelCPEGeneric", cpEstimator);
   if(!cpEstimator.isValid()) return;
   
   const PixelClusterParameterEstimator &cpe(*cpEstimator);
   const TrackerGeometry *tkgeom=&(*tracker);

   /*
   //Access ROC info
   if( !checktrigger(iEvent,iSetup,DCS) ) return;

   if(firstEvent_){
     edm::ESHandle<TrackerGeometry> tmpTkGeometry;
     iSetup.get<TrackerDigiGeometryRecord>().get(tmpTkGeometry);
     trackerGeometry_=&(*tmpTkGeometry);

     edm::ESHandle<SiPixelFedCablingMap> pixelCabling;
     iSetup.get<SiPixelFedCablingMapRcd>().get(pixelCabling);
     cablingMap = pixelCabling.product();

     firstEvent_=false;
   }

   edm::Handle<edmNew::DetSetVector<PixelFEDChannel> > input;

   iEvent.getByToken(pixelFEDChannelCollectionToken_, input);
   if (!input.isValid()) return;
   */



////////////////////////////////////////////////////////////////////////////////////////// 

  // Hp cut
  int TRACK_QUALITY_HIGH_PURITY_BIT = 2;
  int TRACK_QUALITY_HIGH_PURITY_MASK = 1 << TRACK_QUALITY_HIGH_PURITY_BIT;

  // Pt cut
  float TRACK_PT_CUT_VAL = 1.0f;

  // Nstrip cut
  int TRACK_NSTRIP_CUT_VAL = 10;
  
  //D0
  std::array<float, 4> TRACK_D0_CUT_BARREL_VAL = {{0.01f, 0.02f, 0.02f, 0.02f}};
  float TRACK_D0_CUT_FORWARD_VAL = 0.05f;
  
  //Dz
  float TRACK_DZ_CUT_BARREL_VAL = 0.01f;
  float TRACK_DZ_CUT_FORWARD_VAL = 0.5f;
  
  bool isBpixtrack = false, isFpixtrack = false;
  int nStripHits = 0;
  int nBpixL1Hits = 0;
  int nBpixL2Hits = 0;
  int nBpixL3Hits = 0;
  int nBpixL4Hits = 0;
  int nFpixD1Hits = 0;
  int nFpixD2Hits = 0;
  int nFpixD3Hits = 0;
  bool passcuts = true;
  bool passcuts_hit = true;
  
      
  TrajectoryStateOnSurface tsosPXB2;
  bool valid_layerFrom = false;
  
    
  const GeometricSearchTracker * gst_ = trackerMeas->geometricSearchTracker();
  const auto *pxbLayer1_ = gst_->pixelBarrelLayers().front();
  const LayerMeasurements* theLayerMeasurements_ = new LayerMeasurements(*measurementTrackerHandle, *trackerMeas); 
  
  std::vector<TrajectoryMeasurement> expTrajMeasurements;
  std::vector<std::pair<int,bool[3]>> eff_pxb1_vector;

  
  for(const auto &pair: *trajTrackCollectionHandle) {
    const edm::Ref<std::vector<Trajectory>> traj = pair.key;
    const reco::TrackRef track  = pair.val;

    expTrajMeasurements.clear();
    eff_pxb1_vector.clear();
    //this cut is needed to be consisten with residuals calculation
    if (applyVertexCut_ && (track->pt() < 0.75 || std::abs( track->dxy(vertices->at(0).position()) ) > 5*track->dxyError())) continue; 

    isBpixtrack = false, isFpixtrack = false;
    nStripHits = 0;
    nBpixL1Hits = 0;
    nBpixL2Hits = 0;
    nBpixL3Hits = 0;
    nBpixL4Hits = 0;
    nFpixD1Hits = 0;
    nFpixD2Hits = 0;
    nFpixD3Hits = 0;
    passcuts = true;
    passcuts_hit = true;

    // first, look at the full track to see whether it is good
    //auto const & trajParams = track.extra()->trajParams();
    

    
    auto hb = track->recHitsBegin();
    for(unsigned int h=0;h<track->recHitsSize();h++){
      
      auto hit = *(hb+h);
      if(!hit->isValid()) continue;

      DetId id = hit->geographicalId();
      uint32_t subdetid = (id.subdetId());

      //Check the location of valid hit
      if (subdetid == PixelSubdetector::PixelBarrel && hit->isValid())
	{
	  isBpixtrack = true;
	  if(trackerTopology_ -> pxbLayer(id) == 1) nBpixL1Hits++;
	  if(trackerTopology_ -> pxbLayer(id) == 2) nBpixL2Hits++;
	  if(trackerTopology_ -> pxbLayer(id) == 3) nBpixL3Hits++;
          if(trackerTopology_ -> pxbLayer(id) == 4) nBpixL4Hits++;
	}
      if (subdetid == PixelSubdetector::PixelEndcap && hit->isValid())
	{
	  isFpixtrack = true;
	  if(trackerTopology_ -> pxfDisk(id) == 1) nFpixD1Hits++;
	  if(trackerTopology_ -> pxfDisk(id) == 2) nFpixD2Hits++;
	  if(trackerTopology_ -> pxfDisk(id) == 3) nFpixD3Hits++;
	}

      // count strip hits
      if(subdetid==StripSubdetector::TIB) nStripHits++;
      if(subdetid==StripSubdetector::TOB) nStripHits++;
      if(subdetid==StripSubdetector::TID) nStripHits++;
      if(subdetid==StripSubdetector::TEC) nStripHits++;

      // check that we are in the pixel
      //      if (subdetid == PixelSubdetector::PixelBarrel) isBpixtrack = true;
      //      if (subdetid == PixelSubdetector::PixelEndcap) isFpixtrack = true;
    }

    if (!isBpixtrack && !isFpixtrack) continue;
    
    // Hp cut
    if(!((track->qualityMask() & TRACK_QUALITY_HIGH_PURITY_MASK) >> TRACK_QUALITY_HIGH_PURITY_BIT)) passcuts = false;

    // Pt cut
    if(!(TRACK_PT_CUT_VAL < track->pt()))  passcuts = false;

    // Nstrip cut
    if(!(TRACK_NSTRIP_CUT_VAL < nStripHits))    passcuts = false;
      
    // then, look at each hit
    for(unsigned int h=0;h<track->recHitsSize();h++){
        
      passcuts_hit=true;  
      auto hit = *(hb+h);

      DetId id = hit->geographicalId();
      uint32_t subdetid = (id.subdetId());
      if (   subdetid != PixelSubdetector::PixelBarrel 
          && subdetid != PixelSubdetector::PixelEndcap) continue;

      bool isHitValid   = hit->getType()==TrackingRecHit::valid;
      bool isHitMissing = hit->getType()==TrackingRecHit::missing;
      bool isHitInactive = hit->getType()==TrackingRecHit::inactive;
      
      
      //D0
      if(subdetid == PixelSubdetector::PixelBarrel)
        { if(!((std::abs( track->dxy(vertices->at(0).position()) ) * -1.0) < TRACK_D0_CUT_BARREL_VAL[trackerTopology_ -> pxbLayer(id) -1])) passcuts_hit = false;}
      if(subdetid == PixelSubdetector::PixelEndcap)
        { if(!((std::abs( track->dxy(vertices->at(0).position()) ) * -1.0) < TRACK_D0_CUT_FORWARD_VAL)) passcuts_hit = false;}

      
      //Dz
      if(subdetid == PixelSubdetector::PixelBarrel)
        { if(!(std::abs( track->dz(vertices->at(0).position())) < TRACK_DZ_CUT_BARREL_VAL)) passcuts_hit = false;}
      if(subdetid == PixelSubdetector::PixelEndcap)
        { if(!(std::abs( track->dz(vertices->at(0).position())) < TRACK_DZ_CUT_FORWARD_VAL)) passcuts_hit = false;}

	      

      // Pixhit cut
      if(subdetid == PixelSubdetector::PixelBarrel)
	{
	  if(trackerTopology_ -> pxbLayer(id) == 1) if(!(
					   (nBpixL2Hits > 0 && nBpixL3Hits > 0 && nBpixL4Hits > 0) ||
					   (nBpixL2Hits > 0 && nBpixL3Hits > 0 && nFpixD1Hits > 0) ||
					   (nBpixL2Hits > 0 && nFpixD1Hits > 0 && nFpixD2Hits > 0) ||
					   (nFpixD1Hits > 0 && nFpixD2Hits > 0 && nFpixD3Hits > 0))) passcuts_hit = false;
	  if(trackerTopology_ -> pxbLayer(id) == 2) if(!(
					   (nBpixL1Hits > 0 && nBpixL3Hits > 0 && nBpixL4Hits > 0) ||
					   (nBpixL1Hits > 0 && nBpixL3Hits > 0 && nFpixD1Hits > 0) ||
					   (nBpixL1Hits > 0 && nFpixD1Hits > 0 && nFpixD2Hits > 0))) passcuts_hit = false;
	  if(trackerTopology_ -> pxbLayer(id) == 3) if(!(
					   (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nBpixL4Hits > 0) ||
					   (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nFpixD1Hits > 0))) passcuts_hit = false;
	  if(trackerTopology_ -> pxbLayer(id) == 4) if(!(
							 (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nBpixL3Hits > 0))) passcuts_hit = false;
	}
      if(subdetid == PixelSubdetector::PixelEndcap)
	{
	  if(trackerTopology_ -> pxfDisk(id) == 1) if(!(
						    (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nBpixL3Hits > 0) ||
						    (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nFpixD2Hits > 0) ||
						    (nBpixL1Hits > 0 && nFpixD2Hits > 0 && nFpixD3Hits > 0))) passcuts_hit = false;
	  if(trackerTopology_ -> pxfDisk(id) == 2) if(!(
						    (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nFpixD1Hits > 0) ||
						    (nBpixL1Hits > 0 && nFpixD1Hits > 0 && nFpixD3Hits  > 0))) passcuts_hit = false;
	  if(trackerTopology_ -> pxfDisk(id) == 3) if(!(
							(nBpixL1Hits > 0 && nFpixD1Hits > 0 && nFpixD2Hits > 0))) passcuts_hit = false;
	}
      

      //ROC edge cut
      const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit);
      const PixelGeomDetUnit* geomdetunit = dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(id) );
      const PixelTopology& topol = geomdetunit->specificTopology();

      LocalPoint lp;
      if (pixhit) {
        lp = pixhit->localPosition();
      }
      
      MeasurementPoint mp = topol.measurementPosition(lp);
      int row = (int) mp.x() % 80;
      int col = (int) mp.y() % 52;
      
      int centerrow = 40;
      int centercol = 26;
      
      if (!((col < (centercol + 20)) && (col > (centercol - 20)) && (row < (centerrow + 30)) && (row > (centerrow - 30 )))) passcuts = false;
   

      //smart way to define the row and col
      /*
      for(const auto& disabledOnDetId: *input){
	for(const auto& ch: disabledOnDetId) {
	  sipixelobjects::CablingPathToDetUnit path = {ch.fed, ch.link, 0};
	  for (path.roc=1; path.roc<=(ch.roc_last-ch.roc_first)+1; path.roc++){
	    const sipixelobjects::PixelROC *roc = cablingMap->findItem(path);  
	    assert(roc!=nullptr);
	    assert(roc->rawId()==disabledOnDetId.detId());
	    const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (trackerGeometry_->idToDet(roc->rawId()));
	    PixelTopology const * topology = &(theGeomDet->specificTopology());
	    sipixelobjects::LocalPixel::RocRowCol local = {topology->rowsperroc()/2, topology->colsperroc()/2}; //center of ROC
	    sipixelobjects::GlobalPixel global = roc->toGlobal(sipixelobjects::LocalPixel(local));
	    int rowslow = topology->rowsperroc()/2 - 10;
	    int rowshigh =topology->rowsperroc()/2 + 10;
	    int colslow =topology->colsperroc()/2 - 10;
	    int colshigh =topology->colsperroc()/2 + 10;
	     if (!(global.col < colshigh && global.col > colslow && global.row < rowshigh && global.row > rowslow )) passcuts = false;
	  }
	}
      }
      */



      if (passcuts_hit ==true && passcuts){

	if ( !(subdetid == PixelSubdetector::PixelBarrel && trackerTopology_ -> pxbLayer(id) == 1) ){

	  if (isHitValid)   {
	    histo[VALID].fill(id, &iEvent);
	    histo[EFFICIENCY].fill(1, id, &iEvent);
	  }
	  if (isHitMissing) {
	    histo[MISSING].fill(id, &iEvent);
	    histo[EFFICIENCY].fill(0, id, &iEvent);
	  }
	  if (isHitInactive)   {
	    histo[INACTIVE].fill(id, &iEvent);
	  }
	}

    }
	
    }
    
   //layer 1 specific here
   valid_layerFrom = false;
    
   //propagation only from PXB2 and PXD1, more cuts later
    for (const auto &tm : traj->measurements()) {
        if (tm.recHit().get() && tm.recHitR().isValid()) {
                DetId where = tm.recHitR().geographicalId();
                int  source_det = where.subdetId();
                
                if (source_det == PixelSubdetector::SubDetector::PixelBarrel){
                    int  source_layer = trackerTopology_ -> pxbLayer(where);
                    if (source_layer==2){
                        if (tm.updatedState().isValid()) {tsosPXB2 = tm.updatedState(); valid_layerFrom=true;}
                    }
                }
                
                if (source_det == PixelSubdetector::SubDetector::PixelEndcap){
                    int  source_layer = trackerTopology_ -> pxfDisk(where);
                    if (source_layer==1){
                        if (tm.updatedState().isValid()) {tsosPXB2 = tm.updatedState(); valid_layerFrom=true;}
                                   
                    }
                }    
        }
    } //uodated tsosPXB2 here
    
        
    if (!valid_layerFrom) continue;
    if (!tsosPXB2.isValid()) continue;
    
	//propagation B
    expTrajMeasurements = theLayerMeasurements_->measurements(*pxbLayer1_, tsosPXB2, *trackerPropagator_, *chi2MeasurementEstimator_);
    
    std::pair<int,bool[3]> eff_map;
    
    for(uint p=0; p<expTrajMeasurements.size();p++){
        
        passcuts_hit=true;
        
        TrajectoryMeasurement pxb1TM(expTrajMeasurements[p]);
        const auto& pxb1Hit = pxb1TM.recHit();
        
        
        bool valid = (pxb1Hit->getType()==TrackingRecHit::valid);
        bool missing = (pxb1Hit->getType()==TrackingRecHit::missing);
        bool inactive = (pxb1Hit->getType()==TrackingRecHit::inactive);
        
        int detid= pxb1Hit->geographicalId();
        
        if (detid==0) continue;

	//cuts: exactly the same as for other hits but assuming PXB1
	//D0
        if(!((std::abs( track->dxy(vertices->at(0).position()) ) * -1.0) < TRACK_D0_CUT_BARREL_VAL[trackerTopology_ -> pxbLayer(detid) -1])) passcuts_hit = false;
        //Dz
        if(!(std::abs( track->dz(vertices->at(0).position())) < TRACK_DZ_CUT_BARREL_VAL)) passcuts_hit = false;
        // Pixhit cut
        if(!((nBpixL2Hits > 0 && nBpixL3Hits > 0 && nBpixL4Hits > 0) || (nBpixL2Hits > 0 && nBpixL3Hits > 0 && nFpixD1Hits > 0) ||
	     (nBpixL2Hits > 0 && nFpixD1Hits > 0 && nFpixD2Hits > 0) || (nFpixD1Hits > 0 && nFpixD2Hits > 0 && nFpixD3Hits > 0))) passcuts_hit = false;                    
	
	//ROC edge cut                                                                                                                   
	const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(pxb1Hit->hit());
	const PixelGeomDetUnit* geomdetunit = dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(detid) );
	const PixelTopology& topol = geomdetunit->specificTopology();

	LocalPoint lp;
	if (pixhit) {
	  lp = pixhit->localPosition();
	}

	MeasurementPoint mp = topol.measurementPosition(lp);
	int row = (int) mp.x() % 80;
	int col = (int) mp.y() % 52;

	int centerrow = 40;
	int centercol = 26;

	if (!((col < (centercol + 10)) && (col > (centercol - 10)) && (row < (centerrow + 10)) && (row > (centerrow - 10 )))) passcuts = false;




	auto compDets = pxbLayer1_->compatibleDets(tsosPXB2, *trackerPropagator_, *chi2MeasurementEstimator_);
	for (const auto & detAndState : compDets) {
	  const auto & pXb1_lpos = detAndState.second.localPosition();
	  for (edmNew::DetSetVector<SiPixelCluster>::const_iterator iter_cl=siPixelClusters->begin(); iter_cl!=siPixelClusters->end(); iter_cl++ ){
	    DetId detId(iter_cl->id());
	    if(pxb1Hit->geographicalId().rawId()!=detAndState.first->geographicalId().rawId()) continue;   
	    if(detId.rawId()!=detAndState.first->geographicalId().rawId()) continue;

	    const PixelGeomDetUnit *pixdet=(const PixelGeomDetUnit*) tkgeom->idToDetUnit(detId);
	    edmNew::DetSet<SiPixelCluster>::const_iterator itCluster=iter_cl->begin();
	    float minD[2]; minD[0]=minD[1]=10000.;
	    for( ; itCluster!=iter_cl->end(); ++itCluster){

	      LocalPoint lp(itCluster->x(), itCluster->y(), 0.);
	      PixelClusterParameterEstimator::ReturnType params=cpe.getParameters(*itCluster,*pixdet);
	      lp=std::get<0>(params);

	      float Xdist = abs(lp.x()-pXb1_lpos.x());                                              
	      float Ydist = abs(lp.y()-pXb1_lpos.y());
	      if(Xdist<minD[0]){
		minD[0]=Xdist;
	      }else if(Ydist<minD[1]){
		minD[1]=Ydist;
	      }
	    }
	    if (minD[0] != 10000. && minD[1] != 10000.){
	      if (missing && (minD[0] <= 0.02 && minD[0] <= 0.02)) {
		valid = true; missing = false;

	      }
	    }
	  }

	}
               


        
        bool found_det = false;
        
        if (passcuts && passcuts_hit){
            for (unsigned int i_eff=0; i_eff<eff_pxb1_vector.size(); i_eff++){
                //in case found hit in the same det, take only the valid hit
                if (eff_pxb1_vector[i_eff].first==detid){
                    
                    found_det=true;
                    if (eff_pxb1_vector[i_eff].second[0]==false && valid==true){
                        eff_pxb1_vector[i_eff].second[0]=valid;
                        eff_pxb1_vector[i_eff].second[1]=missing;
                        eff_pxb1_vector[i_eff].second[2]=inactive;
                    } 
                }
                
            }
            
            //if no other hit in det
            if (!found_det) {
            
                eff_map.first=detid;
                eff_map.second[0]=valid;
                eff_map.second[1]=missing;
                eff_map.second[2]=inactive;
            
                eff_pxb1_vector.push_back(eff_map);  
                
            }
            
        }
        
    }
    
    if (eff_pxb1_vector.size() == 1) {   
    
    //eff map is filled -> decide what to do for double hits, ie eff_pxb1_vector.size>1 ... if 1 just use MISSING and VALID as usual      
        if (eff_pxb1_vector[0].second[0]) { 
            histo[VALID].fill(eff_pxb1_vector[0].first, &iEvent); 
            histo[EFFICIENCY].fill(1, eff_pxb1_vector[0].first, &iEvent); 
            
        }
        if (eff_pxb1_vector[0].second[1]) { 
            histo[MISSING].fill(eff_pxb1_vector[0].first, &iEvent); 
            histo[EFFICIENCY].fill(0, eff_pxb1_vector[0].first, &iEvent); 
            
        } 
        if (eff_pxb1_vector[0].second[2]) { 
            histo[INACTIVE].fill(eff_pxb1_vector[0].first, &iEvent); 
            
        }

        }
     
  }

    
  histo[VALID   ].executePerEventHarvesting(&iEvent);
  histo[MISSING ].executePerEventHarvesting(&iEvent);
  histo[INACTIVE].executePerEventHarvesting(&iEvent);
}

} // namespace

DEFINE_FWK_MODULE(SiPixelPhase1TrackEfficiency); 
