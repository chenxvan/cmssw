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
  bool applyVertexCut_;
  
  //new
 // const std::string propagatorOpposite_;
  //const std::string estimatorName_;
 // edm::ESHandle<Propagator> thePropagatorOpposite;
  ///edm::ESHandle<Chi2MeasurementEstimatorBase> theEstimator;
  const TrackerTopology*                trackerTopology_;
  const Propagator*                     trackerPropagator_;
  const MeasurementEstimator*           chi2MeasurementEstimator_;
};

SiPixelPhase1TrackEfficiency::SiPixelPhase1TrackEfficiency(const edm::ParameterSet& iConfig) :
  SiPixelPhase1Base(iConfig)//,
  //propagatorOpposite_(iConfig.getParameter<std::string>("PropagatorOpposite")),//new
 // estimatorName_(iConfig.getParameter<std::string>("Chi2MeasurementEstimator"))  //new
{
  tracker_ = consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker")); 
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  vtxToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryvertices"));
  applyVertexCut_=iConfig.getUntrackedParameter<bool>("VertexCut",true);
  trajTrackCollectionToken_ = consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("trajectoryInput"));//new
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

  //new
  edm::Handle<TrajTrackAssociationCollection>  trajTrackCollectionHandle;
  iEvent.getByToken(trajTrackCollectionToken_, trajTrackCollectionHandle);
  
  edm::Handle<MeasurementTrackerEvent> trackerMeas;
  iEvent.getByToken(tracker_, trackerMeas);

  edm::ESHandle<MeasurementTracker> measurementTrackerHandle;
  iSetup.get<CkfComponentsRecord>().get(measurementTrackerHandle);


  //iSetup.get<TrackingComponentsRecord>().get(estimatorName_, theEstimator);
  //iSetup.get<TrackingComponentsRecord>().get(propagatorOpposite_, thePropagatorOpposite);
  

////////////////////////////////////////////////////////////////////////////////////////// 
  std::vector<TrajectoryMeasurement> expTrajMeasurements;
  std::vector<std::pair<int,bool[3]>> eff_pxb1_vector;

  
  for(const auto &pair: *trajTrackCollectionHandle) {
    const edm::Ref<std::vector<Trajectory>> traj = pair.key;
    const reco::TrackRef tracki  = pair.val;



    expTrajMeasurements.clear();
    eff_pxb1_vector.clear();
    if(tracki->pt()<1)continue;
    
   //necessary part here
    
    TrajectoryStateOnSurface tsosPXB2;
    bool valid_layerFrom = false;
    
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
    
    const GeometricSearchTracker * gst_ = trackerMeas->geometricSearchTracker();
    const auto *pxbLayer1_ = gst_->pixelBarrelLayers().front();
    const LayerMeasurements* theLayerMeasurements_ = new LayerMeasurements(*measurementTrackerHandle, *trackerMeas); 
    expTrajMeasurements = theLayerMeasurements_->measurements(*pxbLayer1_, tsosPXB2, *trackerPropagator_, *chi2MeasurementEstimator_);
    
    std::pair<int,bool[3]> eff_map;
    
    for(uint p=0; p<expTrajMeasurements.size();p++){
        TrajectoryMeasurement pxb1TM(expTrajMeasurements[p]);
        const auto& pxb1Hit = pxb1TM.recHit();
        
        bool valid = (pxb1Hit->getType()==TrackingRecHit::valid);
        bool missing = (pxb1Hit->getType()==TrackingRecHit::missing);
        bool inactive = (pxb1Hit->getType()==TrackingRecHit::inactive);
        
        int detid= pxb1Hit->geographicalId();
        
              
               
        std::cout<<" "<<valid<<" "<<missing<<" "<<inactive<<" "<<detid<<std::endl;
        
        bool found_det = false;
        
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
    
    //just checking the hits
    for (unsigned int i_eff=0; i_eff<eff_pxb1_vector.size(); i_eff++){
        std::cout<<"EFFMAP "<<
        eff_pxb1_vector[i_eff].first<<" "<<
        eff_pxb1_vector[i_eff].second[0]<<" valid "<<
        eff_pxb1_vector[i_eff].second[1]<<" missing "<<
        eff_pxb1_vector[i_eff].second[2]<<" inactive zzz"<<std::endl;}
    
    //eff map is filled -> decide what to do for double hits, ie eff_pxb1_vector.size>1 ... if 1 just use MISSING and VALID as usual           
    
    //case eff_pxb1_vector.size()==1 ---> fill
    //case >1 nothing or formula 
    
     expTrajMeasurements.clear();

     //TrajectoryStateOnSurface tsosPXB2;
        for (const auto &tm : traj->measurements()) {
            if (tm.recHit().get() && tm.recHitR().isValid()) {
                DetId where = tm.recHitR().geographicalId();
                int  source_det_ = where.subdetId();
                int  source_layer_ = trackerTopology_ -> pxbLayer(where);
                int source_det2 = trackerTopology_->layer(where);
                /*if (source_det_ != PixelSubdetector::SubDetector::PixelBarrel ||  source_layer_ != 1) {
                    tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                    break;
                }*/
                if(source_det_ == PixelSubdetector::SubDetector::PixelBarrel ) {
 			std::cout << "Pixel Barrel " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl; 
		}
                if(source_det_ == PixelSubdetector::SubDetector::PixelEndcap ) {
                        std::cout << "Pixel Endcap " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TIB ) {
                        std::cout << "TIB " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TOB ) {
                        std::cout << "TOB " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TID ) {
                        std::cout << "TID " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TEC ) {
                        std::cout << "TEC " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if (!(source_det_ == PixelSubdetector::SubDetector::PixelBarrel &&  source_layer_ == 1)) {
                      tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                      //break;
                      std::cout << "starting state on det " << source_det_ << "layer " << source_layer_ << "r = " << tsosPXB2.globalPosition().perp() << "z = " << tsosPXB2.globalPosition().z() << std::endl;
               
		}
 
		if ((source_det_ == PixelSubdetector::SubDetector::PixelBarrel &&  source_layer_ == 1)){

                std::cout << "Global position debug non prorp : x = " << tm.recHit()->globalPosition().x() << 
                     " y = " <<tm.recHit()->globalPosition().y() << " z = " << tm.recHit()->globalPosition().z() <<  std::endl;


		const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(tm.recHit()->hit());
                auto clustref = pixhit->cluster();
		std::cout<<"non propagated cluster position: "<<clustref->x()<<"  "<<clustref->y()<<std::endl; 
                for (const auto &P : clustref->pixels()) {
                    std::cout << "colums "<<P.adc << "  " << P.x << "  " << P.y << std::endl; }}
 

               
            }
        }
        if (!tsosPXB2.isValid()) std::cout << "WARNING: did not find state for PXB2 Hit" << std::endl;
        if (!tsosPXB2.isValid()) continue; // for now
	
	const GeometricSearchTracker * gst = trackerMeas->geometricSearchTracker();
        //const auto & PXBLs = gst->pixelBarrelLayers();
        //        //for (const auto * PXBLayer : PXBLs) { std::cout << "PXB Layer with radius = " << PXBLayer->specificSurface().radius() << std::endl; }
        const auto *pxbLayer1 = gst->pixelBarrelLayers().front();
        auto compDets = pxbLayer1->compatibleDets(tsosPXB2, *trackerPropagator_, *chi2MeasurementEstimator_);

	const LayerMeasurements* theLayerMeasurements = new LayerMeasurements(*measurementTrackerHandle, *trackerMeas); 
	expTrajMeasurements = theLayerMeasurements->measurements(*pxbLayer1, tsosPXB2, *trackerPropagator_, *chi2MeasurementEstimator_);
	for(uint p=0; p<expTrajMeasurements.size();p++){
	TrajectoryMeasurement pxb1TM(expTrajMeasurements[p]);
 	const auto& pxb1Hit = pxb1TM.recHit();
 
	std::cout<<" hit type valid "<<(pxb1Hit->getType()==TrackingRecHit::valid)<<std::endl;
        std::cout<<" hit type missing "<<(pxb1Hit->getType()==TrackingRecHit::missing)<<std::endl;
        std::cout<<" hit type inactive "<<(pxb1Hit->getType()==TrackingRecHit::inactive)<<std::endl;

	if (!pxb1Hit->isValid()) continue;
	std::cout<<"Hit"<<pxb1Hit->isValid()<<std::endl;
	std::cout<<"Hit"<<pxb1Hit->localPosition().x()<<std::endl;
	std::cout<<"Hit"<<pxb1Hit->localPosition().y()<<std::endl;
	std::cout<<" hit type valid "<<(pxb1Hit->getType()==TrackingRecHit::valid)<<std::endl;		
	std::cout<<" hit type missing "<<(pxb1Hit->getType()==TrackingRecHit::missing)<<std::endl;
        std::cout<<" hit type inactive "<<(pxb1Hit->getType()==TrackingRecHit::inactive)<<std::endl;
         printf("\n\n geoid %d ", int(pxb1Hit->geographicalId()) );

        std::cout << "Global position debug hit : x = " << pxb1Hit->globalPosition().x() << " y = " <<pxb1Hit->globalPosition().y() << " z = " << pxb1Hit->globalPosition().z() <<  std::endl;
        

        const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(pxb1Hit->hit());
        auto clustref = pixhit->cluster();
	std::cout<<"propagated cluster position: "<<clustref->x()<<"  "<<clustref->y()<<std::endl;

        for (const auto &P : clustref->pixels()) {
          std::cout << "colums "<<P.adc << "  " << P.x << "  " << P.y << std::endl; }        



	}
        for (const auto & detAndState : compDets) {
	             bool isHitValid   = false;
                     bool isHitMissing = false;
                     bool isHitInactive = false;
                     std::cout<<detAndState.second.localPosition().x()<<"   "<<detAndState.second.localPosition().y()<<std::endl;
                     const auto &mdet = trackerMeas->idToDet(detAndState.first->geographicalId());
                     const auto & allhits = mdet.recHits(detAndState.second);   ////find all the hits in the det
                     std::cout << "Global position debug: x = " << detAndState.second.globalPosition().x() << " y = " <<detAndState.second.globalPosition().y() << " z = " << detAndState.second.globalPosition().z() <<  std::endl;
                     for (const auto & hit : allhits) { 
                    	float distance = std::hypot(std::abs(hit->localPosition().x()-detAndState.second.localPosition().x()), std::abs(hit->localPosition().y()-detAndState.second.localPosition().y()));		  
                        std::cout << "Distance from hit " << distance << std::endl;
		      
                        // Efficiency cut should eta-dependent and maybe with cluster instead of hit
		        if(std::abs(hit->localPosition().x()-detAndState.second.localPosition().x()) < 0.02 && std::abs(hit->localPosition().y()-detAndState.second.localPosition().y()) < 0.03 ){
				std::cout <<  "Hit found!!!" << std::endl;
                                std::cout << hit->localPosition().x() << "  "<<hit->localPosition().y()<<std::endl;
	                        const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit->hit());
				auto clustref = pixhit->cluster();
                                std::cout<<"matched cluster position: "<<clustref->x()<<"  "<<clustref->y()<<std::endl;

                                for (const auto &P : clustref->pixels()) {
                                 std::cout << "colums "<<P.adc << "  " << P.x << "  " << P.y << std::endl; }

                                isHitValid   = hit->getType()==TrackingRecHit::valid;
			}else{
				isHitMissing = true;
			}
		     }
                     if ( (detAndState.first->geographicalId().subdetId() == PixelSubdetector::PixelBarrel && trackerTopology_ -> pxbLayer(detAndState.first->geographicalId()) == 1) ){
    			 if (isHitValid)   {
			        histo[VALID].fill(detAndState.first->geographicalId(), &iEvent);
        			histo[EFFICIENCY].fill(1, detAndState.first->geographicalId(), &iEvent);
      			}
      			if (isHitMissing) {
        			histo[MISSING].fill(detAndState.first->geographicalId(), &iEvent);
        			histo[EFFICIENCY].fill(0, detAndState.first->geographicalId(), &iEvent);
      			}
      			if (isHitInactive)   {
        			histo[INACTIVE].fill(detAndState.first->geographicalId(), &iEvent);
      			}
        	     }
         }

        

}///////////////////////////////////////////////////////////////////////////////////////////////////

  if (!vertices.isValid()) return;

  histo[VERTICES].fill(vertices->size(),DetId(0),&iEvent);

  if (applyVertexCut_ &&  vertices->empty()) return;

  // should be used for weird cuts
  //const auto primaryVertex = vertices->at(0); 

  // get the map
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken( tracksToken_, tracks);
  if (!tracks.isValid()) return;

  for (auto const & track : *tracks) {
   // std::cout<<track.innerDetId()<<std::endl;

    //this cut is needed to be consisten with residuals calculation
    if (applyVertexCut_ && (track.pt() < 0.75 || std::abs( track.dxy(vertices->at(0).position()) ) > 5*track.dxyError())) continue; 

    bool isBpixtrack = false, isFpixtrack = false;
    int nStripHits = 0;
    int nBpixL1Hits = 0;
    int nBpixL2Hits = 0;
    int nBpixL3Hits = 0;
    int nBpixL4Hits = 0;
    int nFpixD1Hits = 0;
    int nFpixD2Hits = 0;
    int nFpixD3Hits = 0;

    // first, look at the full track to see whether it is good
    // auto const & trajParams = track.extra()->trajParams();
    auto hb = track.recHitsBegin();
    for(unsigned int h=0;h<track.recHitsSize();h++){
      
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

    // then, look at each hit
    for(unsigned int h=0;h<track.recHitsSize();h++){
      auto hit = *(hb+h);

      DetId id = hit->geographicalId();
      uint32_t subdetid = (id.subdetId());
      if (   subdetid != PixelSubdetector::PixelBarrel 
          && subdetid != PixelSubdetector::PixelEndcap) continue;

      bool isHitValid   = hit->getType()==TrackingRecHit::valid;
      bool isHitMissing = hit->getType()==TrackingRecHit::missing;
      bool isHitInactive = hit->getType()==TrackingRecHit::inactive;
      
      // Hp cut
      int TRACK_QUALITY_HIGH_PURITY_BIT = 2;
      int TRACK_QUALITY_HIGH_PURITY_MASK = 1 << TRACK_QUALITY_HIGH_PURITY_BIT;
      if(!((track.qualityMask() & TRACK_QUALITY_HIGH_PURITY_MASK) >> TRACK_QUALITY_HIGH_PURITY_BIT)) 
	{
	  isHitValid = false;

	}

      // Pt cut
      float TRACK_PT_CUT_VAL = 1.0f;
      if(!(TRACK_PT_CUT_VAL < track.pt()))
	{
	  isHitValid = false;  
	}

      // Nstrip cut
      int TRACK_NSTRIP_CUT_VAL = 10;
      if(!(TRACK_NSTRIP_CUT_VAL < nStripHits))
	{
	  isHitValid = false;
	}

      // Pixhit cut
      if(subdetid == PixelSubdetector::PixelBarrel)
	{
	  if(trackerTopology_ -> pxbLayer(id) == 1) if(!(
					   (nBpixL2Hits > 0 && nBpixL3Hits > 0 && nBpixL4Hits > 0) ||
					   (nBpixL2Hits > 0 && nBpixL3Hits > 0 && nFpixD1Hits > 0) ||
					   (nBpixL2Hits > 0 && nFpixD1Hits > 0 && nFpixD2Hits > 0) ||
					   (nFpixD1Hits > 0 && nFpixD2Hits > 0 && nFpixD3Hits > 0))) isHitValid = false;
	  if(trackerTopology_ -> pxbLayer(id) == 2) if(!(
					   (nBpixL1Hits > 0 && nBpixL3Hits > 0 && nBpixL4Hits > 0) ||
					   (nBpixL1Hits > 0 && nBpixL3Hits > 0 && nFpixD1Hits > 0) ||
					   (nBpixL1Hits > 0 && nFpixD1Hits > 0 && nFpixD2Hits > 0))) isHitValid = false;
	  if(trackerTopology_ -> pxbLayer(id) == 3) if(!(
					   (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nBpixL4Hits > 0) ||
					   (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nFpixD1Hits > 0))) isHitValid = false;
	  if(trackerTopology_ -> pxbLayer(id) == 4) if(!(
					   (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nBpixL3Hits > 0))) isHitValid = false;
	}
      if(subdetid == PixelSubdetector::PixelBarrel)
	{
	  if(trackerTopology_ -> pxfDisk(id) == 1) if(!(
						    (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nBpixL3Hits > 0) ||
						    (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nFpixD2Hits > 0) ||
						    (nBpixL1Hits > 0 && nFpixD2Hits > 0 && nFpixD3Hits > 0))) isHitValid = false;
	  if(trackerTopology_ -> pxfDisk(id) == 2) if(!(
						    (nBpixL1Hits > 0 && nBpixL2Hits > 0 && nFpixD1Hits > 0) ||
						    (nBpixL1Hits > 0 && nFpixD1Hits > 0 && nFpixD3Hits  > 0))) isHitValid = false;
	  if(trackerTopology_ -> pxfDisk(id) == 3) if(!(
						    (nBpixL1Hits > 0 && nFpixD1Hits > 0 && nFpixD2Hits > 0))) isHitValid = false;
	}


      /*
      const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit);
      const PixelGeomDetUnit* geomdetunit = dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(id) );
      const PixelTopology& topol = geomdetunit->specificTopology();

      // this commented part is useful if one wants ROC level maps of hits, however the local position may fall out of a ROC and the ROC maps will look very strange (with no white cross)
      LocalPoint lp;

      if (pixhit) {
        lp = pixhit->localPosition();
      } else {
        lp = trajParams[h].position();
      }

      MeasurementPoint mp = topol.measurementPosition(lp);
      int row = (int) mp.x();
      int col = (int) mp.y();
      */

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
  histo[VALID   ].executePerEventHarvesting(&iEvent);
  histo[MISSING ].executePerEventHarvesting(&iEvent);
  histo[INACTIVE].executePerEventHarvesting(&iEvent);
}

} // namespace

DEFINE_FWK_MODULE(SiPixelPhase1TrackEfficiency);
