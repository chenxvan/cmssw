def shifttrackinglayout(i, p, *rows): i["00 Shift/Tracking/" + p] = DQMItem(layout=rows)

shifttrackinglayout(dqmitems, "01 - Tracking ReportSummary",
 [{ 'path': "Tracking/EventInfo/reportSummaryMap",
    'description': " Quality Test results plotted for Tracking parameters : Chi2, TrackRate, #of Hits in Track - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "no" }}])
shifttrackinglayout(dqmitems, "02 - Tracks (pp collisions)",
 [{ 'path': "Tracking/TrackParameters/GeneralProperties/GoodTracks/NumberOfGoodTracks_GenTk",
    'description': "Number of Reconstructed Tracks - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/HitProperties/GoodTracks/GoodTrackNumberOfRecHitsPerTrack_GenTk",
    'description': "Number of RecHits per Track - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/GoodTracks/GoodTrackPt_ImpactPoint_GenTk",
    'description': "Pt of Reconstructed Track - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }}],
 [{ 'path': "Tracking/TrackParameters/GeneralProperties/GoodTracks/GoodTrackChi2oNDF_GenTk",
    'description': "Chi Square per DoF -  <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/GoodTracks/GoodTrackPhi_ImpactPoint_GenTk",
    'description': "Phi distribution of Reconstructed Tracks -  <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/GoodTracks/GoodTrackEta_ImpactPoint_GenTk",
    'description': " Eta distribution of Reconstructed Tracks - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }}])
shifttrackinglayout(dqmitems, "03 - Tracks (Cosmic Tracking)",
 [{ 'path': "Tracking/TrackParameters/GeneralProperties/NumberOfTracks_CKFTk",
    'description': "Number of Reconstructed Tracks - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/HitProperties/NumberOfRecHitsPerTrack_CKFTk",
    'description': "Number of RecHits per Track  - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/TrackPt_CKFTk",
    'description': "Pt of Reconstructed Track  - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }}],
 [{ 'path': "Tracking/TrackParameters/GeneralProperties/Chi2oNDF_CKFTk",
    'description': "Chi Sqare per DoF  -  <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/TrackPhi_CKFTk",
    'description': "Phi distribution of Reconstructed Tracks -  <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/TrackEta_CKFTk",
    'description': " Eta distribution of Reconstructed Tracks - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }}])
shifttrackinglayout(dqmitems, "04 - Number of Seeds (pp collisions)",
 [{ 'path': "Tracking/TrackParameters/TrackBuilding/NumberOfSeeds_initialStepSeeds_iter0",
    'description': "Number of Seed in tracking iteration 0 (no entry: ERROR) - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/SiStripOfflineDQMInstructions>SiStripOfflineDQMInstructions</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/TrackBuilding/NumberOfSeeds_lowPtTripletStepSeeds_iter1",
    'description': "Number of Seed in tracking iteration 1 (no entry: ERROR) - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/SiStripOfflineDQMInstructions>SiStripOfflineDQMInstructions</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/TrackBuilding/NumberOfSeeds_pixelPairStepSeeds_iter2",
    'description': "Number of Seed in tracking iteration 2 (no entry: ERROR) - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/SiStripOfflineDQMInstructions>SiStripOfflineDQMInstructions</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/TrackBuilding/NumberOfSeeds_detachedTripletStepSeeds_iter3",
    'description': "Number of Seed in tracking iteration 3 (no entry: ERROR) - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/SiStripOfflineDQMInstructions>SiStripOfflineDQMInstructions</a> ", 'draw': { 'withref': "yes" }}],
 [{ 'path': "Tracking/TrackParameters/TrackBuilding/NumberOfSeeds_mixedTripletStepSeeds_iter4",
    'description': "Number of Seed in tracking iteration 4 (no entry: ERROR) - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/SiStripOfflineDQMInstructions>SiStripOfflineDQMInstructions</a> ", 'draw': { 'withref': "yes" }},
    { 'path': "Tracking/TrackParameters/TrackBuilding/NumberOfSeeds_pixelLessStepSeeds_iter5",
    'description': "Number of Seed in tracking iteration 5 (no entry: ERROR) - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/SiStripOfflineDQMInstructions>SiStripOfflineDQMInstructions</a> ", 'draw': { 'withref': "yes" }},
    { 'path': "Tracking/TrackParameters/TrackBuilding/NumberOfSeeds_tobTecStepSeeds_iter6",
    'description': "Number of Seed in tracking iteration 6 (no entry: ERROR) - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/SiStripOfflineDQMInstructions>SiStripOfflineDQMInstructions</a> ", 'draw': { 'withref': "yes" }}])
shifttrackinglayout(dqmitems, "05 - Tracks (pp collisions) old layout",
 [{ 'path': "Tracking/TrackParameters/GeneralProperties/NumberOfGoodTracks_GenTk",
    'description': "Number of Reconstructed Tracks - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/HitProperties/GoodTrackNumberOfRecHitsPerTrack_GenTk",
    'description': "Number of RecHits per Track - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/GoodTrackPt_ImpactPoint_GenTk",
    'description': "Pt of Reconstructed Track - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }}],
 [{ 'path': "Tracking/TrackParameters/GeneralProperties/GoodTrackChi2oNDF_GenTk",
    'description': "Chi Square per DoF -  <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/GoodTrackPhi_ImpactPoint_GenTk",
    'description': "Phi distribution of Reconstructed Tracks -  <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }},
  { 'path': "Tracking/TrackParameters/GeneralProperties/GoodTrackEta_ImpactPoint_GenTk",
    'description': " Eta distribution of Reconstructed Tracks - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }}])
shifttrackinglayout(dqmitems, "06 - Primary Vertices",
 [{ 'path': "OfflinePV/offlinePrimaryVertices/vtxNbr",
    'description': "Number reconstructed primary vertices - <a href=https://twiki.cern.ch/twiki/bin/view/CMS/DQMShiftOfflineTracking>DQMShiftOfflineTracking</a> ", 'draw': { 'withref': "yes" }}])