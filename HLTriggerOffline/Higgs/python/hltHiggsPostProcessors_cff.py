import FWCore.ParameterSet.Config as cms

from HLTriggerOffline.Higgs.hltHiggsPostProcessor_cfi import *

# Build the standard strings to the DQM
def efficiency_string(objtype,plot_type,triggerpath):
    # --- IMPORTANT: Add here a elif if you are introduce a new collection
    #                (see EVTColContainer::getTypeString) 
    if objtype == "Mu" :
	objtypeLatex="#mu"
    elif objtype == "Photon": 
	objtypeLatex="#gamma"
    elif objtype == "Ele": 
	objtypeLatex="e"
    elif objtype == "MET" :
	objtypeLatex="MET"
    elif objtype == "PFTau": 
	objtypeLatex="#tau"
    else:
	objtypeLatex=objtype

    numer_description = "# gen %s passed the %s" % (objtypeLatex,triggerpath)
    denom_description = "# gen %s " % (objtypeLatex)

    if plot_type == "TurnOn1":
        title = "pT Turn-On"
        xAxis = "p_{T} of Leading Generated %s (GeV/c)" % (objtype)
        input_type = "gen%sMaxPt1" % (objtype)
    if plot_type == "TurnOn2":
        title = "Next-to-Leading pT Turn-On"
        xAxis = "p_{T} of Next-to-Leading Generated %s (GeV/c)" % (objtype)
        input_type = "gen%sMaxPt2" % (objtype)
    if plot_type == "EffEta":
        title = "#eta Efficiency"
        xAxis = "#eta of Generated %s " % (objtype)
        input_type = "gen%sEta" % (objtype)
    if plot_type == "EffPhi":
        title = "#phi Efficiency"
        xAxis = "#phi of Generated %s " % (objtype)
        input_type = "gen%sPhi" % (objtype)

    yAxis = "%s / %s" % (numer_description, denom_description)
    all_titles = "%s for trigger %s; %s; %s" % (title, triggerpath,
                                        xAxis, yAxis)
    return "Eff_%s_%s '%s' %s_%s %s" % (input_type,triggerpath,
		    all_titles,input_type,triggerpath,input_type)

# Adding the reco objects
def add_reco_strings(strings):
    reco_strings = []
    for entry in strings:
        reco_strings.append(entry
                            .replace("Generated", "Reconstructed")
                            .replace("Gen", "Reco")
                            .replace("gen", "rec"))
    strings.extend(reco_strings)


plot_types = ["TurnOn1", "TurnOn2", "EffEta", "EffPhi"]
#--- IMPORTANT: Update this collection whenever you introduce a new object
#               in the code (from EVTColContainer::getTypeString)
obj_types  = ["Mu","Ele","Photon","MET","PFTau"]
#--- IMPORTANT: Trigger are extracted from the hltHiggsValidator_cfi.py module
triggers = [ ] 
efficiency_strings = []

# Extract the triggers used in the hltHiggsValidator 
from HLTriggerOffline.Higgs.hltHiggsValidator_cfi import hltHiggsValidator as _config
triggers = set([])
for an in _config.analysis:
	s = _config.__getattribute__(an)
	vstr = s.__getattribute__("hltPathsToCheck")
	map(lambda x: triggers.add(x.replace("_v","")),vstr)
triggers = list(triggers)
#------------------------------------------------------------

# Generating the list with all the efficiencies
for type in plot_types:
    for obj in obj_types:
	for trig in triggers:
	    efficiency_strings.append(efficiency_string(obj,type,trig))


#add the summary plots
for an in _config.analysis:
    for trig in triggers:
        efficiency_strings.append("Eff_trueVtxDist_"+an+"_gen_"+trig+" ' Efficiency of "+trig+" vs nb of interactions ; nb events passing each path ' trueVtxDist_"+an+"_gen_"+trig+" trueVtxDist_"+an+"_gen")

add_reco_strings(efficiency_strings)



hltHiggsPostHWW = hltHiggsPostProcessor.clone()
hltHiggsPostHWW.subDirs = ['HLT/Higgs/HWW']
hltHiggsPostHWW.efficiencyProfile = efficiency_strings
hltHiggsEffPostHWW = overallEffPlotter.clone()
hltHiggsEffPostHWW.analysisName = cms.untracked.string('HWW')

hltHiggsPostHZZ = hltHiggsPostProcessor.clone()
hltHiggsPostHZZ.subDirs = ['HLT/Higgs/HZZ']
hltHiggsPostHZZ.efficiencyProfile = efficiency_strings
hltHiggsEffPostHZZ = overallEffPlotter.clone()
hltHiggsEffPostHZZ.analysisName = cms.untracked.string('HZZ')

hltHiggsPostHgg = hltHiggsPostProcessor.clone()
hltHiggsPostHgg.subDirs = ['HLT/Higgs/Hgg']
hltHiggsPostHgg.efficiencyProfile = efficiency_strings
hltHiggsEffPostHgg = overallEffPlotter.clone()
hltHiggsEffPostHgg.analysisName = cms.untracked.string('Hgg')

hltHiggsPostH2tau = hltHiggsPostProcessor.clone()
hltHiggsPostH2tau.subDirs = ['HLT/Higgs/H2tau']
hltHiggsPostH2tau.efficiencyProfile = efficiency_strings
hltHiggsEffPostH2tau = overallEffPlotter.clone()
hltHiggsEffPostH2tau.analysisName = cms.untracked.string('H2tau')

hltHiggsPostHtaunu = hltHiggsPostProcessor.clone()
hltHiggsPostHtaunu.subDirs = ['HLT/Higgs/Htaunu']
hltHiggsPostHtaunu.efficiencyProfile = efficiency_strings
hltHiggsEffPostHtaunu = overallEffPlotter.clone()
hltHiggsEffPostHtaunu.analysisName = cms.untracked.string('Htaunu')


hltHiggsPostProcessors = cms.Sequence(
		hltHiggsPostHWW+
		hltHiggsPostHZZ+
		hltHiggsPostHgg+
		hltHiggsPostHtaunu+
		hltHiggsPostH2tau
)

hltEffPost = cms.Sequence(
    hltHiggsEffPostHWW+
    hltHiggsEffPostHZZ+
    hltHiggsEffPostHgg+
    hltHiggsEffPostH2tau+
    hltHiggsEffPostHtaunu
)
