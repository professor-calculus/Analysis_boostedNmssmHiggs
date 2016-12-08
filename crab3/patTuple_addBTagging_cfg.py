# import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# switch to un-scheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
# to run in un-scheduled mode uncomment the following lines
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

# input file for local running (requires AOD input)
process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring('/store/mc/RunIISpring16reHLT80/GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v3/20000/4C3007FC-F83F-E611-8B27-0090FAA577A0.root'),
    fileNames = cms.untracked.vstring('file:/hdfs/dpm/phy.bris.ac.uk/home/cms/store/user/taylor/nmssmSignalCascadeV01_13TeV_mH70p0_mSusy1000p0_ratio0p95_splitting1p0/nmssmSignalCascadeV01_13TeV_processMc03_ed01_mH70p0_mSusy1000p0_ratio0p95_splitting1p0/161206_092718/0000/nmssmSignal_AODSIMstep2of2_1.root'),
)

# max number of events when running locally
process.maxEvents.input = 1000
# name of output .root file
process.out.fileName = 'bTagPatTuple.root'
# to suppress the long output at the end of the job
process.options.wantSummary = False

# global tag information
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# what to keep in the output
process.out.outputCommands = cms.untracked.vstring( ['drop *',
                                                     'keep recoGenParticles_genParticles_*_*',
                                                     'keep *_selectedPatJets_tagInfos_*',
                                                     # 'keep *_selectedPatJetsAK4PF_tagInfos_*',
                                                     'keep *_selectedPatJetsAK8PFCHS_tagInfos_*',
                                                     # 'keep *_selectedPatJetsAK8PFCHSSoftDropSubjets_tagInfos_*',
                                                     # 'keep *_selectedPatJets_caloTowers_*',
                                                     # 'keep *_selectedPatJetsAK4PF_caloTowers_*',
                                                     # 'keep *_selectedPatJetsAK8PFCHS_caloTowers_*',
                                                     # 'keep *_selectedPatJetsAK8PFCHSSoftDropSubjets_caloTowers_*',
                                                     'keep *_selectedPatJets__*',
                                                     # 'keep *_selectedPatJetsAK4PF__*',
                                                     'keep *_selectedPatJetsAK8PFCHS__*',
                                                     # 'keep *_selectedPatJetsAK8PFCHSSoftDropSubjets__*',
                                                     # 'keep *_selectedPatJets_genJets_*',
                                                     # 'keep *_selectedPatJetsAK4PF_genJets_*',
                                                     # 'keep *_selectedPatJetsAK8PFCHS_genJets_*',
                                                     # 'keep *_selectedPatJetsAK8PFCHSSoftDropSubjets_genJets_*',
                                                     # 'keep *_selectedPatJets_pfCandidates_*',
                                                     # 'keep *_selectedPatJetsAK4PF_pfCandidates_*',
                                                     # 'keep *_selectedPatJetsAK8PFCHS_pfCandidates_*',
                                                     # 'keep *_selectedPatJetsAK8PFCHSSoftDropSubjets_pfCandidates_*',
                                                     'keep *_patMETs__*',
                                                     # 'keep *_patMETsAK4PF__*',
                                                     # 'keep *_patMETsAK8PFCHS__*',
                                                     # 'keep *_patMETsAK8PFCHSSoftDropSubjets__*',
                                                     'keep *_selectedPatElectrons__*',
                                                     'keep *_selectedPatMuons__*',
                                                     'keep *_selectedPatTaus__*',
                                                     'keep *_selectedPatPhotons__*',
                                                     ] )

# enables us to add different jet collections to the event content
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

# b-tag discriminators for 'standard b tag' (supported with RECO/AOD/MiniAOD)
btagDiscriminators = [
    # <<<<<BTagging>>>>>
    # 'pfJetBProbabilityBJetTags',
    'pfJetProbabilityBJetTags',
    # 'pfPositiveOnlyJetBProbabilityBJetTags',
    # 'pfPositiveOnlyJetProbabilityBJetTags',
    # 'pfNegativeOnlyJetBProbabilityBJetTags',
    # 'pfNegativeOnlyJetProbabilityBJetTags',
    # 'pfTrackCountingHighPurBJetTags',
    # 'pfTrackCountingHighEffBJetTags',
    # 'pfNegativeTrackCountingHighPurBJetTags',
    # 'pfNegativeTrackCountingHighEffBJetTags',
    # 'pfSimpleSecondaryVertexHighEffBJetTags',
    # 'pfSimpleSecondaryVertexHighPurBJetTags',
    # 'pfNegativeSimpleSecondaryVertexHighEffBJetTags',
    # 'pfNegativeSimpleSecondaryVertexHighPurBJetTags',
    # 'pfSimpleInclusiveSecondaryVertexHighEffBJetTags',
    # 'pfSimpleInclusiveSecondaryVertexHighPurBJetTags',
    # 'pfNegativeSimpleInclusiveSecondaryVertexHighEffBJetTags',
    # 'pfNegativeSimpleInclusiveSecondaryVertexHighPurBJetTags',
    # 'pfCombinedSecondaryVertexV2BJetTags',
    # 'pfPositiveCombinedSecondaryVertexV2BJetTags',
    # 'pfNegativeCombinedSecondaryVertexV2BJetTags',
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    # 'pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags',
    # 'pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags',
    # 'pfGhostTrackBJetTags',
    # 'softPFMuonBJetTags',
    # 'softPFMuonByPtBJetTags',
    # 'softPFMuonByIP3dBJetTags',
    # 'softPFMuonByIP2dBJetTags',
    # 'positiveSoftPFMuonBJetTags',
    # 'positiveSoftPFMuonByPtBJetTags',
    # 'positiveSoftPFMuonByIP3dBJetTags',
    # 'positiveSoftPFMuonByIP2dBJetTags',
    # 'negativeSoftPFMuonBJetTags',
    # 'negativeSoftPFMuonByPtBJetTags',
    # 'negativeSoftPFMuonByIP3dBJetTags',
    # 'negativeSoftPFMuonByIP2dBJetTags',
    # 'softPFElectronBJetTags',
    # 'softPFElectronByPtBJetTags',
    # 'softPFElectronByIP3dBJetTags',
    # 'softPFElectronByIP2dBJetTags',
    # 'positiveSoftPFElectronBJetTags',
    # 'positiveSoftPFElectronByPtBJetTags',
    # 'positiveSoftPFElectronByIP3dBJetTags',
    # 'positiveSoftPFElectronByIP2dBJetTags',
    # 'negativeSoftPFElectronBJetTags',
    # 'negativeSoftPFElectronByPtBJetTags',
    # 'negativeSoftPFElectronByIP3dBJetTags',
    # 'negativeSoftPFElectronByIP2dBJetTags',
    'pfCombinedMVAV2BJetTags',
    # 'pfNegativeCombinedMVAV2BJetTags',
    # 'pfPositiveCombinedMVAV2BJetTags',
    # <<<<<CTagging>>>>>
    # 'pfCombinedCvsLJetTags',
    # 'pfCombinedCvsBJetTags',
    # <<<<<ChargeTagging>>>>>
    # 'pfChargeBJetTags',
]

# adds ak4PFJets with new b-tags to your PAT output
addJetCollection(
   process,
   labelName = 'AK4PF',
   jetSource = cms.InputTag('ak4PFJets'),
   jetCorrections = ('AK4PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
   btagDiscriminators = btagDiscriminators
)
process.patJetsAK4PF.addTagInfos = True

# adds ak8PFJetsCHS with double b-tags to your PAT output
addJetCollection(
   process,
   labelName = 'AK8PFCHS',
   jetSource = cms.InputTag('ak8PFJetsCHS'),
   jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
   algo = 'AK',
   rParam = 0.8,
   btagDiscriminators = ['pfBoostedDoubleSecondaryVertexAK8BJetTags']
)
process.patJetsAK8PFCHS.addTagInfos = True

# adds subjets of ak8PFJetsCHSSoftDrop with new b-tags to your PAT output
addJetCollection(
   process,
   labelName = 'AK8PFCHSSoftDropSubjets',
   jetSource = cms.InputTag('ak8PFJetsCHSSoftDrop','SubJets'),
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'), # Using AK4 JECs for subjets which might not be completely appropriate
   algo = 'AK',  # needed for subjet flavor clustering
   rParam = 0.8, # needed for subjet flavor clustering
   btagDiscriminators = btagDiscriminators,
   explicitJTA = True,  # needed for subjet b tagging
   svClustering = True, # needed for subjet b tagging
   fatJets = cms.InputTag("ak8PFJetsCHS"),               # needed for subjet flavor clustering
   groomedFatJets = cms.InputTag("ak8PFJetsCHSSoftDrop") # needed for subjet flavor clustering
)
process.patJetsAK8PFCHSSoftDropSubjets.addTagInfos = True
