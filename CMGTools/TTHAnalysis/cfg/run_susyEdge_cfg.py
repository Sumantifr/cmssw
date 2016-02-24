##################i#######################################
##       CONFIGURATION FOR JZBEDGE TREES                ##
## skim condition: >= 2 loose leptons, no pt cuts or id ##
##########################################################
# 1.- Common Stuff
# 2.- Lepton Skimming
# 3.- Lepton Analyzer
# 4.- Isolation
# 5.- Photon Analyzer
# 6.- TTHLepEventAnalyzer
# 7.- TTHJetTauAnalyzer
# 8.- Declustering Analyzer
# 9.- TreeProducer SUSY Multilepton
# 10.- Insertion in SUSY Core Sequence
# 11.- Trigger definition an bits
# 12.- Sample imports
# 13.- Sequence
# 14.- Services
# 15.- HOW TO RUN
# 16.- Services
##########################################################
##################Commont stuff###########################
##########################################################
import os
import PhysicsTools.HeppyCore.framework.config as cfg
from CMGTools.TTHAnalysis.analyzers.susyCore_modules_cff import *


##########################################################
##################Lepton Skimming#########################
##########################################################
ttHLepSkim.minLeptons = 2
ttHLepSkim.maxLeptons = 999
#ttHLepSkim.idCut  = ""
#ttHLepSkim.ptCuts = []

runSMS = False

##########################################################
##################Lepton Analyzier########################
##########################################################
lepAna.doMiniIsolation = True
lepAna.packedCandidates = 'packedPFCandidates'
lepAna.miniIsolationPUCorr = 'rhoArea'
lepAna.miniIsolationVetoLeptons = None # use 'inclusive' to veto inclusive leptons and their footprint in all isolation cones
## this is new lepAna.loose_electron_id  = "POG_MVA_ID_Phys14_NonTrig_VLoose"
lepAna.loose_electron_id = "POG_MVA_ID_Spring15_NonTrig_VLooseIdEmu"
lepAna.loose_muon_id      = "POG_ID_Loose"
lepAna.loose_electron_eta = 2.5
lepAna.loose_electron_pt  = 15
lepAna.loose_muon_eta     = 2.4
lepAna.loose_muon_pt      = 15
#These are the nominal cuts
#lepAna.loose_electron_dxy = 0.05
#lepAna.loose_electron_dz  =  0.1
#lepAna.loose_muon_dxy     = 0.05
#lepAna.loose_muon_dz      = 0.1
#But I go more inclusive for further selection in the edgeFriends
lepAna.loose_electron_dxy = 0.15
lepAna.loose_electron_dz  =  0.5
lepAna.loose_muon_dxy     = 0.15
lepAna.loose_muon_dz      =  0.5


jetAna.recalibrateJets = True
#jetAna.calculateSeparateCorrections = True
jetAna.calculateType1METCorrection = True
metAna.recalibrate = "type1"

jetAna.mcGT     = "Summer15_25nsV6_MC"
jetAna.dataGT   = "Summer15_25nsV6_DATA"

##########################################################
######################Isolation###########################
##########################################################
isolation = "miniIso"
#isolation = "ptRel"
if isolation == "ptRel": 
    # delay isolation cut for leptons of pt > 10, for which we do pTrel recovery
    lepAna.loose_muon_isoCut     = lambda muon : muon.relIso03 < 0.5 or muon.pt() > 10
    lepAna.loose_electron_isoCut = lambda elec : elec.relIso03 < 0.5 or elec.pt() > 10
    # in the cleaning, keep the jet if the lepton fails relIso or ptRel
    jetAna.jetLepArbitration = lambda jet,lepton : (
        lepton if (lepton.relIso03 < 0.4 or ptRelv1(lepton.p4(),jet.p4()) > 5) else jet
    )
    ttHCoreEventAna.leptonMVAKindTTH = "SusyWithBoost"
    ttHCoreEventAna.leptonMVAKindSusy = "SusyWithBoost" 
    ttHCoreEventAna.leptonMVAPathTTH = "CMGTools/TTHAnalysis/macros/leptons/trainingPHYS14leptonMVA_PHYS14eleMVA_MiniIso_ttH/weights/%s_BDTG.weights.xml"
    ttHCoreEventAna.leptonMVAPathSusy = "CMGTools/TTHAnalysis/macros/leptons/trainingPHYS14leptonMVA_PHYS14eleMVA_MiniIso_SusyT1/weights/%s_BDTG.weights.xml"
    # insert a second skimmer after the jet cleaning 
    ttHLepSkim2 = cfg.Analyzer(
        ttHLepSkimmer, name='ttHLepSkimmer2',
        minLeptons = 2,
        maxLeptons = 999,
        )
    susyCoreSequence.insert(susyCoreSequence.index(jetAna)+1, ttHLepSkim2)
elif isolation == "miniIso": 
    lepAna.loose_muon_isoCut     = lambda muon : muon.miniRelIso < 0.4
    lepAna.loose_electron_isoCut = lambda elec : elec.miniRelIso < 0.4
elif isolation == None:
    lepAna.loose_muon_isoCut     = lambda muon : True
    lepAna.loose_electron_isoCut = lambda elec : True
else:
    # nothing to do, will use normal relIso03
    pass

##########################################################
##################Photon Analyzer#########################
##########################################################
# Switch off slow photon MC matching
photonAna.do_mc_match = False



##########################################################
###############TTHLepEventAnalyzer########################
##########################################################
## Event Analyzer for susy multi-lepton (at the moment, it's the TTH one)
from CMGTools.TTHAnalysis.analyzers.ttHLepEventAnalyzer import ttHLepEventAnalyzer
ttHEventAna = cfg.Analyzer(
    ttHLepEventAnalyzer, name="ttHLepEventAnalyzer",
    minJets25 = 0,
    )

##########################################################
#################TTHJetTauAnalyzer########################
##########################################################
## JetTau analyzer, to be called (for the moment) once bjetsMedium are produced
from CMGTools.TTHAnalysis.analyzers.ttHJetTauAnalyzer import ttHJetTauAnalyzer
ttHJetTauAna = cfg.Analyzer(
    ttHJetTauAnalyzer, name="ttHJetTauAnalyzer",
    )


##########################################################
############Insertion in susyCoreSequence#################
##########################################################
## Insert the FatJet, SV, HeavyFlavour analyzers in the sequence
susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna), 
                        ttHFatJetAna)
susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna), 
                        ttHSVAna)
susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna), 
                        ttHHeavyFlavourHadronAna)

##########################################################
#################DeClustering Analyzer####################
##########################################################
## Insert declustering analyzer
from CMGTools.TTHAnalysis.analyzers.ttHDeclusterJetsAnalyzer import ttHDeclusterJetsAnalyzer
ttHDecluster = cfg.Analyzer(
    ttHDeclusterJetsAnalyzer, name='ttHDecluster',
    lepCut     = lambda lep,ptrel : lep.pt() > 10,
    maxSubjets = 6, # for exclusive reclustering
    ptMinSubjets = 5, # for inclusive reclustering
    drMin      = 0.2, # minimal deltaR(l,subjet) required for a successful subjet match
    ptRatioMax = 1.5, # maximum pt(l)/pt(subjet) required for a successful match
    ptRatioDiff = 0.1,  # cut on abs( pt(l)/pt(subjet) - 1 ) sufficient to call a match successful
    drMatch     = 0.02, # deltaR(l,subjet) sufficient to call a match successful
    ptRelMin    = 5,    # maximum ptRelV1(l,subjet) sufficient to call a match successful
    prune       = True, # also do pruning of the jets 
    pruneZCut       = 0.1, # pruning parameters (usual value in CMS: 0.1)
    pruneRCutFactor = 0.5, # pruning parameters (usual value in CMS: 0.5)
    verbose     = 0,   # print out the first N leptons
    jetCut = lambda jet : jet.pt() > 20,
    mcPartonPtCut = 20,
    mcLeptonPtCut =  5,
    mcTauPtCut    = 15,
    )

susyCoreSequence.insert(susyCoreSequence.index(ttHFatJetAna)+1, ttHDecluster)

##########################################################
#############TreeProducer SusyMultilepton#################
##########################################################
from CMGTools.TTHAnalysis.analyzers.treeProducerSusyMultilepton import * 
## Tree Producer
treeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='treeProducerSusyEdge',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     PDFWeights = PDFWeights,
     globalVariables = susyMultilepton_globalVariables,
     globalObjects = susyMultilepton_globalObjects,
     collections = susyMultilepton_collections,
)




## histo counter
if not runSMS:
    susyCoreSequence.insert(susyCoreSequence.index(skimAnalyzer),
                            susyCounter)
else:
    susyCoreSequence.insert(susyCoreSequence.index(susyScanAna)+1,susyCounter)

# HBHE new filter
from CMGTools.TTHAnalysis.analyzers.hbheAnalyzer import hbheAnalyzer
hbheAna = cfg.Analyzer(
    hbheAnalyzer, name="hbheAnalyzer", IgnoreTS4TS5ifJetInLowBVRegion=False
    )
if not runSMS:
    susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna),hbheAna)
    treeProducer.globalVariables.append(NTupleVariable("hbheFilterNew50ns", lambda ev: ev.hbheFilterNew50ns, int, help="new HBHE filter for 50 ns"))
    treeProducer.globalVariables.append(NTupleVariable("hbheFilterNew25ns", lambda ev: ev.hbheFilterNew25ns, int, help="new HBHE filter for 25 ns"))
    treeProducer.globalVariables.append(NTupleVariable("hbheFilterIso", lambda ev: ev.hbheFilterIso, int, help="HBHE iso-based noise filter"))


## not needed for Oct05 and prompt-v4doT1METCorr = True
## not needed for Oct05 and prompt-v4if doT1METCorr:
## not needed for Oct05 and prompt-v4    jetAna.calculateType1METCorrection = True
## not needed for Oct05 and prompt-v4    metAna.recalibrate = "type1"
## not needed for Oct05 and prompt-v4    #metAna.old74XMiniAODs = old74XMiniAODs


##########################################################
############Trigger Definition and bits  #################
##########################################################
from CMGTools.RootTools.samples.triggers_13TeV_Spring15 import *
#from CMGTools.RootTools.samples.triggers_8TeV import triggers_1mu_8TeV, triggers_mumu_8TeV, triggers_mue_8TeV, triggers_ee_8TeV;
triggerFlagsAna.triggerBits = {
    'mu17mu8' : triggers_mu17mu8,
    'mu17mu8_dz' : triggers_mu17mu8_dz,
    'mu17tkmu8_dz' : triggers_mu17tkmu8_dz,
    'mu27tkmu8' : triggers_mumu_noniso_50ns,
    'mu17el12' : triggers_mu17el12,
    'mu30ele30' : triggers_mu30ele30,
    'el17el12_dz' : triggers_el17el12_dz,
    'el23el12_dz' : triggers_el23el12_dz,
    'ele33ele33' : triggers_doubleele33,
    'mu8el17' : triggers_mu8el17,
    'mu8el23' : triggers_mu8el23,
    'pfht200' : triggers_pfht200,
    'pfht250' : triggers_pfht250,
    'pfht300' : triggers_pfht300,
    'pfht350' : triggers_HT350,
    'pfht400' : triggers_pfht400,
    'pfht475' : triggers_HT475,
    'pfht600' : triggers_HT600,
    'pfht650' : triggers_pfht650,
    'pfht800' : triggers_HT800,
    'pfht900' : triggers_HT900,
    'at57' : triggers_at57,
    'at55' : triggers_at55,
    'at53' : triggers_at53,
    'at52' : triggers_at52,
    'at51' : triggers_at51,
    'DoubleMu' : triggers_mumu,
    #'DoubleMuNoIso' : triggers_mumu_noniso,
    'DoubleEl' : triggers_ee,
    'MuEG'     : triggers_mue,
    'DoubleMuHT'     : triggers_mumu_ht,
    'DoubleElHT'     : triggers_ee_ht,
    'MuEGHT'     : triggers_mue_ht,
    'SingleMu' : triggers_1mu_iso,
    'SingleEl' : triggers_FR_1e_iso,
    'HTMET' : triggers_HT350_MET120,
    'HTJet' : triggers_htjet,
}

if runSMS:
    susyCoreSequence.remove(triggerFlagsAna)
    susyCoreSequence.remove(triggerAna)
    susyCoreSequence.remove(eventFlagsAna)


##########################################################
################### Sample imports  ######################
##########################################################
from CMGTools.RootTools.samples.samples_13TeV_RunIISpring15MiniAODv2 import *
from CMGTools.RootTools.samples.samples_13TeV_DATA2015 import *
from CMGTools.RootTools.samples.samples_13TeV_signals import *


##########################################################
######################## Sequence ########################
##########################################################
sequence = cfg.Sequence(susyCoreSequence+[
        ttHJetTauAna,
        ttHEventAna,
        treeProducer,
    ])

preprocessor = None
##########################################################
###################### HOW TO RUN#########################
##########################################################
selectedComponents = [] 

selectedComponents = [DYJetsToLL_M50]


for comp in selectedComponents:
    comp.splitFactor = 100
    comp.finesplitFactor = 4

runData = False
if runData:
    ## json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
    ## json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON_v3.txt"
    ## json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258714_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
    ## json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt"

    ## this is the JSON for the full dataset of the year 2015. 2.11 inverse femtobarns
    #json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
    ## this is the JSON for moriond
    json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt"

    ## # Run2015C, 25 ns, 3.8T
    ## processing = "Run2015C-05Oct2015-v1"; short = "Run2015C_miniAODv2"; run_ranges = [ (246908,260627) ]; useAAA=False; is50ns=False; triggerFlagsAna.checkL1Prescale = False;

    # Run2015D, 25 ns, 3.8 T re-miniAOD
    # processing = "Run2015C_25ns-05Oct2015-v1"; short = "Run2015C_25ns_05Oct_v1"; run_ranges = [ (246908,260627) ]; useAAA=False; is50ns=False
    # processing2 = "Run2015C_25ns-05Oct2015-v1"; short2 = "Run2015C_25ns_05Oct_v1"

    # Run2015D, 25 ns, 3.8 T re-miniAOD
    ##processing = "Run2015D-05Oct2015-v1"; short = "Run2015D_05Oct_v1"; run_ranges = [ (246908,260627) ]; useAAA=False; is50ns=False
    ##processing2 = "Run2015D-05Oct2015-v2"; short2 = "Run2015D_05Oct_v2"

    # Run2015D, prompt-v4
    ## processing = "Run2015D-PromptReco-v4"; short = "Run2015D_v4"; run_ranges = [ (246908,260627) ]; useAAA=False; is50ns=False
    ## processing2 = "Run2015D-PromptReco-v4"; short2 = "Run2015D_v4"; 

    compSelection = ""; compVeto = ""
    DatasetsAndTriggers = []
    selectedComponents = []; vetos = []  

    useAAA = True
    is50ns = False
    processings = ["Run2015C_25ns-05Oct2015-v1", "Run2015D-05Oct2015-v1", "Run2015D-PromptReco-v4"]
    shorts      = ["Run2015C_25ns-05Oct_v1"    , "Run2015D-05Oct_v1"    , "Run2015D_v4"]
    run_ranges  = [ (246908,260628) ]
 
    ## hadtriggers = triggers_pfht200 + triggers_pfht250 + triggers_pfht300 + triggers_ht350 + triggers_pfht400 + triggers_ht475 + triggers_ht600 + triggers_HT800 + triggers_HT900 + triggers_at57 + triggers_at55 + triggers_at53 + triggers_at52 + triggers_at51
    ## DatasetsAndTriggers.append( ("DoubleMuon", triggers_mu17mu8 + triggers_mu17mu8_dz + triggers_mu17tkmu8_dz + triggers_mumu_noniso_50ns) )
    ## DatasetsAndTriggers.append( ("DoubleEG",   triggers_el17el12_dz + triggers_el23el12_dz + triggers_doubleele33) )
    ## DatasetsAndTriggers.append( ("MuonEG",     triggers_mu8el17 + triggers_mu8el23 + triggers_mu17el12 + triggers_mu30ele30) )
    ## DatasetsAndTriggers.append( ("JetHT", hadtriggers) )
    ## DatasetsAndTriggers.append( ("HTMHT", hadtriggers) )
    DatasetsAndTriggers.append( ("DoubleMuon",[] ) )
    DatasetsAndTriggers.append( ("DoubleEG",  [] ) )
    DatasetsAndTriggers.append( ("MuonEG",   []  ) )
    DatasetsAndTriggers.append( ("JetHT", []) )
    DatasetsAndTriggers.append( ("HTMHT", []) )

    for pd,triggers in DatasetsAndTriggers:
        for pro in processings:
            i = processings.index(pro)
            short_use = shorts[i]
            if pro == "Run2015D-05Oct2015-v1" and pd == 'MuonEG':
                pro = pro.replace('v1','v2')
                short_use = short_use.replace('v1', 'v2')
            ## short_use = short2 if pd == "MuonEG" else short
            ## processing_use = processing2 if pd == "MuonEG" else processing
            #if pd == "MuonEG": processing = processing.replace("v1", "v2")
            for run_range in run_ranges:
                label = "runs_%d_%d" % run_range if run_range[0] != run_range[1] else "run_%d" % (run_range[0],)
                compname = pd+"_"+short_use+"_"+label
                ## if ((compSelection and not re.search(compSelection, compname)) or
                ##     (compVeto      and     re.search(compVeto,      compname))):
                ##         print "Will skip %s" % (compname)
                ##         continue
                print 'looking for dataset', "/"+pd+"/"+(pro) +"/MINIAOD"
                print 'using short', compname
                comp = kreator.makeDataComponent(compname, 
                                                 "/"+pd+"/"+(pro) +"/MINIAOD", 
                                                 "CMS", ".*root", 
                                                 json=json, 
                                                 run_range=run_range, 
                                                 #triggers=triggers[:], vetoTriggers = vetos[:],
                                                 triggers=[], vetoTriggers = vetos[:],
                                                 useAAA=useAAA)
                print "Will process %s (%d files)" % (comp.name, len(comp.files))
#                print "\ttrigger sel %s, veto %s" % (triggers, vetos)
                comp.splitFactor = len(comp.files)/5
                comp.fineSplitFactor = 1
                selectedComponents.append( comp )
            #vetos += triggers
    if json is None:
        susyCoreSequence.remove(jsonAna)

if runSMS:
    jetAna.mcGT = "MCRUN2_74_V9_FASTSIM_291115"
    jetAnaScaleUp.mcGT = "MCRUN2_74_V9_FASTSIM_291115"
    jetAnaScaleDown.mcGT = "MCRUN2_74_V9_FASTSIM_291115"
    jetAna.applyL2L3Residual = False
    jetAnaScaleUp.applyL2L3Residual = False
    jetAnaScaleDown.applyL2L3Residual = False
    useAAA = False ## accesses only files on EOS or at CERN

    print 'I\'m in the sms test thing here!!'
    selectedComponents = [SMS_T6bbllslepton_mSbottom_400To550_mLSP_200To500_miniAODv2]#, SMS_T6bbllslepton_mSbottom_600To900_mLSP_200To800_miniAODv2]
    comp.files = comp.files[:1]
    #for comp in selectedComponents:
    #    comp.splitFactor = 500
    #    #comp.finesplitFactor = 4
    ##comp.splitFactor = 500
    ##comp.fineSplitFactor = 4

from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption
test = getHeppyOption('test')
#test = '74X-MC'

test= '74X-MC'
if test == 'synch':
    print 'I\'m in the synch test thing here!!'
    #comp = TTLep_pow
    #selectedComponents = [comp]
    comp = selectedComponents[0]
    selectedComponents = [comp]
    #comp.files = comp.files[:1]
    comp.files = comp.files[:1]#[ '/afs/cern.ch/work/m/mdunser/public/synchFiles/004613BA-C46D-E511-9EB6-001E67248732.root' ]
    #comp.finesplitFactor = 10
    #comp.finesplitFactor = 4
elif test == '74X-MC':
    #what = getHeppyOption("sample")
    what = 'TT'
    if what == "TTLep":
        selectedComponents = [ TTLep_pow ]
        comp = selectedComponents[0]
        comp.files = [ '/store/mc/RunIISpring15DR74/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0C1B984D-F408-E511-872E-0002C90B7F2E.root' ]
        tmpfil = os.path.expandvars("/tmp/$USER/0C1B984D-F408-E511-872E-0002C90B7F2E.root")
        if not os.path.exists(tmpfil):
            os.system("xrdcp root://eoscms//eos/cms%s %s" % (comp.files[0],tmpfil))
        comp.files = [ tmpfil ]
    elif what == "TT":
        ttHLepSkim.minLeptons = 0
        selectedComponents = [ TT_bx25 ]
    elif what == "Z":
        selectedComponents = [ ZEE_bx25, ZMM_bx25, ZTT_bx25 ]
    else:
        selectedComponents = RelVals740
    if not getHeppyOption("all"):
        for comp in selectedComponents:
            comp.files = comp.files[:1]
            comp.splitFactor = 1
            comp.fineSplitFactor = 1 if getHeppyOption("single") else 4
    selectedComponents = [ DYJetsToLL_M50]
    for comp in selectedComponents:
        comp.splitFactor = 100
        comp.fineSplitFactor = 4 
elif test == '74X-Data':
    #selectedComponents = [DoubleMuon_Run2015B, MuonEG_Run2015B, DoubleEG_Run2015B, JetHT_Run2015B, HTMHT_Run2015B]
    selectedComponents = [JetHT_Run2015B]#, HTMHT_Run2015B]
    #selectedComponents = [DoubleMuon_Run2015B]
    eventFlagsAna.processName = 'HLT'
    jetAna.recalibrateJets = False
    jetAna.smearJets       = False 
    photonAna.do_mc_match = False
    vertexAna.keepFailingEvents = True 
    jsonAna.useLumiBlocks = True
    #sequence.remove(jsonAna)
    for comp in selectedComponents:
        comp.isMC = False
        comp.isData = True
        #comp.files = ["/afs/cern.ch/work/p/pablom/public/E42FEF61-6E27-E511-B93A-02163E0143C0.root"]
        comp.files = comp.files[:1]

elif test == "express":
    selectedComponents = [ MuEG_740p9 ]
    comp = selectedComponents[0]
    comp.files = [ 'root://eoscms//eos/cms/store/express/Run2015A/ExpressPhysics/FEVT/Express-v1/000/246/908/00000/04B152E7-DE09-E511-8B18-02163E011D4A.root' ]
    comp.name  = 'ExpressPhysics'
    comp.triggers = []
    comp.json     = None
    jetAna.recalibrateJets = False 
    ttHLepSkim.minLeptons = 0
    # preprocessor cfg to be created with
    #    cmsDriver.py miniAOD-data -s PAT --data --runUnscheduled --eventcontent MINIAOD --conditions GR_P_V56 --no_exec
    #    sed -i 's/process.MINIAODoutput_step/process.endpath/' miniAOD-data_PAT.py
    from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
    preprocessor = CmsswPreprocessor("miniAOD-data_PAT.py")

>>>>>>> cms-edgeZ/towards-CMSSW_7_4_16


##########################################################
############# Output Module && services###################
##########################################################
outputService=[]
from PhysicsTools.HeppyCore.framework.services.tfile import TFileService
output_service = cfg.Service(
    TFileService,
    'outputfile',
    name="outputfile",
    fname='treeProducerSusyEdge/tree.root',
    option='recreate'
    )    
outputService.append(output_service)

##########################################################
###################### Services ##########################
##########################################################
# the following is declared in case this cfg is used in input to the heppy.py script
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
from CMGTools.TTHAnalysis.tools.EOSEventsWithDownload import EOSEventsWithDownload
event_class = EOSEventsWithDownload
EOSEventsWithDownload.aggressive = 2 # always fetch if running on Wigner

if getHeppyOption("nofetch"):
    event_class = Events 
config = cfg.Config( components = selectedComponents,
                     sequence = sequence,
                     services = outputService, 
                     preprocessor = preprocessor, 
                     events_class = event_class)



