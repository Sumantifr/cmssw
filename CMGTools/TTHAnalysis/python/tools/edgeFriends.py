## if MT2 doesnt work, put this line in the MT2 file: from ROOT.heppy import Davismt2

from CMGTools.TTHAnalysis.treeReAnalyzer import *
from CMGTools.TTHAnalysis.tools.eventVars_MT2 import *
print 'loading stuff for MT2'
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()
print 'done loading MT2 stuff.'

import copy

class edgeFriends:
    def __init__(self,label,tightLeptonSel,cleanJet,isMC=False):
        self.label = "" if (label in ["",None]) else ("_"+label)
        self.tightLeptonSel = tightLeptonSel
        self.cleanJet = cleanJet
        self.isMC = isMC
        ## with nvtx self.puFile = open("/afs/cern.ch/work/m/mdunser/public/puWeighting/puWeightsVinceLumi1p28.txt","r")
        self.puFile = open("/afs/cern.ch/work/m/mdunser/public/puWeighting/puWeightsOfficialPrescription.txt","r")
        self.pu_dict = eval(self.puFile.read())
        self.puFile.close()
        ## for 1.3 fb-1 self.beamHaloListFile = open("/afs/cern.ch/work/m/mdunser/public/beamHalo/beamHaloEvents_DoubleLep_JetHT_HTMHT.txt","r")
        ## nov14 list self.beamHaloListFile = open("/afs/cern.ch/work/m/mdunser/public/beamHalo/fullDataset/allFullData.txt","r")
        self.beamHaloListFile = open("/afs/cern.ch/user/p/pablom/public/Filters_27_01_2016/csc2015.txt", "r")
        self.fourthBadEESuperCrystalFile = open("/afs/cern.ch/user/p/pablom/public/Filters_27_01_2016/ecalscn1043093.txt","r")
        self.badResolutionTrackTaggerFile = open("/afs/cern.ch/user/p/pablom/public/Filters_27_01_2016/badResolutionTrack.txt","r")
        self.badMuonTrackTaggerFile = open("/afs/cern.ch/user/p/pablom/public/Filters_27_01_2016/muonBadTrack.txt","r")

	if not self.isMC:
            self.beamHaloSet = set()
            self.fourthBadEESuperCrystalSet = set()
            self.badResolutionTrackTaggerSet = set()
            self.badMuonTrackTaggerSet = set()
            for i in list(self.beamHaloListFile):
                self.beamHaloSet.add(i.rstrip('\n'))
            for i in list(self.fourthBadEESuperCrystalFile):
                self.fourthBadEESuperCrystalSet.add(i.rstrip('\n'))
            for i in list(self.badResolutionTrackTaggerFile):
                self.badResolutionTrackTaggerSet.add(i.rstrip('\n'))
            for i in list(self.badMuonTrackTaggerFile):
                self.badMuonTrackTaggerSet.add(i.rstrip('\n'))
        self.beamHaloListFile.close()
        self.fourthBadEESuperCrystalFile.close()
        self.badResolutionTrackTaggerFile.close()
        self.badMuonTrackTaggerFile.close()
        ## load the pdf file for the likelihood
        self.lh_file = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version11.root')
        # these are the kernel pdfs
        self.h_lh_zpt_data = copy.deepcopy(self.lh_file.Get('em_data_pdf_histo_zpt_ds_cuts_of_sr_met150')); self.h_lh_zpt_data.Scale(1./self.h_lh_zpt_data.Integral())
        self.h_lh_met_data = copy.deepcopy(self.lh_file.Get('em_data_pdf_histo_met_ds_cuts_of_sr_met150')); self.h_lh_met_data.Scale(1./self.h_lh_met_data.Integral())
        self.h_lh_mlb_data = copy.deepcopy(self.lh_file.Get('em_data_pdf_histo_mlb_ds_cuts_of_sr_met150')); self.h_lh_mlb_data.Scale(1./self.h_lh_mlb_data.Integral())
        self.h_lh_ldr_data = copy.deepcopy(self.lh_file.Get('em_data_pdf_histo_ldr_ds_cuts_of_sr_met150')); self.h_lh_ldr_data.Scale(1./self.h_lh_ldr_data.Integral())
        self.h_lh_zpt_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_pdf_histo_zpt_ds_cuts_of_sr_met150'  )); self.h_lh_zpt_mc  .Scale(1./self.h_lh_zpt_mc  .Integral())
        self.h_lh_met_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_pdf_histo_met_ds_cuts_of_sr_met150'  )); self.h_lh_met_mc  .Scale(1./self.h_lh_met_mc  .Integral())
        self.h_lh_mlb_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_pdf_histo_mlb_ds_cuts_of_sr_met150'  )); self.h_lh_mlb_mc  .Scale(1./self.h_lh_mlb_mc  .Integral())
        self.h_lh_ldr_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_pdf_histo_ldr_ds_cuts_of_sr_met150'  )); self.h_lh_ldr_mc  .Scale(1./self.h_lh_ldr_mc  .Integral())
        # these are the analytical pdfs
        self.h_lh_ana_zpt_data = copy.deepcopy(self.lh_file.Get('em_data_fit_histo_zpt_ds_cuts_of_sr_met150__lepsZPt_Edge')); self.h_lh_ana_zpt_data.Scale(1./self.h_lh_ana_zpt_data.Integral())
        self.h_lh_ana_met_data = copy.deepcopy(self.lh_file.Get('em_data_fit_histo_met_ds_cuts_of_sr_met150__met_pt'      )); self.h_lh_ana_met_data.Scale(1./self.h_lh_ana_met_data.Integral())
        self.h_lh_ana_mlb_data = copy.deepcopy(self.lh_file.Get('em_data_fit_histo_mlb_ds_cuts_of_sr_met150__sum_mlb_Edge')); self.h_lh_ana_mlb_data.Scale(1./self.h_lh_ana_mlb_data.Integral())
        self.h_lh_ana_ldr_data = copy.deepcopy(self.lh_file.Get('em_data_fit_histo_ldr_ds_cuts_of_sr_met150__lepsDR_Edge' )); self.h_lh_ana_ldr_data.Scale(1./self.h_lh_ana_ldr_data.Integral())
        self.h_lh_ana_zpt_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_fit_histo_zpt_ds_cuts_of_sr_met150__lepsZPt_Edge'  )); self.h_lh_ana_zpt_mc  .Scale(1./self.h_lh_ana_zpt_mc  .Integral())
        self.h_lh_ana_met_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_fit_histo_met_ds_cuts_of_sr_met150__met_pt'        )); self.h_lh_ana_met_mc  .Scale(1./self.h_lh_ana_met_mc  .Integral())
        self.h_lh_ana_mlb_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_fit_histo_mlb_ds_cuts_of_sr_met150__sum_mlb_Edge'  )); self.h_lh_ana_mlb_mc  .Scale(1./self.h_lh_ana_mlb_mc  .Integral())
        self.h_lh_ana_ldr_mc   = copy.deepcopy(self.lh_file.Get('tt_mc_fit_histo_ldr_ds_cuts_of_sr_met150__lepsDR_Edge'   )); self.h_lh_ana_ldr_mc  .Scale(1./self.h_lh_ana_ldr_mc  .Integral())
        self.lh_file.Close()

    def listBranches(self):
        label = self.label
        biglist = [ ("nLepTight"+label, "I"), ("nJetSel"+label, "I"), ("nPairLep"+label, "I"),
                 ("iLT"+label,"I",20,"nLepTight"+label), 
                 ("iJ"+label,"I",20,"nJetSel"+label), # index >= 0 if in Jet; -1-index (<0) if in DiscJet
                 ("nLepGood20"+label, "I"), ("nLepGood20T"+label, "I"),
                 ("nJet35"+label, "I"), ("htJet35j"+label), ("nBJetLoose35"+label, "I"), ("nBJetMedium35"+label, "I"), 
                 ("iL1T"+label, "I"), ("iL2T"+label, "I"), 
                 ("lepsMll"+label, "F"), ("lepsJZB"+label, "F"), ("lepsJZB_raw"+label, "F"), ("lepsJZB_recoil"+label, "F"), ("lepsDR"+label, "F"), ("lepsMETRec"+label, "F"), ("lepsZPt"+label, "F"), ("metl1DPhi"+label, "F"), ("metl2DPhi"+label, "F"), ("lepsDPhi"+label, "F"),
                 ("Lep1_pt"+label, "F"), ("Lep1_eta"+label, "F"), ("Lep1_phi"+label, "F"), ("Lep1_miniRelIso"+label, "F"), ("Lep1_pdgId"+label, "I"), ("Lep1_mvaIdSpring15"+label, "F"), ("Lep1_minTauDR"+label, "F"),
                 ("Lep2_pt"+label, "F"), ("Lep2_eta"+label, "F"), ("Lep2_phi"+label, "F"), ("Lep2_miniRelIso"+label, "F"), ("Lep2_pdgId"+label, "I"), ("Lep2_mvaIdSpring15"+label, "F"), ("Lep2_minTauDR"+label, "F"),
                 ("PileupW"+label, "F"), ("min_mlb1"+label, "F"), ("min_mlb2"+label, "F"), ("sum_mlb"+label, "F"), ("st"+label,"F"), ("srID"+label, "I"), ("mt2"+label, "F"),
                 ("lh_zpt_data"+label, "F") , ("lh_met_data"+label, "F") , ("lh_mlb_data"+label, "F") , ("lh_ldr_data"+label, "F") ,
                 ("lh_zpt_mc"+label  , "F") , ("lh_met_mc"+label  , "F") , ("lh_mlb_mc"+label  , "F") , ("lh_ldr_mc"+label  , "F") ,
                 ("lh_ana_zpt_data"+label, "F") , ("lh_ana_met_data"+label, "F") , ("lh_ana_mlb_data"+label, "F") , ("lh_ana_ldr_data"+label, "F") ,
                 ("lh_ana_zpt_mc"+label  , "F") , ("lh_ana_met_mc"+label  , "F") , ("lh_ana_mlb_mc"+label  , "F") , ("lh_ana_ldr_mc"+label  , "F") ,

                 ("cum_zpt_data"+label, "F") , ("cum_met_data"+label, "F") , ("cum_mlb_data"+label, "F") , ("cum_ldr_data"+label, "F") ,
                 ("cum_zpt_mc"+label  , "F") , ("cum_met_mc"+label  , "F") , ("cum_mlb_mc"+label  , "F") , ("cum_ldr_mc"+label  , "F") ,
                 ("cum_ana_zpt_data"+label, "F") , ("cum_ana_met_data"+label, "F") , ("cum_ana_mlb_data"+label, "F") , ("cum_ana_ldr_data"+label, "F") ,
                 ("cum_ana_zpt_mc"+label  , "F") , ("cum_ana_met_mc"+label  , "F") , ("cum_ana_mlb_mc"+label  , "F") , ("cum_ana_ldr_mc"+label  , "F")
                 
                 ]
        ## for lfloat in 'pt eta phi miniRelIso pdgId'.split():
        ##     if lfloat == 'pdgId':
        ##         biglist.append( ("Lep"+label+"_"+lfloat,"I", 10, "nPairLep"+label) )
        ##     else:
        ##         biglist.append( ("Lep"+label+"_"+lfloat,"F", 10, "nPairLep"+label) )
        for jfloat in "pt eta phi mass btagCSV rawPt".split():
            biglist.append( ("JetSel"+label+"_"+jfloat,"F",20,"nJetSel"+label) )
        if self.isMC:
            biglist.append( ("JetSel"+label+"_mcPt",     "F",20,"nJetSel"+label) )
            biglist.append( ("JetSel"+label+"_mcFlavour","I",20,"nJetSel"+label) )
            biglist.append( ("JetSel"+label+"_mcMatchId","I",20,"nJetSel"+label) )
        return biglist
    def __call__(self,event):
        leps  = [l for l in Collection(event,"LepGood","nLepGood")]
        lepso = [l for l in Collection(event,"LepOther","nLepOther")]
        jetsc = [j for j in Collection(event,"Jet","nJet")]
        jetsd = [j for j in Collection(event,"DiscJet","nDiscJet")]
        metco = [m for m in Collection(event,"metcJet","nDiscJet")]
        (met, metphi)  = event.met_pt, event.met_phi
        (met_raw, metphi_raw)  = event.met_rawPt, event.met_rawPhi
        if self.isMC:
            gentaus  = [t for t in Collection(event,"genTau","ngenTau")]
            ntrue = event.nTrueInt
        ## nvtx = event.nVert
        metp4 = ROOT.TLorentzVector()
        metp4.SetPtEtaPhiM(met,0,metphi,0)
        metp4_raw = ROOT.TLorentzVector()
        metp4_raw.SetPtEtaPhiM(met_raw,0,metphi_raw,0)
        ret = {}; jetret = {}; 
        lepret = {}
        
        #
        ### Define tight leptons
        ret["iLT"] = []; ret["nLepGood20T"] = 0

        # ====================
        # do pileupReweighting
        # ====================
        puWt = self.pu_dict[ntrue] if self.isMC else 1.
        #if puWt > 10: puWt = 10.
        ret["PileupW"] = puWt

        # ===============================
        # new, simpler sorting of leptons
        # ===============================
        for il,lep in enumerate(leps):
            clean = True
            if self.tightLeptonSel(lep) and clean:
                ret["iLT"].append(il)
                if lep.pt > 20: ret["nLepGood20T"] += 1
        # other leptons, negative indices
        for il,lep in enumerate(lepso):
            clean = True
            if self.tightLeptonSel(lep) and clean:
                ret["iLT"].append(-1-il)
                if lep.pt > 20: ret["nLepGood20T"] += 1
        ret["nLepTight"] = len(ret["iLT"])
        #
        # sort the leptons by pT:
        ret["iLT"].sort(key = lambda idx : leps[idx].pt if idx >= 0 else lepso[-1-idx].pt, reverse = True)

        ## search for the lepton pair
        #lepst  = [ leps [il] for il in ret["iLT"] ]

        lepst = []
        for il in ret['iLT']:
            if il >=0: 
                lepst.append(leps[il])
            else: 
                lepst.append(lepso[-1-il])
        #
      
        iL1iL2 = self.getPairVariables(lepst, metp4, metp4_raw)
        ret['iL1T'] = ret["iLT"][ iL1iL2[0] ] if (len(ret["iLT"]) >=1 and iL1iL2[0] != -999) else -999
        ret['iL2T'] = ret["iLT"][ iL1iL2[1] ] if (len(ret["iLT"]) >=2 and iL1iL2[1] != -999) else -999
        ret['lepsMll'] = iL1iL2[2] 
        ret['lepsJZB'] = iL1iL2[3] 
        ret['lepsJZB_raw'] = iL1iL2[4] 
        ret['lepsDR'] = iL1iL2[5] 
        ret['lepsMETRec'] = iL1iL2[6] 
        ret['lepsZPt'] = iL1iL2[7] 
        ret['lepsDPhi'] = iL1iL2[8]


        #print 'new event =================================================='
        l1 = ROOT.TLorentzVector()
        l2 = ROOT.TLorentzVector()
        ltlvs = [l1, l2]
        lepvectors = []

        for lfloat in 'pt eta phi miniRelIso pdgId mvaIdSpring15'.split():
            if lfloat == 'pdgId':
                lepret["Lep1_"+lfloat+self.label] = -99
                lepret["Lep2_"+lfloat+self.label] = -99
            else:
                lepret["Lep1_"+lfloat+self.label] = -42.
                lepret["Lep2_"+lfloat+self.label] = -42.
        if ret['iL1T'] != -999 and ret['iL2T'] != -999:
            ret['nPairLep'] = 2
            # compute the variables for the two leptons in the pair
            lcount = 1
            for idx in [ret['iL1T'], ret['iL2T']]:
                lep = leps[idx] if idx >= 0 else lepso[-1-idx]
                minDRTau = 99.
                if self.isMC:
                    for tau in gentaus:
                        tmp_dr = deltaR(lep, tau)
                        if tmp_dr < minDRTau:
                            minDRTau = tmp_dr
                for lfloat in 'pt eta phi miniRelIso pdgId mvaIdSpring15'.split():
                    lepret["Lep"+str(lcount)+"_"+lfloat+self.label] = getattr(lep,lfloat)
                lepvectors.append(lep)
                lepret['metl'+str(lcount)+'DPhi'+self.label] = abs( deltaPhi( getattr(lep, 'phi'), metphi ))
                lepret["Lep"+str(lcount)+"_"+"minTauDR"+self.label] = minDRTau
                ltlvs[lcount-1].SetPtEtaPhiM(lep.pt, lep.eta, lep.phi, 0.0005 if lep.pdgId == 11 else 0.106)
                lcount += 1
                #print 'good lepton', getattr(lep,'pt'), getattr(lep,'eta'), getattr(lep,'phi'), getattr(lep,'pdgId') 
        else:
            ret['nPairLep'] = 0

        mt2 = -1.
        if ret['nPairLep'] == 2:
            l1mt2 = ROOT.reco.Particle.LorentzVector(lepvectors[0].p4().Px(), lepvectors[0].p4().Py(),lepvectors[0].p4().Pz(),lepvectors[0].p4().Energy())
            l2mt2 = ROOT.reco.Particle.LorentzVector(lepvectors[1].p4().Px(), lepvectors[1].p4().Py(),lepvectors[1].p4().Pz(),lepvectors[1].p4().Energy())
            metp4obj = ROOT.reco.Particle.LorentzVector(met*cos(metphi),met*sin(metphi),0,met)
            mt2 = computeMT2(l1mt2, l2mt2, metp4obj)
        ret['mt2'] = mt2
            
        ### Define jets
        ret["iJ"] = []
        # 0. mark each jet as clean
        for j in jetsc+jetsd: j._clean = True
        # set _clean flag of bad jets to False
        for j in jetsc+jetsd:
            if abs(j.eta) > 2.4 or j.pt < 35:
                j._clean = False
                continue
            for l in lepst:
                #lep = leps[l]
                if deltaR(l,j) < 0.4:
                    j._clean = False

        # 2. compute the jet list
        
        for ijc,j in enumerate(jetsc):
            if not j._clean: continue
            ret["iJ"].append(ijc)
        for ijd,j in enumerate(jetsd):
            if not j._clean: continue
            ret["iJ"].append(-1-ijd)
        ret['nJetSel'] = len(ret["iJ"])

        # 3. sort the jets by pt
        
        ret["iJ"].sort(key = lambda idx : jetsc[idx].pt if idx >= 0 else jetsd[-1-idx].pt, reverse = True)

        # 4. compute the variables
        
        for jfloat in "pt eta phi mass btagCSV rawPt".split():
            jetret[jfloat] = []
        if self.isMC:
            for jmc in "mcPt mcFlavour mcMatchId".split():
                jetret[jmc] = []
        for idx in ret["iJ"]:
            jet = jetsc[idx] if idx >= 0 else jetsd[-1-idx]
            for jfloat in "pt eta phi mass btagCSV rawPt".split():
                jetret[jfloat].append( getattr(jet,jfloat) )
            if self.isMC:
                for jmc in "mcPt mcFlavour mcMatchId".split():
                    jetret[jmc].append( getattr(jet,jmc) )
        
        # 5. compute the sums
        
        ret["nJet35"] = 0; ret["htJet35j"] = 0; ret["nBJetLoose35"] = 0; ret["nBJetMedium35"] = 0
        totalRecoil = ROOT.TLorentzVector()
        for j in jetsc+jetsd:
            if not j._clean: continue
            ret["nJet35"] += 1; ret["htJet35j"] += j.pt; 
            if j.btagCSV>0.423: ret["nBJetLoose35"] += 1
            if j.btagCSV>0.890: ret["nBJetMedium35"] += 1
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
            totalRecoil = totalRecoil + jet
          ## compute mlb for the two lepton  
        ret['lepsJZB_recoil'] = totalRecoil.Pt() - ret['lepsZPt']
	
        jet = ROOT.TLorentzVector()
        min_mlb = 1e6
        max_mlb = 1e6
        _lind, _jind = -99, -99
        leplist = [l1, l2]
        # find the global minimum mlb (or mlj)
        jetIsB = False
        for lepton in leplist:
            if ret['nPairLep'] < 2: continue
            for j in jetsc+jetsd:
                if not j._clean: continue
                jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                tmp = (lepton+jet).M()
                if j.btagCSV>0.814:   
                   if tmp < min_mlb: 
                         min_mlb  = tmp
                         jetisB = True 
                         _lind = leplist.index(lepton)
                         _jind = j
                else:
                     if tmp < min_mlb and jetIsB == False:
                         min_mlb = tmp
                         _lind = leplist.index(lepton)
                         _jind = j
        
        # compute the minimum mlb (or mlj) for the other lepton
        jetIsB = False
        for lepton in leplist: 
            if ret['nPairLep'] < 2: continue 
            for j in jetsc+jetsd: 
                if not j._clean: continue
                if j == _jind: continue
                jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                tmp = ( (l1 if _lind == 1 else l2) +jet).M()
                if j.btagCSV>0.814:  
                    if tmp < max_mlb: 
                        max_mlb  = tmp
                        jetIsB = True
                else:
                    if tmp < max_mlb and jetIsB == False:
                        max_mlb = tmp   
        ret["min_mlb1"] = min_mlb if min_mlb < 1e6  else -1.
        ret["min_mlb2"] = max_mlb if max_mlb < 1e6  else -1.
        ret["sum_mlb"] = ret["min_mlb1"] + ret["min_mlb2"]
        ret["st"] = met+lepret["Lep1_pt"+self.label]+lepret["Lep2_pt"+self.label]

        ## beam halo filter list file:
        ## do this only for data
        if not self.isMC:
            evt_str = '%d:%d:%d'%(event.run, event.lumi, event.evt)
            if evt_str in self.beamHaloSet:
                ret['nPairLep'] = -1
            if evt_str in self.fourthBadEESuperCrystalSet:
                ret['nPairLep'] = -1
            if evt_str in self.badResolutionTrackTaggerSet:
                ret['nPairLep'] = -1
            if evt_str in self.badMuonTrackTaggerSet:
                ret['nPairLep'] = -1
        ## ====== done with beam halo and other filters check
        
        ## get the SR id which is 1xx for central and 2xx for forward. the 10 digit is the number of 
        ## b-tags and the signle digit is the mll region going from 1-5
        isBasicSREvent = (ret['nPairLep'] > 0 and ret["lepsDR"] > 0.3 and lepret["Lep1_pt"+self.label] > 20. and lepret["Lep2_pt"+self.label] > 20. and ret['lepsMll'] > 20.)
        isBasicSREvent = isBasicSREvent * (abs(lepret["Lep1_eta"+self.label] - 1.5) > 0.1 and abs(lepret["Lep2_eta"+self.label] - 1.5) > 0.1)
        isBasicSREvent = isBasicSREvent * ((met > 150 and ret['nJetSel'] >= 2 ) or (met > 100. and ret['nJetSel'] >=3))

        if isBasicSREvent:
            srID = self.getSRID(ret['lepsMll'], lepret["Lep1_eta"+self.label], lepret["Lep2_eta"+self.label], ret["nBJetMedium35"])
            ret["srID"] = srID
            for t in ['data', 'mc']:
                for u in ['_ana', '']:
                    ret["lh%s_zpt_%s"%(u,t)] = getattr(self, 'h_lh%s_zpt_%s'%(u,t)).GetBinContent( getattr(self, 'h_lh%s_zpt_%s'%(u,t)).FindBin( ret["lepsZPt"]  ) ) 
                    ret["lh%s_met_%s"%(u,t)] = getattr(self, 'h_lh%s_met_%s'%(u,t)).GetBinContent( getattr(self, 'h_lh%s_met_%s'%(u,t)).FindBin( met             ) ) 
                    ret["lh%s_mlb_%s"%(u,t)] = getattr(self, 'h_lh%s_mlb_%s'%(u,t)).GetBinContent( getattr(self, 'h_lh%s_mlb_%s'%(u,t)).FindBin( ret["sum_mlb"]  ) ) 
                    ret["lh%s_ldr_%s"%(u,t)] = getattr(self, 'h_lh%s_ldr_%s'%(u,t)).GetBinContent( getattr(self, 'h_lh%s_ldr_%s'%(u,t)).FindBin( ret["lepsDR" ]  ) ) 
                    if not ret["lh%s_mlb_%s"%(u,t)]: ret["lh%s_mlb_%s"%(u,t)] = 1e-6
                    if not ret["lh%s_ldr_%s"%(u,t)]: ret["lh%s_ldr_%s"%(u,t)] = 1e-6
                    if not ret["lh%s_met_%s"%(u,t)]: ret["lh%s_met_%s"%(u,t)] = 1e-6
                    if not ret["lh%s_zpt_%s"%(u,t)]: ret["lh%s_zpt_%s"%(u,t)] = 1e-6

                    ret["cum%s_zpt_%s"%(u,t)] = getattr(self, 'h_lh%s_zpt_%s'%(u,t)).Integral( 1, getattr(self, 'h_lh%s_zpt_%s'%(u,t)).FindBin( ret["lepsZPt"]  ) ) 
                    ret["cum%s_met_%s"%(u,t)] = getattr(self, 'h_lh%s_met_%s'%(u,t)).Integral( 1, getattr(self, 'h_lh%s_met_%s'%(u,t)).FindBin( met             ) ) 
                    ret["cum%s_mlb_%s"%(u,t)] = getattr(self, 'h_lh%s_mlb_%s'%(u,t)).Integral( 1, getattr(self, 'h_lh%s_mlb_%s'%(u,t)).FindBin( ret["sum_mlb"]  ) ) 
                    ret["cum%s_ldr_%s"%(u,t)] = getattr(self, 'h_lh%s_ldr_%s'%(u,t)).Integral( 1, getattr(self, 'h_lh%s_ldr_%s'%(u,t)).FindBin( ret["lepsDR" ]  ) ) 

        else:
            ret["srID"]      = -99
            for t in ['data', 'mc']:
                for u in ['_ana', '']:
                    ret["lh%s_mlb_%s"%(u,t)] = -999.
                    ret["lh%s_ldr_%s"%(u,t)] = -999.
                    ret["lh%s_met_%s"%(u,t)] = -999.
                    ret["lh%s_zpt_%s"%(u,t)] = -999.

                    ret["cum%s_mlb_%s"%(u,t)] = -999.
                    ret["cum%s_ldr_%s"%(u,t)] = -999.
                    ret["cum%s_met_%s"%(u,t)] = -999.
                    ret["cum%s_zpt_%s"%(u,t)] = -999.

        
        fullret = {}
        for k,v in ret.iteritems(): 
            fullret[k+self.label] = v
        for k,v in jetret.iteritems(): 
            fullret["JetSel%s_%s" % (self.label,k)] = v
        #for k,v in lepret.iteritems(): 
        #    fullret["Lep%s_%s" % (self.label,k)] = v
        for k,v in lepret.iteritems(): 
            fullret[k] = v
        return fullret

    def getMll_JZB(self, l1, l2, met, met_raw):
        metrecoil = (met + l1 + l2).Pt()
        metrawrecoil = (met_raw + l1 + l2).Pt() 
        zpt = (l1 + l2).Pt()
        jzb = metrecoil - zpt
        jzb_raw = metrawrecoil - zpt
        return (l1+l2).M(), jzb, jzb_raw, l1.DeltaR(l2), metrecoil, zpt, abs( deltaPhi( l1.Phi(), l2.Phi() ) )

    def getPairVariables(self,lepst, metp4, metp4_raw):
        ret = (-999,-999,-99., -9000., -9000, -99., -99., -99., -99.)
        if len(lepst) >= 2:
            [mll, jzb, jzb_raw, dr, metrec, zpt, dphi] = self.getMll_JZB(lepst[0].p4(), lepst[1].p4(), metp4, metp4_raw)
            ret = (0, 1, mll, jzb, jzb_raw, dr, metrec, zpt, dphi)
        return ret

    def getSRID(self, mll, eta1, eta2, nb):
        mllid, bid, etaid = -1, -1, -1
        if    20. < mll <  70.:
            mllid = 1
        elif  70. < mll <  81.:
            mllid = 2
        elif  81. < mll < 101.:
            mllid = 3
        elif 101. < mll < 120.:
            mllid = 4
        elif 120. < mll:
            mllid = 5
            
        if abs(eta1) < 1.4 and abs(eta2) < 1.4:
            etaid = 1
        else:
            etaid = 2

        return (100*etaid + 10*nb + mllid)

##  ## ===============================================================
##  ## ===== bunch of b-tagging stuff. ===============================
##  ## ===============================================================
##      def init_btagMediumScaleFactor(self,CSVbtagFileName,EFFbtagFileName,CSVbtagFileNameFastSim):
##          self.do_btagSF = True
##          self.btagMediumCalib = ROOT.BTagCalibration("CSVv2", CSVbtagFileName)
##          if CSVbtagFileNameFastSim: self.btagMediumCalibFastSim = ROOT.BTagCalibration("CSV_FastSim", CSVbtagFileNameFastSim)
##          self.btagMediumReader=[]
##          self.btagMediumReader.append(ROOT.BTagCalibrationReader(self.btagMediumCalib, 1, "mujets", "down"))
##          self.btagMediumReader.append(ROOT.BTagCalibrationReader(self.btagMediumCalib, 1, "mujets", "central"))
##          self.btagMediumReader.append(ROOT.BTagCalibrationReader(self.btagMediumCalib, 1, "mujets", "up"))
##          if CSVbtagFileNameFastSim:
##              self.btagMediumReaderFastSim=[]
##              self.btagMediumReaderFastSim.append(ROOT.BTagCalibrationReader(self.btagMediumCalibFastSim, 1, "fastsim", "down"))
##              self.btagMediumReaderFastSim.append(ROOT.BTagCalibrationReader(self.btagMediumCalibFastSim, 1, "fastsim", "central"))
##              self.btagMediumReaderFastSim.append(ROOT.BTagCalibrationReader(self.btagMediumCalibFastSim, 1, "fastsim", "up"))
##          self.btagMediumReaderLight=[]
##          self.btagMediumReaderLight.append(ROOT.BTagCalibrationReader(self.btagMediumCalib, 1, "comb", "down"))
##          self.btagMediumReaderLight.append(ROOT.BTagCalibrationReader(self.btagMediumCalib, 1, "comb", "central"))
##          self.btagMediumReaderLight.append(ROOT.BTagCalibrationReader(self.btagMediumCalib, 1, "comb", "up"))
##          if CSVbtagFileNameFastSim:
##              self.btagMediumReaderLightFastSim=[]
##              self.btagMediumReaderLightFastSim.append(ROOT.BTagCalibrationReader(self.btagMediumCalibFastSim, 1, "fastsim", "down"))
##              self.btagMediumReaderLightFastSim.append(ROOT.BTagCalibrationReader(self.btagMediumCalibFastSim, 1, "fastsim", "central"))
##              self.btagMediumReaderLightFastSim.append(ROOT.BTagCalibrationReader(self.btagMediumCalibFastSim, 1, "fastsim", "up"))
##          self.btagEffFile = ROOT.TFile(EFFbtagFileName,"read")
##          self.btagEffHistos = (self.btagEffFile.Get("h2_BTaggingEff_csv_med_Eff_b"),self.btagEffFile.Get("h2_BTaggingEff_csv_med_Eff_c"),self.btagEffFile.Get("h2_BTaggingEff_csv_med_Eff_udsg"))
##      def read_btagMediumScaleFactor(self,readersBC,readersLight,jet,flavor,shift=0,croplowpt=0,crophighpt=1e6):
##          if abs(shift)>2: raise RuntimeError, 'Unsupported b-tag shift value was passed: %d'%shift
##          # agreed upon: include jets in under/overflow in last bins
##          pt = min(max(jet.pt,croplowpt+0.001),crophighpt-0.001)
##          eta = min(max(jet.eta,-2.399),2.399)
##          if abs(flavor)==5: fcode = 0
##          elif abs(flavor)==4: fcode = 1
##          else: fcode = 2
##          res = 0
##          if fcode<2: # correlate systs of B and C
##              _s = shift if abs(shift)<2 else 0
##              res = readersBC[_s+1].eval(fcode,eta,pt)
##          else:
##              _s = shift/2 if abs(shift)!=1 else 0
##              res = readersLight[_s+1].eval(fcode,eta,pt)
##          if res==0: raise RuntimeError,'Btag SF returned zero, something is not correct: flavor=%d, eta=%f, pt=%f'%(flavor,jet.eta,jet.pt)
##          return res
##      def read_btagMediumEfficiency(self,jet,flavor):
##          if abs(flavor)==5: fcode = 0
##          elif abs(flavor)==4: fcode = 1
##          else: fcode = 2
##          h = self.btagEffHistos[fcode]
##          ptbin = max(1,min(h.GetNbinsX(),h.GetXaxis().FindBin(jet.pt)))
##          etabin = max(1,min(h.GetNbinsY(),h.GetYaxis().FindBin(abs(jet.eta))))
##          return h.GetBinContent(ptbin,etabin)
##  
##      def btagMediumScaleFactor(self,event,bjets,alljets,shift=0):
##          if event.isData or (not self.do_btagSF): 
##              return 1.0
##          pmc = 1.0; pdata = 1.0
##          for j in alljets:
##              sf = self.read_btagMediumScaleFactor(self.btagMediumReader,self.btagMediumReaderLight,j,j.mcFlavour,shift if abs(shift)<3 else 0,croplowpt=30,crophighpt=670)
##              if self.isFastSim: 
##                  sf = sf * self.read_btagMediumScaleFactor(self.btagMediumReaderFastSim,self.btagMediumReaderLightFastSim,j,j.mcFlavour,0 if abs(shift)<3 else int(copysign(abs(shift)-2,shift)),croplowpt=20,crophighpt=800)
##              eff = self.read_btagMediumEfficiency(j,j.mcFlavour)
##              if j in bjets:
##                  pmc   = pmc * eff
##                  pdata = pdata * eff * sf
##              else:
##                  pmc = pmc * (1-eff)
##                  pdata = pdata * (1-eff*sf)
##          res = pdata/pmc if pmc!=0 else 1.
##          return res
    
def _susyEdge(lep):
        if lep.pt <= 10.: return False
        if abs(lep.eta) > 2.4: return False
        if abs(lep.dxy) > 0.05: return False
        if abs(lep.dz ) > 0.10: return False
        if abs(lep.eta) > 1.4 and abs(lep.eta) < 1.6: return False
        # marc aug07 if abs(lep.pdgId) == 13 and lep.muonMediumId != 1: return False
        if abs(lep.pdgId) == 13:
          if lep.mediumMuonId != 1: return False
          if lep.miniRelIso > 0.2: return False
        #if abs(lep.pdgId) == 11 and (lep.tightId < 1 or (abs(lep.etaSc) > 1.4442 and abs(lep.etaSc) < 1.566)) : return False
        if abs(lep.pdgId) == 11:
          if (abs(lep.etaSc) > 1.4442 and abs(lep.etaSc) < 1.566) : return False
          if (lep.convVeto == 0) or (lep.lostHits > 0) : return False
          ## phys14 training if (abs(lep.eta) < 0.8 and lep.mvaIdPhys14 < 0.73) : return False
          ## phys14 training if (abs(lep.eta) > 0.8 and abs(lep.eta) < 1.479 and lep.mvaIdPhys14 < 0.57) : return False
          ## phys14 training if (abs(lep.eta) > 1.479 and lep.mvaIdPhys14 < 0.05) : return False
          if (abs(lep.eta) < 0.8 and lep.mvaIdSpring15 < 0.87) : return False
          if (abs(lep.eta) > 0.8 and abs(lep.eta) < 1.479 and lep.mvaIdSpring15 < 0.60) : return False
          if (abs(lep.eta) > 1.479 and lep.mvaIdSpring15 < 0.17) : return False
          if lep.miniRelIso > 0.1: return False
        return True

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf1 = edgeFriends("Edge", 
                lambda lep : _susyEdge(lep),
                cleanJet = lambda lep,jet,dr : (jet.pt < 35 and dr < 0.4 and abs(jet.eta) > 2.4))
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            print self.sf1(ev)
            print self.sf2(ev)
            print self.sf3(ev)
            print self.sf4(ev)
            print self.sf5(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)
