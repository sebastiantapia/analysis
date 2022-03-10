#include "MyJetAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <g4eval/JetEvalStack.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <g4jets/FastJetAlgo.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <g4jets/JetMap.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
//#include <memory> 

//#include "/sphenix/user/stapiaar/jet_correlator/EnergyEnergyCorrelators/eec/include/EEC.hh"

using namespace std;

MyJetAnalysis::MyJetAnalysis(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename)
  : SubsysReco("MyJetAnalysis_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_trackJetMatchingRadius(.7)
  , m_hInclusiveE(nullptr)
  , m_hInclusiveEta(nullptr)
  , m_hInclusivePhi(nullptr)
  , m_T(nullptr)
  , m_event(-1)
  , m_id(-1)
  , m_nComponent(-1)
  , m_eta(numeric_limits<float>::signaling_NaN())
  , m_phi(numeric_limits<float>::signaling_NaN())
  , m_e(numeric_limits<float>::signaling_NaN())
  , m_pt(numeric_limits<float>::signaling_NaN())
  , m_truthID(-1)
  , m_truthNComponent(-1)
  , m_truthEta(numeric_limits<float>::signaling_NaN())
  , m_truthPhi(numeric_limits<float>::signaling_NaN())
  , m_truthE(numeric_limits<float>::signaling_NaN())
  , m_truthPt(numeric_limits<float>::signaling_NaN())
  , m_nMatchedTrack(-1)
{
  m_trackdR.fill(numeric_limits<float>::signaling_NaN());
  m_trackpT.fill(numeric_limits<float>::signaling_NaN());
  m_trackPID.fill(numeric_limits<float>::signaling_NaN());
}

MyJetAnalysis::~MyJetAnalysis()
{
}

int MyJetAnalysis::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::Init - Outoput to " << m_outputFileName << endl;

  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  cout << "MyJetAnalysis::Init - Outoput to " << m_outputFileName << endl;

  // Histograms
  m_hInclusiveE = new TH1F(
      "hInclusive_E",  //
      TString(m_recoJetName) + " inclusive jet E;Total jet energy (GeV)", 100, 0, 100);

  m_hInclusiveEta =
      new TH1F(
          "hInclusive_eta",  //
          TString(m_recoJetName) + " inclusive jet #eta;#eta;Jet energy density", 50, -1, 1);
  m_hInclusivePhi =
      new TH1F(
          "hInclusive_phi",  //
          TString(m_recoJetName) + " inclusive jet #phi;#phi;Jet energy density", 50, -M_PI, M_PI);

  initializeTrees();

  //eec_longestside(2, 75, {1e-5, 1});
  //fastjet::contrib::eec::EECLongestSideLog eec_longestside(2, 75, {1e-5, 1});
  eec_longestside = new fastjet::contrib::eec::EECLongestSideLog(2, 75, {1e-5, 1});
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::End(PHCompositeNode* topNode)
{
  cout << "MyJetAnalysis::End - Output to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);

  m_hInclusiveE->Write();
  m_hInclusiveEta->Write();
  m_hInclusivePhi->Write();
  m_T->Write();
  m_Treclus->Write();
  //m_TreEcoor->Write();
  if(m_do_ecc_print) cout << *eec_longestside;

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::InitRun(PHCompositeNode* topNode)
{
  cout << "MyJetAnalysis::InitRun" << endl;
  m_jetEvalStack = shared_ptr<JetEvalStack>(new JetEvalStack(topNode, m_recoJetName, m_truthJetName));
  m_jetEvalStack->get_stvx_eval_stack()->set_use_initial_vertex(initial_vertex);
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyJetAnalysis::process_event(PHCompositeNode* topNode)
{
  cout << "MyJetAnalysis::process_event" << endl;
  if (Verbosity() >= MyJetAnalysis::VERBOSITY_SOME)
    cout << "MyJetAnalysis::process_event() entered" << endl;

  //m_jetEvalStack->next_event(topNode);
  //JetRecoEval* recoeval = m_jetEvalStack->get_reco_eval();
  ++m_event;

  if( (m_event % 10) == 0 ) cout << "Event number = "<< m_event << endl;


  // interface to jets
  JetMap* jets = findNode::getClass<JetMap>(topNode, m_recoJetName);
  if (!jets)
  {
    cout
        << "MyJetAnalysis::process_event - Error can not find DST JetMap node "
        << m_recoJetName << endl;
    exit(-1);
  }

  JetMap* jetsMC = findNode::getClass<JetMap>(topNode, m_truthJetName);
  if (!jetsMC)
  {
    cout
        << "MyJetAnalysis::process_event - Error can not find DST Truth JetMap node "
        << m_truthJetName << endl;
    exit(-1);
  }

  // interface to tracks
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
  {
    trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
    if (!trackmap)
    {
      cout
          << "MyJetAnalysis::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
      exit(-1);
    }
  }

 PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
 {
   cout << " ERROR: Can't find G4TruthInfo" << endl;
   exit(-1);
 }

  
  m_inputs_smallR.clear();
  for (JetMap::Iter iter = jetsMC->begin(); iter != jetsMC->end(); ++iter)
  {
    Jet* truthjet = iter->second;
   
    //assert(jet);
    bool eta_cut = (truthjet->get_eta() >= m_etaRange.first) and (truthjet->get_eta() <= m_etaRange.second);
    bool pt_cut = (truthjet->get_pt() >= m_ptRange.first) and (truthjet->get_pt() <= m_ptRange.second);
    if ((not eta_cut) or (not pt_cut)) continue;    
    
    m_inputs_smallR.push_back(truthjet);


/*    if ((not eta_cut) or (not pt_cut))
    {
      if (Verbosity() >= MyJetAnalysis::VERBOSITY_MORE)
      {
        cout << "MyJetAnalysis::process_event() - jet failed acceptance cut: ";
        cout << "eta cut: " << eta_cut << ", ptcut: " << pt_cut << endl;
        cout << "jet eta: " << jet->get_eta() << ", jet pt: " << jet->get_pt() << endl;
        jet->identify();
      }
      continue;
    }

    // fill histograms
    assert(m_hInclusiveE);
    m_hInclusiveE->Fill(jet->get_e());
    assert(m_hInclusiveEta);
    m_hInclusiveEta->Fill(jet->get_eta());
    assert(m_hInclusivePhi);
    m_hInclusivePhi->Fill(jet->get_phi());
*/
    // fill trees - jet spectrum
    //Jet* truthjet = recoeval->max_truth_jet_by_energy(jet);

    m_id = -1;
    m_nComponent = -1;
    m_eta = NAN;
    m_phi = NAN;
    m_e = NAN;
    m_pt = NAN;

    m_truthID = truthjet->get_id();
    m_truthNComponent = truthjet->size_comp();
    m_truthEta = truthjet->get_eta();
    m_truthPhi = truthjet->get_phi();
    m_truthE = truthjet->get_e();
    m_truthPt = truthjet->get_pt();
    m_truthdR = 10;

  m_nMatchedTrack = 0;


  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
     for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
          iter != range.second;
          ++iter)
     {
       PHG4Particle *part = iter->second;
   
	//KEEP WORKING on THIS to reject embed ( HIJING ) particles
         //const int this_embed_id = truthinfo->isEmbeded(part->get_track_id());
         //cout << "embed id: " << this_embed_id << endl;   

//remove some particles (muons, taus, neutrinos)...
//12 == nu_e
//13 == muons
//14 == nu_mu
//15 == taus
//16 == nu_tau

       if ((abs(part->get_pid()) >= 12) && (abs(part->get_pid()) <= 16)) continue;

	TVector3 v(part->get_px(), part->get_py(), part->get_pz());
        const double dEta = v.Eta() - m_truthEta;
        const double dPhi = v.Phi() - m_truthPhi;
        const double dR = sqrt(dEta * dEta + dPhi * dPhi);

	if(dR>m_trackJetMatchingRadius) continue;
	//cout << "truth conts PID, pt: " << v.Pt() << ", " << part->get_pid() << endl;
	
	m_trackdR[m_nMatchedTrack] = dR;
	m_trackPID[m_nMatchedTrack] = part->get_pid();
        m_trackpT[m_nMatchedTrack] = v.Pt();

        ++m_nMatchedTrack;

     }//end truth particles
    
    //if (truthjet)
    for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
    {
	Jet* jet = iter->second;

	if(jet->get_pt() < 1) continue; // to remove noise jets
	TVector3 v(jet->get_px(), jet->get_py(), jet->get_pz());

	const double dEta = v.Eta() - m_truthEta;
	const double dPhi = v.Phi() - m_truthPhi;
	const double dR = sqrt(dEta * dEta + dPhi * dPhi);
	if(m_truthdR<dR) continue;
      m_id = jet->get_id();
      m_nComponent = jet->size_comp();
      m_eta = jet->get_eta();
      m_phi = jet->get_phi();
      m_e = jet->get_e();
      m_pt = jet->get_pt();
      m_truthdR = dR;
      if(jet->has_property(Jet::PROPERTY::prop_SeedD)) m_unsub_pt = jet->get_property(Jet::PROPERTY::prop_SeedD);
    } //end reco jets

   // m_inputs_reco_smallR.clear();
   // for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
   // {
   //      Jet* jet = iter->second;
   //      if(jet->get_pt()<10) continue;
  ///	 m_inputs_reco_smallR.push_back(jet);
  //  }

   // fill trees - jet track matching
    //m_nMatchedTrack = 0;
/*
    for (SvtxTrackMap::Iter iter = trackmap->begin();
         iter != trackmap->end();
         ++iter)
    {
      SvtxTrack* track = iter->second;

      TVector3 v(track->get_px(), track->get_py(), track->get_pz());
      const double dEta = v.Eta() - m_eta;
      const double dPhi = v.Phi() - m_phi;
      const double dR = sqrt(dEta * dEta + dPhi * dPhi);

      if (dR < m_trackJetMatchingRadius)
      {
        //matched track to jet

        assert(m_nMatchedTrack < kMaxMatchedTrack);

        m_trackdR[m_nMatchedTrack] = dR;
        m_trackpT[m_nMatchedTrack] = v.Perp();

        ++m_nMatchedTrack;
      }

      if (m_nMatchedTrack >= kMaxMatchedTrack)
      {
        cout << "MyJetAnalysis::process_event() - reached max track that matching a jet. Quit iterating tracks" << endl;
        break;
      }

    }  //    for (SvtxTrackMap::Iter iter = trackmap->begin();
*/
    m_T->Fill();
  }  //   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)

     m_inputs_reco_smallR.clear();
     for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
     {
          Jet* jet = iter->second;
          if(jet->get_pt()<10) continue;
          m_inputs_reco_smallR.push_back(jet);
     }

// fastjet standalone

     int pseudo_idx = 0;
     std::vector<fastjet::PseudoJet> pseudoTruthjets;
     PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
     for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
          iter != range.second;
          ++iter)
     {
       PHG4Particle *part = iter->second; 

       if ((abs(part->get_pid()) >= 12) && (abs(part->get_pid()) <= 16)) continue;

       if ((part->get_px() == 0.0) && (part->get_py() == 0.0)) continue;  // avoid pt=0
       if (part->get_e() == 0.) continue;
	float eta = asinh(part->get_pz() / sqrt(pow(part->get_px(), 2) + pow(part->get_py(), 2)));
       if (eta < -1.) continue;
       if (eta > 1.) continue;
       if(sqrt(pow(part->get_px(),2)+pow(part->get_py(),2)) < 1.) continue; // is not the best cut, Ideally should remove truth particles with small life time

       //cout << "PID particle: " << part->get_pid() << endl;
       //cout << "pT particle: " << sqrt(pow(part->get_px(),2)+pow(part->get_py(),2)) << endl;

             fastjet::PseudoJet pseudojet(part->get_px(),
                                          part->get_py(),
                                          part->get_pz(),
                                          part->get_e());
             pseudojet.set_user_index(pseudo_idx);
             pseudoTruthjets.push_back(pseudojet);       
             ++pseudo_idx;

     } 

   //cout << "N truth particles " << pseudoTruthjets.size() << endl;

   fastjet::JetDefinition jetDef_ecor(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
   fastjet::ClusterSequence jetFinder_ecor(pseudoTruthjets, jetDef_ecor);

   std::vector<fastjet::PseudoJet> fastjets_ecor = jetFinder_ecor.inclusive_jets(5);

   //cout << "N Jets for E correlation: " << fastjets_ecor.size() << endl;

   //fastjet::contrib::eec::EECLongestSideLog eec_longestside(2, 75, {1e-5, 1});

   for (int j = 0; j < (int)fastjets_ecor.size(); ++j)// loop over reclustered large R jets
   {
         vector<fastjet::PseudoJet> consts = jetFinder_ecor.constituents(fastjets_ecor[j]);
         //cout << "truth jets for Ecoor pT: " << fastjets_ecor[j].perp() << endl; 
         //cout << "N of cont: " << consts.size() << endl;

         eec_longestside->compute(consts);

         //cout << "Number of events: " << eec_longestside->event_count() << endl;
         //cout << "Sum of " << eec_longestside->compname() << ": " << eec_longestside->sum() << endl;

   } 

   //cout << *eec_longestside;

   std::vector<fastjet::PseudoJet> pseudojets;
   for (int j = 0; j < (int)m_inputs_smallR.size(); ++j)
   {
	if (m_inputs_smallR.at(j)->get_e() == 0.) continue;
	if (m_inputs_smallR.at(j)->get_pt() < 15.) continue;
		fastjet::PseudoJet pseudojet(m_inputs_smallR.at(j)->get_px(),
                                     	     m_inputs_smallR.at(j)->get_py(),
                                   	     m_inputs_smallR.at(j)->get_pz(),
                                 	     m_inputs_smallR.at(j)->get_e());
		pseudojet.set_user_index(j);
		pseudojets.push_back(pseudojet);
     //cout << "filling small jet in fastjet N: " << j << endl;
   }


    std::vector<fastjet::PseudoJet> pseudojets_reco;
    for (int j = 0; j < (int)m_inputs_reco_smallR.size(); ++j)
    {
    	 if (m_inputs_reco_smallR.at(j)->get_e() == 0.) continue;
         if (m_inputs_reco_smallR.at(j)->get_pt() < 10.) continue;
             fastjet::PseudoJet pseudojet(m_inputs_reco_smallR.at(j)->get_px(),
                                          m_inputs_reco_smallR.at(j)->get_py(),
                                          m_inputs_reco_smallR.at(j)->get_pz(),
                                          m_inputs_reco_smallR.at(j)->get_e());
             pseudojet.set_user_index(j);
             pseudojets_reco.push_back(pseudojet);
    }


   fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, 0.8, fastjet::E_scheme, fastjet::Best);// JEt definition
   fastjet::ClusterSequence jetFinder(pseudojets, jetDef);// vector of jet + jet definition
   std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();// This are the Jets find by FastJet <--- TRUTH

   fastjet::ClusterSequence jetFinder_reco(pseudojets_reco, jetDef);
   std::vector<fastjet::PseudoJet> fastjets_reco = jetFinder_reco.inclusive_jets(); // This are the Jets find by FastJet <--- RECO

   for (int j = 0; j < (int)fastjets.size(); ++j)// loop over reclustered large R jets
   {
	vector<fastjet::PseudoJet> consts = jetFinder.constituents(fastjets[j]);
        //for (int i = 0; i < (int)consts.size(); ++i) cout << "Small-R const pT: " << consts[i].perp() << endl;       
 
        fastjet::JetDefinition jetDef2(fastjet::kt_algorithm, 1.0, fastjet::E_scheme, fastjet::Best);
        fastjet::ClusterSequence kt_clust_seq(consts, jetDef2);
        //std::vector<fastjet::PseudoJet> fastjets_kt = kt_clust_seq.inclusive_jets();
        fastjet::PseudoJet kt_jet = fastjet::sorted_by_pt(kt_clust_seq.inclusive_jets()).front();
        double d12 = 1.0*sqrt(kt_clust_seq.exclusive_subdmerge(kt_jet, 1));

	m_JetLargeR_truth_pt.push_back(fastjets[j].perp());
	m_JetLargeR_truth_eta.push_back(fastjets[j].eta());
	m_JetLargeR_truth_phi.push_back(fastjets[j].phi());
	m_JetLargeR_truth_e.push_back(fastjets[j].e());
	m_JetLargeR_truth_d12.push_back(d12);
	m_JetLargeR_truth_nComp.push_back((int)consts.size());

        //for (int i = 0; i < (int)fastjets_kt.size(); ++i) cout << "Large-R kt pT: " << fastjets_kt[i].perp() << " d12: " << 1.0*sqrt(jetFinder2.exclusive_subdmerge(fastjets_kt[i], 1)) <<endl;
	for (int i = 0; i < (int)fastjets_reco.size(); ++i)
	{

		const double dEta = fastjets_reco[i].eta() - fastjets[j].eta();
		const double dPhi = fastjets_reco[i].phi() - fastjets[j].phi();
		const double dR = sqrt(dEta * dEta + dPhi * dPhi);

		if( dR > 0.6 ) continue;

                vector<fastjet::PseudoJet> consts_reco = jetFinder_reco.constituents(fastjets_reco[i]);
		fastjet::JetDefinition jetDef3(fastjet::kt_algorithm, 1.0, fastjet::E_scheme, fastjet::Best);
                fastjet::ClusterSequence kt_clust_seq_reco(consts_reco, jetDef3);
                fastjet::PseudoJet kt_reco_jet = fastjet::sorted_by_pt(kt_clust_seq_reco.inclusive_jets()).front();
                double d12_reco = 1.0*sqrt(kt_clust_seq_reco.exclusive_subdmerge(kt_reco_jet, 1));
		
		m_JetLargeR_reco_pt.push_back(fastjets_reco[i].perp());
    	        m_JetLargeR_reco_eta.push_back(fastjets_reco[i].eta());
  	        m_JetLargeR_reco_phi.push_back(fastjets_reco[i].phi());
   	        m_JetLargeR_reco_e.push_back(fastjets_reco[i].e());
       		m_JetLargeR_reco_d12.push_back(d12_reco);
	        m_JetLargeR_reco_nComp.push_back((int)consts_reco.size());
	}
   }

   //cout << "vector of large-recluster Truth size fastjet: " << fastjets.size() << endl;

 m_Treclus->Fill();

  Clean();

  //delete m_fastjet;

  return Fun4AllReturnCodes::EVENT_OK;
}



void MyJetAnalysis::initializeTrees()
{

m_T = new TTree("T", "MyJetAnalysis Tree");
m_T->Branch("m_event", &m_event, "event/I");
m_T->Branch("id", &m_id, "id/I");
m_T->Branch("nComponent", &m_nComponent, "nComponent/I");
m_T->Branch("eta", &m_eta, "eta/F");
m_T->Branch("phi", &m_phi, "phi/F");
m_T->Branch("e", &m_e, "e/F");
m_T->Branch("pt", &m_pt, "pt/F");
m_T->Branch("pt_unsub", &m_unsub_pt, "pt/F");
m_T->Branch("truthID", &m_truthID, "truthID/I");
m_T->Branch("truthNComponent", &m_truthNComponent, "truthNComponent/I");
m_T->Branch("truthEta", &m_truthEta, "truthEta/F");
m_T->Branch("truthPhi", &m_truthPhi, "truthPhi/F");
m_T->Branch("truthE", &m_truthE, "truthE/F");
m_T->Branch("truthPt", &m_truthPt, "truthPt/F");
m_T->Branch("truthdR", &m_truthdR, "truthdR/F");
m_T->Branch("nMatchedTrack", &m_nMatchedTrack, "nMatchedTrack/I");
m_T->Branch("AssoParPID", m_trackPID.data(), "trackPID[nMatchedTrack]/F");
m_T->Branch("AssoPardR", m_trackdR.data(), "trackdR[nMatchedTrack]/F");
m_T->Branch("AssoParPt", m_trackpT.data(), "trackpT[nMatchedTrack]/F");


m_Treclus = new TTree("Treclus", "MyJetAnalysis recluster");
m_Treclus->Branch("m_event", &m_event, "event/I");
m_Treclus->Branch("JetLargeR_truth_pt", &m_JetLargeR_truth_pt);
m_Treclus->Branch("JetLargeR_truth_eta", &m_JetLargeR_truth_eta);
m_Treclus->Branch("JetLargeR_truth_phi", &m_JetLargeR_truth_phi);
m_Treclus->Branch("JetLargeR_truth_e", &m_JetLargeR_truth_e);
m_Treclus->Branch("JetLargeR_truth_d12", &m_JetLargeR_truth_d12);
m_Treclus->Branch("JetLargeR_truth_nComp", &m_JetLargeR_truth_nComp);

m_Treclus->Branch("JetLargeR_reco_pt", &m_JetLargeR_reco_pt);
m_Treclus->Branch("JetLargeR_reco_eta", &m_JetLargeR_reco_eta);
m_Treclus->Branch("JetLargeR_reco_phi", &m_JetLargeR_reco_phi);
m_Treclus->Branch("JetLargeR_reco_e", &m_JetLargeR_reco_e);
m_Treclus->Branch("JetLargeR_reco_d12", &m_JetLargeR_reco_d12);
m_Treclus->Branch("JetLargeR_reco_nComp", &m_JetLargeR_reco_nComp);

//m_TreEcoor = new TTree("m_TreEcoor", "MyJetAnalysis Energy correlator");
//m_TreEcoor->Branch("truth_pseudojet", &m_truth_pseudojet);


}

 void MyJetAnalysis::Clean()
{

m_JetLargeR_truth_pt.clear();
m_JetLargeR_truth_eta.clear();
m_JetLargeR_truth_phi.clear();
m_JetLargeR_truth_e.clear();
m_JetLargeR_truth_d12.clear();
m_JetLargeR_truth_nComp.clear();

m_JetLargeR_reco_pt.clear();
m_JetLargeR_reco_eta.clear();
m_JetLargeR_reco_phi.clear();
m_JetLargeR_reco_e.clear();
m_JetLargeR_reco_d12.clear();
m_JetLargeR_reco_nComp.clear();

}







