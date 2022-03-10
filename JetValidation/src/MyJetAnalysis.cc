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

  
  for (JetMap::Iter iter = jetsMC->begin(); iter != jetsMC->end(); ++iter)
  {
    Jet* truthjet = iter->second;
   
    //assert(jet);
    bool eta_cut = (truthjet->get_eta() >= m_etaRange.first) and (truthjet->get_eta() <= m_etaRange.second);
    bool pt_cut = (truthjet->get_pt() >= m_ptRange.first) and (truthjet->get_pt() <= m_ptRange.second);
    if ((not eta_cut) or (not pt_cut)) continue;    
    

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
    } //end reco jets

    m_T->Fill();
  }  //   for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)


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

}


