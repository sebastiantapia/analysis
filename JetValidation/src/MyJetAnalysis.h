#ifndef MYJETANALYSIS_H
#define MYJETANALYSIS_H

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>
#include <vector>
#include <utility>  // std::pair, std::make_pair

#include <array>
#include <g4jets/Jet.h>

#include <fastjet/PseudoJet.hh>

#include "/sphenix/user/stapiaar/jet_correlator/EnergyEnergyCorrelators/eec/include/EEC.hh"

class PHCompositeNode;
class JetEvalStack;
class TTree;
class TH1;
class FastJetAlgo;

/// \class MyJetAnalysis
class MyJetAnalysis : public SubsysReco
{
 public:
  MyJetAnalysis(
      const std::string &recojetname = "AntiKt_Tower_r04",
      const std::string &truthjetname = "AntiKt_Truth_r04",
      const std::string &outputfilename = "myjetanalysis.root");

  virtual ~MyJetAnalysis();

  //! set eta range
  void
  setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }
  //! set eta range
  void
  setPtRange(double low, double high)
  {
    m_ptRange.first = low;
    m_ptRange.second = high;
  }
 
  void
  setMindR(double jetradius)
  {
      m_trackJetMatchingRadius = jetradius;
  }

  void doECCprint( bool doprint ) 
	{
		m_do_ecc_print = doprint;
	}

  void use_initial_vertex(const bool b = true) {initial_vertex = b;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void initializeTrees();
  void Clean();

 private:
  //! cache the jet evaluation modules
  std::shared_ptr<JetEvalStack> m_jetEvalStack;
  FastJetAlgo * m_fastjet;
  std::vector<Jet*> m_inputs_smallR;
  std::vector<Jet*> m_inputs_reco_smallR;
  std::vector<Jet*> m_inputs_recoMatch_smallR;
  std::vector<Jet*> m_output_largeR;

  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  std::vector<fastjet::PseudoJet> m_truth_pseudojet;

  //! eta range
  std::pair<double, double> m_etaRange;

  //! pT range
  std::pair<double, double> m_ptRange;

  //! flag to use initial vertex in track evaluator
  bool initial_vertex = false;

  //! max track-jet matching radius
  double m_trackJetMatchingRadius;

  //! Output histograms
  TH1 *m_hInclusiveE;
  TH1 *m_hInclusiveEta;
  TH1 *m_hInclusivePhi;

  //! Output Tree variables
  TTree *m_T;
  TTree *m_Treclus;
  TTree *m_TreEcoor;

  int m_event;
  int m_id;
  int m_nComponent;
  float m_eta;
  float m_phi;
  float m_e;
  float m_pt;
  float m_unsub_pt;

  int m_truthID;
  int m_truthNComponent;
  float m_truthEta;
  float m_truthPhi;
  float m_truthE;
  float m_truthPt;
  float m_truthdR;

bool m_do_ecc_print = true;
fastjet::contrib::eec::EECLongestSideLog *eec_longestside;

//large R
std::vector<double> m_JetLargeR_truth_pt;
std::vector<double> m_JetLargeR_truth_eta;
std::vector<double> m_JetLargeR_truth_phi;
std::vector<double> m_JetLargeR_truth_e;
std::vector<double> m_JetLargeR_truth_d12;
std::vector<int>    m_JetLargeR_truth_nComp;

std::vector<double> m_JetLargeR_reco_pt;
std::vector<double> m_JetLargeR_reco_eta;
std::vector<double> m_JetLargeR_reco_phi;
std::vector<double> m_JetLargeR_reco_e;
std::vector<double> m_JetLargeR_reco_d12;
std::vector<int>    m_JetLargeR_reco_nComp;

//small R
std::vector<std::vector<double>> m_JetSmallR_truth_pt;
std::vector<std::vector<double>> m_JetSmallR_truth_eta;
std::vector<std::vector<double>> m_JetSmallR_truth_phi;
std::vector<std::vector<double>> m_JetSmallR_truth_e;

  //! number of matched tracks
  int m_nMatchedTrack;

  enum
  {
    //! max number of tracks
    kMaxMatchedTrack = 1000
  };
  std::array<float, kMaxMatchedTrack> m_trackdR;
  std::array<float, kMaxMatchedTrack> m_trackpT;
  std::array<float, kMaxMatchedTrack> m_trackPID;
};

#endif  // MYJETANALYSIS_H
