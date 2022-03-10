#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <g4jets/FastJetAlgo.h>
#include <g4jets/JetReco.h>
#include <g4jets/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>

#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>

#include <myjetanalysis/MyJetAnalysis.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libmyjetanalysis.so)

namespace Enable
{
  bool HIJETS = true;
  int HIJETS_VERBOSITY = 1;
}  // namespace Enable

namespace G4HIJETS
{
  bool do_flow = true;
  bool do_CS = false;
}  // namespace G4HIJETS
#endif


void Fun4All_CaloAna(const char *filelisttruth = "dst_truth_test.list", 
		     const char *filelist = "dst_track_test.list", 
                     const char *filelistcalo = "dst_calo_test.list",  
                     const char *outname = "outputest.root")
{

  gSystem->Load("libg4dst");
  gSystem->Load("libmyjetanalysis");

  
  Fun4AllServer *se = Fun4AllServer::instance();
  int verbosity = 0;

  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  //rc->set_IntFlag("RUNNUMBER",1);

  JetReco *truthjetreco = new JetReco();
  TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
  tji->add_embedding_flag(1);  // changes depending on signal vs. embedded
  truthjetreco->add_input(tji);
  truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.2), "AntiKt_Truth_r02");
  truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.3), "AntiKt_Truth_r03");
  truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
  truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.5), "AntiKt_Truth_r05");
  truthjetreco->set_algo_node("ANTIKT");
  truthjetreco->set_input_node("TRUTH");
  truthjetreco->Verbosity(verbosity); 
  se->registerSubsystem(truthjetreco);

  RetowerCEMC *rcemc = new RetowerCEMC(); 
  rcemc->Verbosity(verbosity); 
  se->registerSubsystem(rcemc);

  JetReco *towerjetreco = new JetReco();
  towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER_RETOWER));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWER));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWER));
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.2), "AntiKt_Tower_HIRecoSeedsRaw_r02");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  towerjetreco->Verbosity(verbosity); 
  se->registerSubsystem(towerjetreco);

  DetermineTowerBackground *dtb = new DetermineTowerBackground();
  dtb->SetBackgroundOutputName("TowerBackground_Sub1");
  dtb->SetFlow(G4HIJETS::do_flow);
  dtb->SetSeedType(0);
  dtb->SetSeedJetD(3);
  dtb->Verbosity(verbosity); 
  se->registerSubsystem(dtb);

  CopyAndSubtractJets *casj = new CopyAndSubtractJets();
  casj->SetFlowModulation(G4HIJETS::do_flow);
  casj->Verbosity(verbosity); 
  se->registerSubsystem(casj);

  DetermineTowerBackground *dtb2 = new DetermineTowerBackground();
  dtb2->SetBackgroundOutputName("TowerBackground_Sub2");
  dtb2->SetFlow(G4HIJETS::do_flow);
  dtb2->SetSeedType(1);
  dtb2->SetSeedJetPt(7);
  dtb2->Verbosity(verbosity); 
  se->registerSubsystem(dtb2);

  SubtractTowers *st = new SubtractTowers();
  st->SetFlowModulation(G4HIJETS::do_flow);
  st->Verbosity(verbosity);
  se->registerSubsystem(st);

  towerjetreco = new JetReco();
  towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER_SUB1));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWER_SUB1));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWER_SUB1));
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.2, 1), "AntiKt_Tower_r02_Sub1");
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.3, 1), "AntiKt_Tower_r03_Sub1");
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4, 1), "AntiKt_Tower_r04_Sub1");
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.5, 1), "AntiKt_Tower_r05_Sub1");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);


// My Analysis
  MyJetAnalysis *myJetAnalysis = new MyJetAnalysis("AntiKt_Tower_r02_Sub1", "AntiKt_Truth_r02", "myjet_ppOverlay_TEST.root");
  myJetAnalysis->setPtRange(5, 100);
  myJetAnalysis->setEtaRange(-1.1, 1.1);
  se->registerSubsystem(myJetAnalysis);
  
  Fun4AllInputManager *intrue2 = new Fun4AllDstInputManager("DSTtruth2");
  intrue2->AddListFile(filelisttruth,1);
  se->registerInputManager(intrue2);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTtrack");
  in->AddListFile(filelist,1);
  se->registerInputManager(in);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
  in2->AddListFile(filelistcalo,1);
  se->registerInputManager(in2);

  
  se->run(100);
  se->End();

gSystem->Exit(0);
return 0;

}
