#include <TChain.h>
#include "RestFrames/RestFrames.hh"

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

// FrameWork includes
#include "AsgTools/ToolHandle.h"
#include "AsgTools/AsgTool.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticleHelpers.h"

#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETAssociationMap.h"
#include "xAODMissingET/MissingETContainer.h"

#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODTau/TauJetContainer.h"

#include "assert.h"
#include "TFile.h"
#include "TRandom3.h"

#include "JetCalibTools/JetCalibrationTool.h"



#include "METUtilities/METMaker.h"

// For RJigsaw Frame Objects
// Background-like tree
RestFrames::LabRecoFrame * LAB_B;
RestFrames::SelfAssemblingRecoFrame * S_B;
RestFrames::VisibleRecoFrame * V_B;
RestFrames::InvisibleRecoFrame * I_B;
RestFrames::InvisibleGroup * INV_B;
RestFrames::CombinatoricGroup * VIS_B;
RestFrames::SetMassInvJigsaw * MinMass_B;
RestFrames::SetRapidityInvJigsaw * Rapidity_B;

// Signal-like tree
RestFrames::LabRecoFrame * LAB_R;
RestFrames::DecayRecoFrame * GG_R;
RestFrames::DecayRecoFrame * Ga_R;
RestFrames::DecayRecoFrame * Gb_R;
RestFrames::DecayRecoFrame * Ca_R;
RestFrames::DecayRecoFrame * Cb_R;
RestFrames::VisibleRecoFrame * V1a_R;
RestFrames::VisibleRecoFrame * V2a_R;
RestFrames::InvisibleRecoFrame * Xa_R;
RestFrames::VisibleRecoFrame * V1b_R;
RestFrames::VisibleRecoFrame * V2b_R;
RestFrames::InvisibleRecoFrame * Xb_R;
RestFrames::InvisibleGroup * INV_R;
RestFrames::CombinatoricGroup * VIS_R;
RestFrames::SetMassInvJigsaw * MinMassJigsaw_R;
RestFrames::SetRapidityInvJigsaw * RapidityJigsaw_R;
RestFrames::ContraBoostInvJigsaw * ContraBoostJigsaw_R;
RestFrames::MinMassesCombJigsaw * HemiJigsaw_R;
RestFrames::MinMassesCombJigsaw * CaHemiJigsaw_R;
RestFrames::MinMassesCombJigsaw * CbHemiJigsaw_R;

met::METMaker * metMaker;

int setupRestFrames();
int calculateMET(double & metx,
		 double & mety
		  );
int calculateRJigsawVariables(const std::vector<xAOD::Jet>& jets,
			     Double_t metx,
			     Double_t mety,
			     std::map<TString,float>& RJigsawVariables,
			     Double_t jetPtCut);


int segFaultRestFrames(){



  std::vector<std::string> sampleName = {
    "/data/users/rsmith/mc15_13TeV.361060.Sherpa_CT10_SinglePhotonPt4000_CVetoBVeto.merge.DAOD_TRUTH1.e3587_p2375/DAOD_TRUTH1.05969136._000001.pool.root.1"
  };

  setupRestFrames();

  for(int nsample = 0; nsample < sampleName.size(); ++nsample){
    //this test file should work.  Feel free to contact me if there is a problem with the file.
    TString const fileName = sampleName.at(nsample);
    std::auto_ptr< TFile > ifile( TFile::Open( fileName, "READ" ) );
    assert( ifile.get() );

      // Create a TEvent object.
    xAOD::TEvent * event = new xAOD::TEvent( xAOD::TEvent::kClassAccess );
    assert( event->readFrom( ifile.get(), xAOD::TEvent::kClassAccess  ) );

    // Create a transient object store. Needed for the tools.
    xAOD::TStore store;

    metMaker= new met::METMaker("maker");

    for(Long64_t ientry = 0; ientry < event.getEntries(); ++ientry){

      xAOD::JetContainer const * jets = nullptr;
      assert( metMaker->evtStore()->retrieve(jets , "AntiKt4EMTopoJets") );

      xAOD::JetContainer* goodJets = new xAOD::JetContainer();
      xAOD::AuxContainerBase* goodJetsAux = new xAOD::AuxContainerBase();
      goodJets->setStore( goodJetsAux ); //< Connect the two

      double const jetptcut = 20000.;

      for(xAOD::JetContainer::const_iterator it = jets->begin();
	  it != jets->end();
	  ++it
	  ){//selection from zero lepton code

	if ( (*it)->pt() <= jetptcut ) continue;
	//   if ( (*it)->auxdecor<char>("passOR") == 0) continue;
	//   if ( (*it)->auxdecor<char>("bad") == 0  ) {
	if ( std::abs((*it)->eta()) < 2.8 ) {

	  xAOD::Jet* jet = new xAOD::Jet();
	  jet->makePrivateStore( **it );
	  goodJets->push_back( jet );
	}
      }

      std::map<TString, float > RJigsawVariables;

      double metx = 0;
      double mety = 0;

      calculateMET(metx,mety);

      calculateRJigsawVaribles(goodjets,
			       metx,
			       mety,
			       RJigsawVariables,
			       jetptcut
			       );


    }//tfile
  }


  return 0;
}

int calculateMET(double & metx,
		 double & mety
		 ){
    //retrieve the original containers
    const xAOD::MissingETContainer* coreMet  = nullptr;
    std::string coreMetKey = "MET_Core_AntiKt4EMTopo";
    assert( metMaker->evtStore()->retrieve(coreMet, coreMetKey) );

    const xAOD::ElectronContainer* electrons = nullptr;
    assert( metMaker->evtStore()->retrieve(electrons, "Electrons") );

    const xAOD::MuonContainer* muons = nullptr;
    assert( metMaker->evtStore()->retrieve(muons, "Muons") );

    const xAOD::PhotonContainer* photons = nullptr;
    assert( metMaker->evtStore()->retrieve(photons, "Photons"));

    const xAOD::TauJetContainer* taus = nullptr;
    assert( metMaker->evtStore()->retrieve(taus, "TauJets"));

    const xAOD::JetContainer* jets = nullptr;
    assert( metMaker->evtStore()->retrieve(jets, "AntiKt4EMTopoJets"));//this retrieves and applies the correction

    //retrieve the MET association map
    const xAOD::MissingETAssociationMap* metMap = nullptr;
    std::string metAssocKey = "METAssoc_AntiKt4EMTopo";
    assert( metMaker->evtStore()->retrieve(metMap, metAssocKey) );

    xAOD::MissingETContainer*    newMetContainer    = new xAOD::MissingETContainer();
    xAOD::MissingETAuxContainer* newMetAuxContainer = new xAOD::MissingETAuxContainer();
    newMetContainer->setStore(newMetAuxContainer);

    // It is necessary to reset the selected objects before every MET calculation
    metMap->resetObjSelectionFlags();

    //here we apply some basic cuts and rebuild the met at each step
    //Electrons
    if(!electrons->empty()){
      ConstDataVector<xAOD::ElectronContainer> metElectrons(SG::VIEW_ELEMENTS);
      for(const auto& el : *electrons) {
	//if(CutsMETMaker::accept(el)) metElectrons.push_back(el);
      }
      assert(metMaker->rebuildMET("RefEle",                   //name of metElectrons in metContainer
				 xAOD::Type::Electron,       //telling the rebuilder that this is electron met
				 newMetContainer,            //filling this met container
				 metElectrons.asDataVector(),//using these metElectrons that accepted our cuts
				 metMap)                     //and this association map
	     );
    }//Photons
    if(!photons->empty()){
      ConstDataVector<xAOD::PhotonContainer> metPhotons(SG::VIEW_ELEMENTS);
      for(const auto& ph : *photons) {
	//if(CutsMETMaker::accept(ph)) metPhotons.push_back(ph);
      }
      assert(metMaker->rebuildMET("RefPhoton",
				 xAOD::Type::Photon,
				 newMetContainer,
				 metPhotons.asDataVector(),
				 metMap)
	     );
    }//Taus
    if(!taus->empty()){
      ConstDataVector<xAOD::TauJetContainer> metTaus(SG::VIEW_ELEMENTS);
      for(const auto& tau : *taus) {
	//if(CutsMETMaker::accept(tau)) metTaus.push_back(tau);
      }
      assert(metMaker->rebuildMET("RefTau",
				 xAOD::Type::Tau,
				 newMetContainer,
				 metTaus.asDataVector(),
				 metMap)
	     );
    }
    //Muons
    if(!muons->empty()){
      ConstDataVector<xAOD::MuonContainer> metMuons(SG::VIEW_ELEMENTS);
      for(const auto& mu : *muons) {
	//if(CutsMETMaker::accept(mu)) metMuons.push_back(mu);
      }
      assert(metMaker->rebuildMET("RefMuon",
				 xAOD::Type::Muon,
				 newMetContainer,
				 metMuons.asDataVector(),
				 metMap)
	     );
    }

    //Now time to rebuild jetMet and get the soft term
    //This adds the necessary soft term for both CST and TST
    //these functions create an xAODMissingET object with the given names inside the container
    assert(  metMaker->rebuildJetMET("RefJet",        //name of jet met
				    "SoftClus",      //name of soft cluster term met
				    "PVSoftTrk",     //name of soft track term met
				    newMetContainer, //adding to this new met container
				    jets,       //using this jet collection to calculate jet met
				    coreMet,         //core met container
				    metMap,          //with this association map
				    false            //don't apply jet jvt cut
				    )
	     );


    assert( metMaker->buildMETSum("FinalTrk" , newMetContainer, MissingETBase::Source::Track ) );

    xAOD::MissingETContainer::const_iterator met_it = newMetContainer->find("FinalTrk");

    metx = (*met_it)->mpx();
    mety = (*met_it)->mpy();

    return 0;
}

int calculateRJigsawVariables(const std::vector<xAOD::Jet>& jets,
			     Double_t metx,
			     Double_t mety,
			     std::map<TString,float>& RJigsawVariables,
			     Double_t jetPtCut)
{
  LAB_R->ClearEvent();
  LAB_B->ClearEvent();


  vector<RestFrames::RFKey> jetID_R;                    // ID for tracking jets in tree

  //std::cout << "number of jets is " << jets.size() << std::endl;

  // Still need to add jets to frames ///////////////
  std::vector<TLorentzVector> myjets;
  for(size_t ijet=0; ijet<jets.size(); ijet++)
    {
      TLorentzVector jet;
      jet.SetPtEtaPhiM(jets[ijet].pt(),
		       jets[ijet].eta(),
		       jets[ijet].phi(),
		       jets[ijet].m());
      myjets.push_back(jet);
    }

  for(size_t ijet=0; ijet<jets.size(); ijet++)
    {
      if(myjets[ijet].Pt()<jetPtCut) continue;
      jetID_R.push_back( VIS_R->AddLabFrameFourVector( myjets[ijet] )  );
      TLorentzVector temp = myjets[ijet];
      temp.SetPtEtaPhiM(temp.Pt(),0.,temp.Phi(),temp.M());
      VIS_B->AddLabFrameFourVector( temp );
    }


  if(jetID_R.size() < 2){
    RJigsawVariables = std::map<TString, float>();
    return -1;
  }


  TVector3 MET_TV3;

  MET_TV3.SetZ(0.);
  MET_TV3.SetX(metx+.01);
  MET_TV3.SetY(mety+.01);


  // std::cout << "njets calc " << jetID_R.size() << std::endl;
  // for(size_t ijet=0; ijet<myjets.size(); ijet++){
  //   myjets.at(ijet).Print();
  // }

  // std::cout << "myjets " << jetID_R.size() << std::endl;
  // std::cout << "METx " << metx << std::endl;
  // std::cout << "METy " << mety << std::endl;

  INV_B->SetLabFrameThreeVector(MET_TV3);
  LAB_B->AnalyzeEvent();

  INV_R->SetLabFrameThreeVector(MET_TV3);
  LAB_R->AnalyzeEvent();



  ////////////////////////////////////////////////////////////////////////////////
  // 1st order vars


  RJigsawVariables[ "RJVars_PP_Mass"           ] = GG_R->GetMass();
  RJigsawVariables[ "RJVars_PP_InvGamma"       ] = 1./GG_R->GetGammaInParentFrame();
  RJigsawVariables[ "RJVars_PP_dPhiBetaR"      ] = GG_R->GetDeltaPhiBoostVisible();
  RJigsawVariables[ "RJVars_PP_dPhiVis"        ] = GG_R->GetDeltaPhiVisible();
  RJigsawVariables[ "RJVars_PP_CosTheta"       ] = GG_R->GetCosDecayAngle();
  RJigsawVariables[ "RJVars_PP_dPhiDecayAngle" ] = GG_R->GetDeltaPhiDecayAngle();
  RJigsawVariables[ "RJVars_PP_VisShape"       ] = GG_R->GetVisibleShape();
  RJigsawVariables[ "RJVars_PP_MDeltaR"        ] = GG_R->GetVisibleShape() * GG_R->GetMass() ;
  RJigsawVariables[ "RJVars_P1_Mass"           ] = Ga_R->GetMass();
  RJigsawVariables[ "RJVars_P1_CosTheta"       ] = Ga_R->GetCosDecayAngle();
  RJigsawVariables[ "RJVars_P2_Mass"           ] = Gb_R->GetMass();
  RJigsawVariables[ "RJVars_P2_CosTheta"       ] = Gb_R->GetCosDecayAngle();
  RJigsawVariables[ "RJVars_I1_Depth"          ] = Ga_R->GetFrameDepth(*Xa_R);
  RJigsawVariables[ "RJVars_I2_Depth"          ] = Gb_R->GetFrameDepth(*Xb_R);

  // end
  //////////////////////////////////////////////////////////////////////////////////





  ////////////////////////////////////////////////////////////////////////////////
  // 2nd order "gluino-like" vars

  RestFrames::DecayRecoFrame* G[2];
  RestFrames::DecayRecoFrame* C[2];
  RestFrames::VisibleRecoFrame* VS[2];
  RestFrames::VisibleRecoFrame* VC[2];
  RestFrames::InvisibleRecoFrame* X[2];
  // Randomize the two hemispheres
  TRandom3 random;

  int flip = (random.Rndm() > 0.5);
  G[flip] = Ga_R;
  G[(flip+1)%2] = Gb_R;
  C[flip] = Ca_R;
  C[(flip+1)%2] = Cb_R;
  VS[flip] = V1a_R;
  VS[(flip+1)%2] = V1b_R;
  VC[flip] = V2a_R;
  VC[(flip+1)%2] = V2b_R;
  X[flip] = Xa_R;
  X[(flip+1)%2] = Xb_R;


  double NV[2];
  double jet1PT[2];
  double jet2PT[2];


  for(int i = 0; i < 2; i++){

    NV[i] =  VIS_R->GetNElementsInFrame(*VS[i]);
    NV[i] += VIS_R->GetNElementsInFrame(*VC[i]);

    int N = jetID_R.size();
    // std::cout << "In SklimmerAnalysis:  N Jets " << N << std::endl;

    double pTmax[2]; pTmax[0] = -1.; pTmax[1] = -1.;
    for(int j = 0; j < N; j++){
      const RestFrames::RestFrame& frame = VIS_R->GetFrame(jetID_R[j]);

      if(VS[i]->IsSame(frame) || VC[i]->IsSame(frame)){
        double pT = VIS_R->GetLabFrameFourVector(jetID_R[j]).Pt();
        //std::cout << "In SklimmerAnalysis: ijet pT " << pT << std::endl;

        if(pT > pTmax[0]){
          pTmax[1] = pTmax[0];
          pTmax[0] = pT;
        } else {
          if(pT > pTmax[1]) pTmax[1] = pT;
        }
      }
    }

    jet1PT[i] = pTmax[0];
    jet2PT[i] = pTmax[1];


    if(NV[i] > 1){
      RJigsawVariables[Form("RJVars_C_%d_CosTheta",i)     ] = C[i]->GetCosDecayAngle();
      RJigsawVariables[Form("RJVars_P_%d_dPhiGC",i)       ] = G[i]->GetDeltaPhiDecayPlanes(*C[i]);
      RJigsawVariables[Form("RJVars_P_%d_MassRatioGC",i)  ] = (C[i]->GetMass()-X[i]->GetMass())/(G[i]->GetMass()-X[i]->GetMass());
    } else {
      RJigsawVariables[Form("RJVars_C_%d_CosTheta",i)     ] = -10.;
      RJigsawVariables[Form("RJVars_P_%d_dPhiGC",i)       ] = -10.;
      RJigsawVariables[Form("RJVars_P_%d_MassRatioGC",i)  ] = -10.;
    }

    RJigsawVariables[ Form("RJVars_P_%d_CosTheta",i)    ] = G[i]->GetCosDecayAngle();
    RJigsawVariables[ Form("RJVars_P_%d_Jet1_pT",i)     ] = jet1PT[i];
    RJigsawVariables[ Form("RJVars_P_%d_Jet2_pT",i)     ] = jet2PT[i];


    TVector3 P1_G = VS[i]->GetFourVector(*G[i]).Vect();
    TVector3 P2_G = VC[i]->GetFourVector(*G[i]).Vect();


    float Pinv = (P1_G+P2_G).Mag();
    float P1 = P1_G.Mag();
    float P2 = P2_G.Mag();
    RJigsawVariables[ Form("RJVars_P_%d_PInvHS",i) ] = 2*Pinv/(P1+P2+Pinv);

  }


  RJigsawVariables[ "RJVars_V1_N" ] = NV[0];
  RJigsawVariables[ "RJVars_V2_N" ] = NV[1];


  TLorentzVector vV1 = G[0]->GetVisibleFourVector(*G[0]);
  TLorentzVector vV2 = G[1]->GetVisibleFourVector(*G[1]);
  float MG = (vV1.M2()-vV2.M2())/(2.*(vV1.E()-vV2.E()));

  float PG = G[0]->GetMomentum(*GG_R);
  float MGG = 2.*sqrt(PG*PG + MG*MG);
  float gaminvGG = 2.*MG/MGG;
  float gaminv = GG_R->GetVisibleShape();
  float beta = sqrt(1.- gaminv*gaminv);
  float betaGG = sqrt(1.- gaminvGG*gaminvGG);

  //*** velocity difference between 'massive' and 'mass-less'
  float DeltaBetaGG = -(betaGG-beta)/(1.-betaGG*beta);

  //*** delta phi between GG visible decay products and GG decay axis
  float dphiVG = GG_R->GetDeltaPhiDecayVisible();


  RJigsawVariables[ "RJVars_MG"          ] = MG;
  RJigsawVariables[ "RJVars_DeltaBetaGG" ] = DeltaBetaGG;
  RJigsawVariables[ "RJVars_dphiVG"      ] = dphiVG;


  // Signal-like variables end
  ////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  // QCD Variables


  // dphiR and Rptshat (formerly cosPT)
  // for QCD rejection
  double dphiR = GG_R->GetDeltaPhiBoostVisible();
  double PTCM = GG_R->GetFourVector(*LAB_R).Pt();
  double Rptshat = PTCM / (PTCM + GG_R->GetMass()/4.);

  // QCD rejection using the 'background tree'
  // MET 'sibling' in background tree auxillary calculations
  TLorentzVector Psib = I_B->GetSiblingFrame().GetFourVector(*LAB_B);
  TLorentzVector Pmet = I_B->GetFourVector(*LAB_B);
  double Psib_dot_METhat = max(0., Psib.Vect().Dot(MET_TV3.Unit()));
  double Mpar2 = Psib.E()*Psib.E()-Psib.Vect().Dot(MET_TV3.Unit())*Psib.Vect().Dot(MET_TV3.Unit());
  double Msib2 = Psib.M2();
  double MB2 = 2.*(Pmet.E()*Psib.E()-MET_TV3.Dot(Psib.Vect()));
  TVector3 boostPsibM = (Pmet+Psib).BoostVector();


  // QCD rejection variables from 'background tree'
  //double DepthBKG = S_B->GetFrameDepth(I_B);
  //int Nsib = I_B->GetSiblingFrame()->GetNDescendants();
  //double cosBKG = I_B->GetParentFrame()->GetCosDecayAngle();
  //double dphiMsib = fabs(MET_TV3.DeltaPhi(Psib.Vect()));
  double RpsibM = Psib_dot_METhat / (Psib_dot_METhat + MET_TV3.Mag());
  double RmsibM = 1. / ( MB2/(Mpar2-Msib2) + 1.);
  Psib.Boost(-boostPsibM);
  double cosPsibM = -1.*Psib.Vect().Unit().Dot(boostPsibM.Unit());
  cosPsibM = (1.-cosPsibM)/2.;
  double DeltaQCD1 = (cosPsibM-RpsibM)/(cosPsibM+RpsibM);
  double DeltaQCD2 = (cosPsibM-RmsibM)/(cosPsibM+RmsibM);

  RJigsawVariables[ "RJVars_QCD_dPhiR"    ] = dphiR;
  RJigsawVariables[ "RJVars_QCD_Rpt"      ] = Rptshat;
  RJigsawVariables[ "RJVars_QCD_Rmsib"    ] = RmsibM;
  RJigsawVariables[ "RJVars_QCD_Delta2"   ]  = DeltaQCD2;
  RJigsawVariables[ "RJVars_QCD_Rpsib"    ] = RpsibM;
  RJigsawVariables[ "RJVars_QCD_Delta1"   ]  = DeltaQCD1;

  // end
  ////////////////////////////////////////////////////////////////////////////////


  return 0 ;

}

int setupRestFrames(){

  // cleanup previously computed variables
  if ( !LAB_B) delete LAB_B;
  if ( !S_B) delete S_B;
  if ( !V_B) delete V_B;
  if ( !I_B) delete I_B;
  if ( !INV_B) delete INV_B;
  if ( !VIS_B) delete VIS_B;
  if ( !MinMass_B) delete MinMass_B;
  if ( !Rapidity_B) delete Rapidity_B;
  if ( !LAB_R) delete LAB_R;
  if ( !GG_R) delete GG_R;
  if ( !Ga_R) delete Ga_R;
  if ( !Gb_R) delete Gb_R;
  if ( !Ca_R) delete Ca_R;
  if ( !Cb_R) delete Cb_R;
  if ( !V1a_R) delete V1a_R;
  if ( !V2a_R) delete V2a_R;
  if ( !Xa_R) delete Xa_R;
  if ( !V1b_R) delete V1b_R;
  if ( !V2b_R) delete V2b_R;
  if ( !Xb_R) delete Xb_R;
  if ( !INV_R) delete INV_R;
  if ( !VIS_R) delete VIS_R;
  if ( !MinMassJigsaw_R) delete MinMassJigsaw_R;
  if ( !RapidityJigsaw_R) delete RapidityJigsaw_R;
  if ( !ContraBoostJigsaw_R) delete ContraBoostJigsaw_R;
  if ( !HemiJigsaw_R) delete HemiJigsaw_R;
  if ( !CaHemiJigsaw_R) delete CaHemiJigsaw_R;
  if ( !CbHemiJigsaw_R) delete CbHemiJigsaw_R;

  LAB_B = new RestFrames::LabRecoFrame("LAB_B","LAB_B");

  S_B = new RestFrames::SelfAssemblingRecoFrame("CM_B","CM_B");
  V_B = new RestFrames::VisibleRecoFrame("V_B","Vis_B");
  I_B = new RestFrames::InvisibleRecoFrame("I_B","Iinv_B");
  INV_B = new RestFrames::InvisibleGroup ("INV_B","Invisible State Jigsaws");
  VIS_B = new RestFrames::CombinatoricGroup("VIS_B","Visible Object Jigsaws");

  MinMass_B = new RestFrames::SetMassInvJigsaw("MINMASS_JIGSAW_B", "Invisible system mass Jigsaw");
  Rapidity_B = new RestFrames::SetRapidityInvJigsaw("RAPIDITY_JIGSAW_B", "Invisible system rapidity Jigsaw");

  LAB_R = new RestFrames::LabRecoFrame("LAB_R","LAB");
  GG_R = new RestFrames::DecayRecoFrame("GG_R","#tilde{g}#tilde{g}");
  Ga_R = new RestFrames::DecayRecoFrame("Ga_R","#tilde{g}_{a}");
  Gb_R = new RestFrames::DecayRecoFrame("Gb_R","#tilde{g}_{b}");
  Ca_R = new RestFrames::DecayRecoFrame("Ca_R","C_{a}");
  Cb_R = new RestFrames::DecayRecoFrame("Cb_R","C_{b}");
  V1a_R = new RestFrames::VisibleRecoFrame("V1a_R","j_{1a}");
  V2a_R = new RestFrames::VisibleRecoFrame("V2a_R","j_{2a}");
  Xa_R = new RestFrames::InvisibleRecoFrame("Xa_R","#tilde{#chi}_{a}");
  V1b_R = new RestFrames::VisibleRecoFrame("V1b_R","j_{1b}");
  V2b_R = new RestFrames::VisibleRecoFrame("V2b_R","j_{2b}");
  Xb_R = new RestFrames::InvisibleRecoFrame("Xb_R","#tilde{#chi}_{b}");
  INV_R = new RestFrames::InvisibleGroup ("INV_R","WIMP Jigsaws");
  VIS_R = new RestFrames::CombinatoricGroup("VIS","Visible Object Jigsaws");
  MinMassJigsaw_R = new RestFrames::SetMassInvJigsaw("MINMASS_R", "Invisible system mass Jigsaw");
  RapidityJigsaw_R = new RestFrames::SetRapidityInvJigsaw("RAPIDITY_R", "Invisible system rapidity Jigsaw");
  ContraBoostJigsaw_R = new RestFrames::ContraBoostInvJigsaw("CONTRA_R","Contraboost invariant Jigsaw");
  HemiJigsaw_R = new RestFrames::MinMassesCombJigsaw ("HEM_JIGSAW_R","Minimize m _{V_{a,b}} Jigsaw");
  CaHemiJigsaw_R = new RestFrames::MinMassesCombJigsaw("CbHEM_JIGSAW_R","Minimize m _{C_{a}} Jigsaw");
  CbHemiJigsaw_R = new RestFrames::MinMassesCombJigsaw("CaHEM_JIGSAW_R","Minimize m _{C_{b}} Jigsaw");

  INV_B->AddFrame(*I_B);
  VIS_B->AddFrame(*V_B);
  VIS_B->SetNElementsForFrame(*V_B,1,false);

  LAB_B->SetChildFrame(*S_B);
  S_B->AddChildFrame(*V_B);
  S_B->AddChildFrame(*I_B);

  LAB_B->InitializeTree();

// Will just set invisible mass to zero
  INV_B->AddJigsaw(*MinMass_B);

// will set rapidity to zero
  INV_B->AddJigsaw(*Rapidity_B);
  Rapidity_B->AddVisibleFrames( (LAB_B->GetListVisibleFrames()) );

  LAB_B->InitializeAnalysis();

  //
  //
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // SPARTICLE TREE //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //
  //

  // Set up 'signal-like' analysis tree
  LAB_R->SetChildFrame(*GG_R);
  GG_R->AddChildFrame(*Ga_R);
  GG_R->AddChildFrame(*Gb_R);
  Ga_R->AddChildFrame(*V1a_R);
  Ga_R->AddChildFrame(*Ca_R);
  Ca_R->AddChildFrame(*V2a_R);
  Ca_R->AddChildFrame(*Xa_R);
  Gb_R->AddChildFrame(*V1b_R);
  Gb_R->AddChildFrame(*Cb_R);
  Cb_R->AddChildFrame(*V2b_R);
  Cb_R->AddChildFrame(*Xb_R);


  //if(!LAB_R->InitializeTree()) cout << "Problem with signal-like reconstruction tree" << endl;
  LAB_R->InitializeTree();

  INV_R->AddFrame(*Xa_R);
  INV_R->AddFrame(*Xb_R);
  // visible frames in first decay step must always have at least one element
  VIS_R->AddFrame(*V1a_R);
  VIS_R->AddFrame(*V1b_R);
  VIS_R->SetNElementsForFrame(*V1a_R,1,false);
  VIS_R->SetNElementsForFrame(*V1b_R,1,false);
  // visible frames in second decay step can have zero elements
  VIS_R->AddFrame(*V2a_R);
  VIS_R->AddFrame(*V2b_R);
  VIS_R->SetNElementsForFrame(*V2a_R,0,false);
  VIS_R->SetNElementsForFrame(*V2b_R,0,false);

  INV_R->AddJigsaw(*MinMassJigsaw_R);
  INV_R->AddJigsaw(*RapidityJigsaw_R);
  RapidityJigsaw_R->AddVisibleFrames((LAB_R->GetListVisibleFrames()));
  INV_R->AddJigsaw(*ContraBoostJigsaw_R);
  ContraBoostJigsaw_R->AddVisibleFrames((Ga_R->GetListVisibleFrames()), 0);
  ContraBoostJigsaw_R->AddVisibleFrames((Gb_R->GetListVisibleFrames()), 1);
  ContraBoostJigsaw_R->AddInvisibleFrames((Ga_R->GetListInvisibleFrames()), 0);
  ContraBoostJigsaw_R->AddInvisibleFrames((Gb_R->GetListInvisibleFrames()), 1);
  VIS_R->AddJigsaw(*HemiJigsaw_R);
  HemiJigsaw_R->AddFrame(*V1a_R,0);
  HemiJigsaw_R->AddFrame(*V1b_R,1);
  HemiJigsaw_R->AddFrame(*V2a_R,0);
  HemiJigsaw_R->AddFrame(*V2b_R,1);
  VIS_R->AddJigsaw(*CaHemiJigsaw_R);
  CaHemiJigsaw_R->AddFrame(*V1a_R,0);
  CaHemiJigsaw_R->AddFrame(*V2a_R,1);
  // CaHemiJigsaw_R->AddFrame(*Xa_R,1); //This and the next line removed because they allow for inv particle to be the only thing left in second order frame.
  VIS_R->AddJigsaw(*CbHemiJigsaw_R);
  CbHemiJigsaw_R->AddFrame(*V1b_R,0);
  CbHemiJigsaw_R->AddFrame(*V2b_R,1);
  // CbHemiJigsaw_R->AddFrame(*Xb_R,1);

  //if(!LAB_R->InitializeAnalysis()) cout << "Problem with signal-tree jigsaws" << endl;
  LAB_R->InitializeAnalysis();

  return 0;
}
