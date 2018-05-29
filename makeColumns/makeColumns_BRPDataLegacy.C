#include "../interface/cutsAndBin_bgk.h"
#include "../interface/commonUtility.h"
#include "TProfile.h"
#include "TProfile2D.h"

void makeColumns(){
	//This script is for plotting the V2 vs PT. The event mixing effect is subtracted. 
	//Both BRP and FRP are summed up. GEN-level RAPGAP is also histogramed.

	Int_t forMixingQuota = 5;
	TString treeFolder = "/home/samboren/Workspace/upcDijAnaCode/data/skimTree/";

	TFile *fMCTrig = new TFile(treeFolder+"TrkCutsppreco_TrueMC_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0.root");
	TFile *fDATA = new TFile(treeFolder+"TrkCutsppreco_upcDiJetSkim2ndVer171107_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0_2017y_11m_7d_21h_16m.root");
	TFile *fGEN = new TFile(treeFolder+"TrkCutsppreco_GenMC_trig_jetCollectionak4PFJetAnalyzer_minJetPt0.root");
	//"TrkCutsppreco_GenMC_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0.root");
	//"TrkCutsppreco_GenMC_trig_jetCollectionak4PFJetAnalyzer_minJetPt0.root");

	TTree *dijetTree_DATA   = (TTree*)fDATA->Get("dijet");
	TTree *trkTree_DATA     = (TTree*)fDATA->Get("fullTrkTree");
	TTree *calTree_DATA = (TTree*)fDATA->Get("Cal");
	TTree *evtTree_DATA = (TTree*)fDATA->Get("evt");
	dijetTree_DATA->AddFriend(trkTree_DATA);
	dijetTree_DATA->AddFriend(calTree_DATA);
	dijetTree_DATA->AddFriend(evtTree_DATA);

	Float_t BRP;
	Float_t FRP;

	UPCdiJet djObj;
	Float_t HFplusmax;
	Float_t HFminusmax;
	TLorentzVector jet1a;
	TLorentzVector jet2a;

	TLorentzVector jet1b;
	TLorentzVector jet2b;

	TLorentzVector jetplus;
	TLorentzVector jetminus;

	TLorentzVector jetplusGen;
	TLorentzVector jetminusGen;

	Float_t zee = 0.5;
	Float_t anglePerp;
	Float_t anglePerpGen;

	Float_t d_phi;


	dijetTree_DATA->SetBranchAddress("dj",&djObj);
	calTree_DATA->SetBranchAddress("HFplusmax",&HFplusmax);
	calTree_DATA->SetBranchAddress("HFminusmax",&HFminusmax);

	Float_t trkPt[1000];
	Float_t trkEta[1000];
	Int_t nTrks;
	Float_t etaEdg = 2.5;
	Float_t ptThreshold = .400;

	Int_t nEvents;
	Int_t matched = 0;
	Int_t tooLoop = dijetTree_DATA->GetEntries();

	dijetTree_DATA->SetBranchAddress("dj",&djObj);
	calTree_DATA->SetBranchAddress("HFplusmax",&HFplusmax);
	calTree_DATA->SetBranchAddress("HFminusmax",&HFminusmax);
	trkTree_DATA->SetBranchAddress("pT",&trkPt);
	trkTree_DATA->SetBranchAddress("Eta",&trkEta);
	trkTree_DATA->SetBranchAddress("ntrk",&nTrks);

	UPCEvent event;
	evtTree_DATA->SetBranchAddress("event",&event);	

	nEvents = dijetTree_DATA->GetEntries();
	Int_t mixingNum = forMixingQuota;
	Float_t vz1 = 0;
	Float_t vz2 = 0;

	Bool_t leftSide = false;
	Bool_t rightSide = false;

	TLorentzVector v01, v02, vv1, vv2, vqt, vpt, m, w, g, h, p, q;
	TVector2 n;

	//Define variables
	Double_t pi = TMath::Pi();
	float v1_norm, v2_norm, n_norm, m_norm;
	float c12, s12, v1v2, a12;
	Float_t pt1, pt2, eta1, eta2, phi1, phi2, e1, e2;
	Float_t px1, py1, pz1, ee1, px2, py2, pz2, ee2;
	float ptvqt, ptvpt;
	float sign=0.0;

	Float_t vtx1,vtx2,frp1,brp1,frp2,brp2;

	gROOT->Reset();

	TFile *fmd = new TFile("columns.root","RECREATE");
	TNtuple *mixedDataNtuple = new TNtuple("mixedDataNtuple","mixed data from tree","e1:px1:py1:pz1:e2:px2:py2:pz2:frp1:brp1:frp2:brp2:vtx1:vtx2");

	//FILE *fp = fopen("ntuple.dat","r"); 
/*
	//Mixed data, nominal cuts
	nEvents=100;
	for(Long64_t i=0; i<nEvents; i++){
		cout<<"mixed data lrg: "<<1.0*i/nEvents<<endl;
		leftSide = false;
		rightSide = false;
		dijetTree_DATA->GetEntry(i);

		if(djObj.e1>djObj.e2){
			jet1a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet2a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);
		}else{
			jet2a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet1a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);
		}

		if(djObj.nJet!=2 || abs(djObj.eta1)>1.8 || djObj.pt1<20 || abs(djObj.eta2)>1.8 || djObj.pt2<15 || djObj.mass<35 || djObj.dphi<2) continue;
		jetplus=jet1a+jet2a;
		jetminus=(1-zee)*jet1a - zee*jet2a;
		BRP =  rapidityGapAlgorithm("left",etaEdg,trkPt,trkEta,ptThreshold,nTrks);
		FRP =  rapidityGapAlgorithm("right",etaEdg,trkPt,trkEta,ptThreshold,nTrks);		                           
		vz1 = event.vz;

		if( !(BRP>1.2 && BRP>FRP) && !(FRP>1.2 && FRP>BRP)  ) continue;
		if( FRP>1.2 && FRP>BRP) rightSide = true;
		if( BRP>1.2 && BRP>FRP) leftSide = true;

		frp1=FRP;
		brp1=BRP;

		for(Long64_t j=i+1; j<nEvents-20; j++){
			if(mixingNum==0) continue;
			dijetTree_DATA->GetEntry(j);
			if(i==j || djObj.nJet!=2 || abs(djObj.eta1)>1.8 || djObj.pt1<20 || abs(djObj.eta2)>1.8 || djObj.pt2<15 || djObj.mass<35 || djObj.dphi<2) continue;

			if(djObj.e1>djObj.e2){
				jet1b.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
				jet2b.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);
			}else{
				jet2b.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
				jet1b.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);
			}

			BRP =  rapidityGapAlgorithm("left",etaEdg,trkPt,trkEta,ptThreshold,nTrks);
			FRP =  rapidityGapAlgorithm("right",etaEdg,trkPt,trkEta,ptThreshold,nTrks);
			vz2 = event.vz;
			if( !(BRP>1.2 && BRP>FRP && leftSide) && !(FRP>1.2 && FRP>BRP && rightSide)  ) continue;
			if( abs(vz1-vz2)>3 ) continue;

			frp2=FRP;
			brp2=BRP;

			vtx1 = vz1;
			vtx2 = vz2;

			mixingNum--;
			jetplus=jet1a+jet2b;
			jetminus=(1-zee)*jet1a - zee*jet2b;

			TVector2 p, v1;
			TVector2 q, v2;
			p.Set(jet1a[0],jet1a[1]);
			q.Set(jet2b[0],jet2b[1]);
			//Computing Qt and Pt as 2-Vectors
			v1.Set(p.X() + q.X(),p.Y()+q.Y());
			v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
			//computing the norm of Qt and Pt
			v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
			v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
			//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
			TVector2 v1unit, v2unit;
			v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
			v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
			//Computing the dot product of Qt-hat and Pt-hat
			v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
			//The dot product is the cosine of the angle
			c12 = v1v2  ;
			//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
			n.Set(v1unit.Y(),-v1unit.X()) ;
			n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
			//Sine of the angle
			s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
			//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
			a12 = atan2(s12, c12);			
			if (a12>=0) a12 = a12;
			if (a12<0) a12 = a12 + 2*pi;
			//Computing the cos(2phi) using trigonometry expression
			c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
			anglePerp=a12;

			mixedDataNtuple->Fill(jet1a.E(),jet1a.Px(),jet1a.Py(),jet1a.Pz(),jet2b.E(),jet2b.Px(),jet2b.Py(),jet2b.Pz(),frp1,brp1,frp2,brp2,vtx1,vtx2);

		}
		mixingNum=forMixingQuota;
	}
*/

	//The pure data, nominal cuts
	nEvents=dijetTree_DATA->GetEntries();

	Int_t superCount = 0;

	TNtuple *dataNtuple = new TNtuple("dataNtuple","data from tree","e1:px1:py1:pz1:e2:px2:py2:pz2:frp1:brp1:vtx1");

	for(Long64_t i=0; i<nEvents; i++){
		matched=0;
//		cout<<"data lrg "<<1.0*i/nEvents<<endl;
		dijetTree_DATA->GetEntry(i);

		if(djObj.e1>djObj.e2){
			jet1a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet2a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);
		}else{
			jet2a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet1a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);
		}

		if(djObj.nJet!=2 || abs(djObj.eta1)>1.8 || djObj.pt1<20 || abs(djObj.eta2)>1.8 || djObj.pt2<15 || djObj.mass<35 || djObj.dphi<2) continue;
		BRP =  rapidityGapAlgorithm("left",etaEdg,trkPt,trkEta,ptThreshold,nTrks);
		FRP =  rapidityGapAlgorithm("right",etaEdg,trkPt,trkEta,ptThreshold,nTrks);
		jetplus=jet1a+jet2a;
		jetminus=(1-zee)*jet1a - zee*jet2a;

		if( (BRP<1.2 || BRP<FRP) /* !(FRP>1.2 && FRP>BRP)*/  ) continue;
		vz1 = event.vz;

		jetplus=jet1a+jet2a;
		jetminus=(1-zee)*jet1a - zee*jet2a;

		TVector2 p, v1;
		TVector2 q, v2;
		p.Set(jet1a[0],jet1a[1]);
		q.Set(jet2a[0],jet2a[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		TVector2 v1unit, v2unit;
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);			
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		anglePerp=a12;
		
		superCount++;
		cout<<superCount<<endl;
		dataNtuple->Fill(jet1a.E(),jet1a.Px(),jet1a.Py(),jet1a.Pz(),jet2a.E(),jet2a.Px(),jet2a.Py(),jet2a.Pz(),FRP,BRP,vz1);
	}

/*
	TTree *dijetTree_MC  	= (TTree*)fGEN->Get("dijet");
	TTree *trkTree_MC     	= (TTree*)fGEN->Get("fullTrkTree");
	TTree *calTree_MC 		= (TTree*)fGEN->Get("Cal");
	TTree *evtTree_MC 		= (TTree*)fGEN->Get("evt");
	TTree *dijetTree_GEN 	= (TTree*)fGEN->Get("djGen");
	TTree *genTree      	= (TTree*)fGEN->Get("hi");
	dijetTree_GEN->AddFriend(genTree);
	genTree->SetMakeClass(1);
	dijetTree_MC->AddFriend(trkTree_MC);
	dijetTree_MC->AddFriend(calTree_MC);
	dijetTree_MC->AddFriend(evtTree_MC);
	dijetTree_MC->AddFriend(dijetTree_GEN);
	vector<float>   *etaGen = new vector<float>;
	vector<float>   *ptGen  = new vector<float>;
	UPCdiJet djObjGen;
	dijetTree_MC->SetBranchAddress("dj",&djObj);
	dijetTree_GEN->SetBranchAddress("djGen",&djObjGen);
	calTree_MC->SetBranchAddress("HFplusmax",&HFplusmax);
	calTree_MC->SetBranchAddress("HFminusmax",&HFminusmax);
	trkTree_MC->SetBranchAddress("pT",&trkPt);
	trkTree_MC->SetBranchAddress("Eta",&trkEta);
	trkTree_MC->SetBranchAddress("ntrk",&nTrks);
	genTree->SetBranchAddress("pt",&ptGen);
	genTree->SetBranchAddress("eta",&etaGen);
	evtTree_MC->SetBranchAddress("event",&event);
	mixingNum = forMixingQuota;
	//nEvents = dijetTree_MC->GetEntries();

	nEvents = 100;
	TNtuple *mixedMCNtuple = new TNtuple("mixedMCNtuple","mixed mc from tree","e1:px1:py1:pz1:e2:px2:py2:pz2:frp1:brp1:frp2:brp2:vtx1:vtx2");
	for(Long64_t i=0; i<nEvents; i++){
		dijetTree_GEN->GetEntry(i);
		cout<<"mixed mc:"<<1.0*i/nEvents<<endl;

		if(djObjGen.pt1>djObjGen.pt2){
			jet1a.SetPtEtaPhiE(djObjGen.pt1,djObjGen.eta1,djObjGen.phi1,djObjGen.e1);
			jet2a.SetPtEtaPhiE(djObjGen.pt2,djObjGen.eta2,djObjGen.phi2,djObjGen.e2);

			if( abs(djObjGen.eta1)>1.8 || djObjGen.pt1<20 || abs(djObjGen.eta2)>1.8 || djObjGen.pt2<15 || djObjGen.mass<35 || djObjGen.dphi<2) continue;
			BRP =  rapidityGapAlgorithmGEN("left",etaEdg,ptGen,etaGen,ptThreshold);
			FRP =  rapidityGapAlgorithmGEN("right",etaEdg,ptGen,etaGen,ptThreshold);

			if( BRP<1.2 ) continue;
			if( BRP<FRP ) continue;

			frp1=FRP;
			brp1=BRP;

			vz1 = event.vz;
			for(Long64_t j=i+1; j<nEvents; j++){
				if(mixingNum==0) continue;
				dijetTree_GEN->GetEntry(j);
				if( abs(djObjGen.eta1)>1.8 || djObjGen.pt1<20 || abs(djObjGen.eta2)>1.8 || djObjGen.pt2<15 || djObjGen.mass<35 || djObjGen.dphi<2) continue;
				BRP =  rapidityGapAlgorithmGEN("left",etaEdg,ptGen,etaGen,ptThreshold);
				FRP =  rapidityGapAlgorithmGEN("right",etaEdg,ptGen,etaGen,ptThreshold);
				if( BRP<1.2 ) continue;
				if( BRP<FRP ) continue;

				frp2=FRP;
				brp2=BRP;

				vz2 = event.vz;
				if( abs(vz1-vz2)>3 ) continue;
				mixingNum--;

				dijetTree_GEN->GetEntry(i);
				jet1a.SetPtEtaPhiE(djObjGen.pt1,djObjGen.eta1,djObjGen.phi1,djObjGen.e1);

				dijetTree_GEN->GetEntry(j);
				jet2b.SetPtEtaPhiE(djObjGen.pt2,djObjGen.eta2,djObjGen.phi2,djObjGen.e2);

				if(jet1a.E()<jet2b.E()){
					TLorentzVector v;
					v = jet2b;
					jet2b = jet1a;
					jet1a = v;
				}

				jetplusGen=jet1a+jet2b;
				jetminusGen=(1-zee)*jet1a - zee*jet2b;

				TVector2 p, v1;
				TVector2 q, v2;
				p.Set(jet1a[0],jet1a[1]);
				q.Set(jet2b[0],jet2b[1]);
				//Computing Qt and Pt as 2-Vectors
				v1.Set(p.X() + q.X(),p.Y()+q.Y());
				v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
				//computing the norm of Qt and Pt
				v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
				v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
				//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
				TVector2 v1unit, v2unit;
				v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
				v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
				//Computing the dot product of Qt-hat and Pt-hat
				v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
				//The dot product is the cosine of the angle
				c12 = v1v2  ;
				//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
				n.Set(v1unit.Y(),-v1unit.X()) ;
				n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
				//Sine of the angle
				s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
				//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
				a12 = atan2(s12, c12);			
				if (a12>=0) a12 = a12;
				if (a12<0) a12 = a12 + 2*pi;
				//Computing the cos(2phi) using trigonometry expression
				c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
				anglePerpGen=a12;

				mixedMCNtuple->Fill(jet1a.E(),jet1a.Px(),jet1a.Py(),jet1a.Pz(),jet2b.E(),jet2b.Px(),jet2b.Py(),jet2b.Pz(),frp1,brp1,frp2,brp2,vz1,vz2);
			}
			mixingNum=forMixingQuota;
		}
	}


        TNtuple *MCNtuple = new TNtuple("MCNtuple","mc from tree","e1:px1:py1:pz1:e2:px2:py2:pz2:frp1:brp1:vtx1");

	for(Long64_t i=0; i<nEvents; i++){
		dijetTree_GEN->GetEntry(i);
		cout<<"pure mc:"<<1.0*i/nEvents<<endl;
		if( abs(djObjGen.eta1)>1.8 || djObjGen.pt1<20 || abs(djObjGen.eta2)>1.8 || djObjGen.pt2<15 || djObjGen.mass<35 || djObjGen.dphi<2) continue;
		BRP =  rapidityGapAlgorithmGEN("left",etaEdg,ptGen,etaGen,ptThreshold);
		FRP =  rapidityGapAlgorithmGEN("right",etaEdg,ptGen,etaGen,ptThreshold);
		if( !(BRP>1.2) ) continue;
		if( !(BRP>FRP) ) continue;


		dijetTree_GEN->GetEntry(i);
		jet1a.SetPtEtaPhiE(djObjGen.pt1,djObjGen.eta1,djObjGen.phi1,djObjGen.e1);
		jet2a.SetPtEtaPhiE(djObjGen.pt2,djObjGen.eta2,djObjGen.phi2,djObjGen.e2);

                vz1 = event.vz;

		if(jet1a.E()<jet2a.E()){
			TLorentzVector v;
			v = jet2a;
			jet2a = jet1a;
			jet1a = v;
		}

		jetplusGen=jet1a+jet2a;
		jetminusGen=(1-zee)*jet1a - zee*jet2a;

		TVector2 p, v1;
		TVector2 q, v2;
		p.Set(jet1a[0],jet1a[1]);
		q.Set(jet2a[0],jet2a[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		TVector2 v1unit, v2unit;
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);			
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		anglePerpGen=a12;

                MCNtuple->Fill(jet1a.E(),jet1a.Px(),jet1a.Py(),jet1a.Pz(),jet2a.E(),jet2a.Px(),jet2a.Py(),jet2a.Pz(),FRP,BRP,vz1);


	}
*/
	fmd->Write();
}
