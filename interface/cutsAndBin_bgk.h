#ifndef CutAndBinCollection_C
#define CutAndBinCollection_C

#include <TF1.h>
#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>

struct HydjetEvent{

	Int_t event;
	Float_t b;
	Float_t npart;
	Float_t ncoll;
	Float_t nhard;
	Float_t phi0;
	Float_t scale;

	Int_t n[3];
	Float_t ptav[3];

	Int_t mult;
	std::vector<Float_t> E;
	std::vector<Float_t> mass;
	std::vector<Float_t> pt;
	std::vector<Float_t> eta;
	std::vector<Float_t> phi;
	std::vector<Int_t> pdg;
	std::vector<Int_t> chg;
	std::vector<Int_t> sube;
	std::vector<Int_t> sta;
	std::vector<Int_t> matchingID;
	std::vector<Int_t> nMothers;
	std::vector<std::vector<Int_t> > motherIndex;
	std::vector<Int_t> nDaughters;
	std::vector<std::vector<Int_t> > daughterIndex;

	Float_t vx;
	Float_t vy;
	Float_t vz;
	Float_t vr;

};

struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlus, KaPlus; };
ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.139570, 0.49367 };


int kPPDATA = 0 ;
int kPADATA = 1 ;
int kAADATA = 2 ; // L1 doubleMu 0
int kPPMC = 3 ;
int kPAMC = 4 ;
int kAAMC = 5 ;
int kAADATAPeri = 6 ;
int kAADATACentL3 = 7 ;

TString getCollID( int collid ) {
	if ( collid == kPPDATA ) return "PP_DATA";
	else if ( collid == kPADATA ) return "PA_DATA";
	else if ( collid == kAADATA ) return "AA_DATA";
	else if ( collid == kPPMC ) return "PP_MC";
	else if ( collid == kPAMC ) return "PA_MC";
	else if ( collid == kAAMC ) return "AA_MC";
	else if ( collid == kAADATAPeri ) return "AA_DATA_PeriL1";
	else if ( collid == kAADATACentL3 ) return "AA_DATA_CentL3";
	else return "none";
}

int kEPl2HF = 0;
int kEPOppositeHF = 1;
int kEPSameSideHF = 2;


TString getEPSel( int eventPln) {
	if ( eventPln == kEPl2HF)  return "BothHFs";
	else if ( eventPln == kEPOppositeHF ) return "OppositeHF" ;
	else if ( eventPln == kEPSameSideHF ) return "SameSideHF" ;
	else return "none";
}

class diJet { 
	public:
		diJet() :
			run(0), lumi(0), event(0), hfsum(0), hfplus(0), hfminus(0),  // 4
			vz(-99),  mass(-1), pt(-1), y(999), phi(999), eta(999), dphi(-99), // 7
			pt1(-1), eta1(-1), phi1(-1), pu1(-1),         // 3
			pt2(-1), eta2(-1), phi2(-1), pu2(-1)        // 3 
	{}

		int run;
		int lumi;
		int event;
		float hfsum;
		float hfplus;
		float hfminus;
		float vz;
		float mass;
		float pt;
		float y;
		float phi;    
		float eta;
		float dphi;
		float pt1; 
		float eta1;
		float phi1;
		float pu1; 
		float pt2;
		float eta2;
		float phi2;    
		float pu2; 

		void clear() {
			run = -99; lumi = -99; event = -99; hfsum = -99; hfplus=-99; hfminus=-99;  // 4
			vz = -99;  mass =-99; pt = -99; y = -99 ; phi = -99; eta = -99; dphi=-99; // 7
			pt1 = -99; eta1 = -99; phi1 = -99; pu1=-99;      // 4
			pt2 = -99; eta2 = -99; phi2 = -99; pu2=-99;      // 4
		}

};

class UPCdiJet {
	public:
		UPCdiJet() :
			nJet(-1), mass(-1), pt(-1), y(999), phi(999), eta(999), dphi(-99), dpt(-99), deta(-99), aj(-99),// 7
			pt1(-99), eta1(-99), phi1(-99), e1(-1),       // 3
			pt2(-99), eta2(-99), phi2(-99), e2(-1),         // 3
			PT(-99), QT(-99), aniPhi(-99), cos2Phi(-99),
			aniPhi_c6(-99), cos2Phi_c6(-99),
                        aniPhi_c7(-99), cos2Phi_c7(-99),
                        aniPhi_c8(-99), cos2Phi_c8(-99),
			aniPhi_c10(-99), cos2Phi_c10(-99)
	{}

		int nJet;
		float mass;
		float pt;
		float y;
		float phi;
		float eta;
		float dphi;
		float dpt;
		float deta;
		float aj;
		float pt1;
		float eta1;
		float phi1;
		float e1;
		float pt2;
		float eta2;
		float phi2;
		float e2;
		float PT;
		float QT;
		float aniPhi;
		float cos2Phi;
                float aniPhi_c6;
                float cos2Phi_c6;
                float aniPhi_c7;
                float cos2Phi_c7;
                float aniPhi_c8;
                float cos2Phi_c8;
                float aniPhi_c10;
                float cos2Phi_c10;

		void clear() {
			nJet = -1; mass = -99; pt = -99; y = -99; phi = -99; eta = -99; dphi=-99; dpt=-99; deta=-99; aj=-99;// 7
			pt1 = -99; eta1 = -99; phi1 = -99; e1 = -99;    // 4
			pt2 = -99; eta2 = -99; phi2 = -99; e2 = -99;    // 4
			PT = -99; QT = -99; aniPhi = -99; cos2Phi = -99;
			aniPhi_c6 = -99; cos2Phi_c6 = -99;
			aniPhi_c7 = -99; cos2Phi_c7 = -99;
			aniPhi_c8 = -99; cos2Phi_c8 = -99;
                        aniPhi_c10= -99; cos2Phi_c10= -99;
		}

};
TString djBranchString = "nJet/I:mass/F:pt:y:phi:eta:dphi:dpt:deta:aj:pt1:eta1:phi1:e1:pt2:eta2:phi2:e2:PT:QT:aniPhi:cos2Phi:aniPhi_c6:cos2Phi_c6:aniPhi_c7:cos2Phi_c7:aniPhi_c8:cos2Phi_c8:aniPhi_c10:cos2Phi_c10";

class UPCnTrk {
	public:
		UPCnTrk() :
			nTrack(-99),
			nTrkabsEtaover1p5(-99), nTrkabsEtaunder1p5(-99),
			nTrketa0to0p5(-99), nTrketa0p5to1(-99), nTrketa1to1p5(-99), nTrketa1p5to2(-99), nTrketa2to2p5(-99),
			nTrketam0to0p5(-99), nTrketam0p5to1(-99), nTrketam1to1p5(-99), nTrketam1p5to2(-99), nTrketam2to2p5(-99),
			TrkPhi(-99), TrkEta(-99)
	{}

		int nTrack;
		int nTrkabsEtaover1p5;
		int nTrkabsEtaunder1p5;
		int nTrketa0to0p5;
		int nTrketa0p5to1;
		int nTrketa1to1p5;
		int nTrketa1p5to2;
		int nTrketa2to2p5;
		int nTrketam0to0p5;
		int nTrketam0p5to1;
		int nTrketam1to1p5;
		int nTrketam1p5to2;
		int nTrketam2to2p5;
		float TrkPhi;
		float TrkEta;

		void clear() {
			nTrack = 0;
			nTrkabsEtaover1p5 = 0; nTrkabsEtaunder1p5 = 0;
			nTrketa0to0p5 = 0; nTrketa0p5to1 = 0; nTrketa1to1p5 = 0; nTrketa1p5to2 = 0; nTrketa2to2p5 = 0;
			nTrketam0to0p5 = 0; nTrketam0p5to1 = 0; nTrketam1to1p5 = 0; nTrketam1p5to2 = 0; nTrketam2to2p5 = 0;
			TrkPhi = -99; TrkEta = -99;
		}

};
TString nTrkString = "nTrk/I:nTrkabsEtaover1p5:nTrkabsEtaunder1p5:nTrketa0to0p5:nTrketa0p5to1:nTrketa1to1p5:nTrketa1p5to2:nTrketa2to2p5:nTrketam0to0p5:nTrketam0p5to1:nTrketam1to1p5:nTrketam1p5to2:nTrketam2to2p5/F:TrkPhi:TrkEta"; 

class UPCEvent { 
	public:
		UPCEvent() :
			run(0), lumi(0), event(0) , nPho(-99), nTrk(-99) , hfsum(0),  hfplus(0), hfminus(0), vz(-99)
	{}

		int run;
		int lumi;
		int event;
		int nPho;
		float nTrk;
		float hfsum;
		float hfplus;
		float hfminus;
		float vz;

		void clear() {
			run = -99; lumi = -99; event = -99; nPho=-99; nTrk=-99; hfsum = -99; vz = -99; hfplus=-99; hfminus=-99; 
		}

};
TString eventBranchString = "run/I:lumi:event:nPho:nTrk:hfsum/F:hfplus:hfminus:vz";

class DiPhoton { 
	public:
		DiPhoton() :
			mass(-1), pt(-1), y(999), phi(999), eta(999), dphi(-99),
			pt1(-1), eta1(-1), phi1(-1), sii1(-1), hoe1(-1), 
			pt2(-1), eta2(-1), phi2(-1), sii2(-1), hoe2(-1),
			eiso1(-1), hiso1(-1), tiso1(-1),
			eiso2(-1), hiso2(-1), tiso2(-1)

	{}

		float mass;
		float pt;
		float y;
		float phi;    
		float eta;
		float dphi;

		float pt1; 
		float eta1;
		float phi1;
		float sii1;
		float hoe1;

		float pt2;
		float eta2;
		float phi2;    
		float sii2;
		float hoe2;

		float eiso1; 
		float hiso1; 
		float tiso1; 
		float eiso2; 
		float hiso2; 
		float tiso2; 


		void clear() {
			mass =-99; pt = -99; y = -99 ; phi = -99; eta = -99; dphi=-99;
			pt1 = -99; eta1 = -99; phi1 = -99; sii1= 99; hoe1=99;
			pt2 = -99; eta2 = -99; phi2 = -99; sii2= 99; hoe2=99;
			eiso1 = -99; hiso1 = -99; tiso1=-99;
			eiso2 = -99; hiso2 = -99; tiso2=-99;
		}

};
TString dpBranchString = "mass/F:pt:y:phi:eta:dphi:pt1:eta1:phi1:sii1:hoe1:pt2:eta2:phi2:sii2:hoe2:eiso1:hiso1:tiso1:eiso2:hiso2:tiso2";



class DiPhoPF { 
	public:
		DiPhoPF() :
			mass(-1), pt(-1), y(999), phi(999), eta(999), dphi(-99),
			pt1(-1), eta1(-1), phi1(-1),
			pt2(-1), eta2(-1), phi2(-1)

	{}

		float mass;
		float pt;
		float y;
		float phi;    
		float eta;
		float dphi;

		float pt1; 
		float eta1;
		float phi1;

		float pt2;
		float eta2;
		float phi2;    


		void clear() {
			mass =-99; pt = -99; y = -99 ; phi = -99; eta = -99; dphi=-99;
			pt1 = -99; eta1 = -99; phi1 = -99; 
			pt2 = -99; eta2 = -99; phi2 = -99; 
		}

};
TString dpBranchStringPF = "mass/F:pt:y:phi:eta:dphi:pt1:eta1:phi1:pt2:eta2:phi2";


#endif
