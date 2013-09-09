#ifndef JET_SMEAR_STUDY_HH
#define JET_SMEAR_STUDY_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TStnJetBlock.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnJetBlock.hh"
#include <iostream>

class JetSmearStudy: public TStnModule
{

	protected:
		TStnHeaderBlock* fHeaderBlock;
		TStnJetBlock* fJetBlock;

	public:
		JetSmearStudy(const char* name="JetSmearStudy", const char* title="JetSmearStudy");
		~JetSmearStudy();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnJetBlock*   GetJetBlock()  { return fJetBlock;  }

		double Jer_normG(double jet_E, double eta_det, int stat_code, int syst_code);
		double Jer_sigmaL(double jet_E, double eta_det, int stat_code, int syst_code);
		double Jer_mpvL(double jet_E, double eta_det, int stat_code, int syst_code);
		double Jer_sigmaG(double jet_E, double eta_det, int stat_code, int syst_code);
		double Jer_meanG(double jet_E, double eta_det, int stat_code, int syst_code);
		double IntegralUpLimit(double jet_E, double eta_det, int stat_code, int syst_code);
		double JERparam(int code, int i, int j);
		double JERCorrCoeff(int i, int j);
		int WhatJEReta(double eta_det);


		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();


		void BookHistograms();
		unsigned GetNptsToGenerate() const { return iNptsGenerate; }
		void NptsToGenerate(const unsigned int pt) { iNptsGenerate = pt; } 
		void SmearJet(const TLorentzVector jetvec);

	private:
		unsigned int iNptsGenerate;		// number of random pts generated per jet
		TH1F *hJetPt_Orig;
		TH1F *hJetPt_Smear;

	ClassDef(JetSmearStudy,1);
};
#endif
