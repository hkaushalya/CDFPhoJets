#ifndef TPDFREWEIGHT_HH
#define TPDFREWEIGHT_HH

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include <Stntuple/loop/TStnModule.hh>
#include <Stntuple/obj/TStnVertexBlock.hh>
#include <Stntuple/obj/TGenpBlock.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <string>
#include <vector>

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif


// TPDFReweight.hh contains the description for the module which computes the reweighting 
//                 due to the PDFs (following Stephen Miller's method).
// Written by O. Gonzalez (21-V-2004)

// *************************************************************
/* Minor modifications for pho+jets analysis by sam. 
 *
 * $Id: TPDFReweight.hh,v 1.3 2009/11/19 22:38:52 samantha Exp $
 * $Log: TPDFReweight.hh,v $
 * Revision 1.3  2009/11/19 22:38:52  samantha
 * ADDED: GetConsiderAlphasEffect() to see if this effect is used.
 *
 * Revision 1.2  2009/11/07 17:33:33  samantha
 * Added few accessors to get the settings used for the job that will be printed at
 * the end of the job.
 *
 * Revision 1.1  2009/10/23 16:49:40  samantha
 * Slight modifications to original to fit my needs. I needed store the weights
 * for later use and this does. Added accessors and commented out some
 * stuff I do not need.
 *
 *
 */



class TPDFReweight: public TStnModule {
	public:
		struct Hist_t {

			TH1F*     fEtGenParticles;
			TH1F*     fPtGenParticles;
			TH1F*     fEtaGenParticles;
			TH1F*     fEtDecParticles;
			TH1F*     fPtDecParticles;
			TH1F*     fEtaDecParticles;

			TH1F*     fEtGenParticles_orig;
			TH1F*     fPtGenParticles_orig;
			TH1F*     fEtaGenParticles_orig;
			TH1F*     fEtDecParticles_orig;
			TH1F*     fPtDecParticles_orig;
			TH1F*     fEtaDecParticles_orig;

			void Write (void) {
				fEtGenParticles->Write();
				fEtaGenParticles->Write();
				fEtDecParticles->Write();
				fEtaDecParticles->Write();
			}

		};

	protected:
		// pointers to the data blocks used
		TGenpBlock* fGenpBlock;
		TStnHeaderBlock *fHeaderBlock;


		// histograms filled
		Hist_t  fHist;
		// could be used as a filter 
		// by deafult everything is of
		/*
			Int_t    ZVertexVClass;  // z-vertex quality flag
			Float_t Z_max_cut;
			Int_t    _fail;  
			Int_t    _fail_zcut;  
			Int_t    _pass;
			*/

		std::string _usingMCgen;
		// Available Options:
		//        Isajet    
		//        Pythia

		float _usingScale;   // Using this Q instead of the calculated (useful for resonances)

		int _pdfref_grp;
		int _pdfref_set;

		std::string _usingPDFset;

		int _nEvents;				//total number of events processed

		int _nalternativepdfs;
		bool _hasstandarderrors;

		int _pdfalt_grp[50]; 
		int _pdfalt_set[50];

		bool _consideralphaseffect;
		int _iordalphas;  // For pythia it should be -1. For Isajet (if QCD value was changed)
		// it should be also LO (although Isajet is tricky)


		double tot_weight[50];
		double tot_weightsqr[50];

		int _indexHistogram;   // Indicates which Reweighted PDF is selected for the histograms

	public:

		TPDFReweight(const char* name="PDFReweight", const char* title="PDFReweight");
		~TPDFReweight();

		// ****** accessors
		Hist_t*        GetHist     () { return &fHist;    }
		double GetTotWeight (const unsigned ipdf=0) const { return ((ipdf<50)?tot_weight[ipdf]:0);}
		double GetTotWeightSqr (const unsigned ipdf=0) const { return ((ipdf<50)?tot_weightsqr[ipdf]:0);}
		double GetAvrWeight (const unsigned ipdf=0) const { return ((ipdf<50)?tot_weight[ipdf]/_nEvents:0);}
		double GetErrorOnAvrWeight (const unsigned ipdf=0) const;
		std::string GetMCGenerator () const {return _usingMCgen;}
		int GetNalternativePDFs() const { return _nalternativepdfs; }
		void GetReferencePDF (int &igroup,int &ipdfset) const {igroup = _pdfref_grp; ipdfset = _pdfref_set;}
		std::string GetReweightedwithPDF () const { return _usingPDFset; }
		float GetMCScale () const { return _usingScale;}
		bool GetConsiderAlphasEffect () const { return _consideralphaseffect; }

		// ****** setters
		void SetMCGenerator (std::string gen) {_usingMCgen=gen;}
		void SetMCScale (float q) {_usingScale=q;}
		void SetReferencePDF (int group,int pdfset) {_pdfref_grp=group;_pdfref_set=pdfset;}
		void SetReweightedwithPDF (std::string pdfname);
		void SetConsiderAlphasEffect (bool consideralphaseffect)
		{
			_consideralphaseffect = consideralphaseffect;
		}

		void SetAlphaSOrder (int iord) {_iordalphas=iord;}
		void SetIndexForControlHistogramming (const int idx) {_indexHistogram=idx;}

		std::vector<float> GetEventPDFWeights() const { return vPDFweights; }
		float GetEventQ() const { return Q; }
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

		/*  void  SetZVertexVClass(int vclass) {ZVertexVClass=vclass;} 

			 void SetZVertexZCut (float zcut) {Z_max_cut=zcut;}
			 */ 
		// ****** overloaded methods of 
		// TStnModule
		int     BeginJob();
		int     BeginRun();
		int     Event   (int ientry);
		int     EndJob  ();
		// ****** other methods
		void    BookHistograms();
		//void    PlotHistograms(int run_number, int slide);
		//  void    PlotHistograms(char *run_number, int slide);


	
	private:
		std::vector<float> vPDFweights; //holds PDF weight for the current event
		float Q;				//Q value of the event
		bool bNoSummary;			//control the end job summary print
		unsigned int iPassed;   // number of events pass this module
		unsigned int iProcess;  // nummber of events processed
		std::string ToStr(const double number); //converts any number to a string

		ClassDef(TPDFReweight,0)
};
#endif
