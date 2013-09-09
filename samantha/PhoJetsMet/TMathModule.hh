/* this is SAM's stuff */
#ifndef TMATHMODULE_HH
#define TMATHMODULE_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <iostream>
#include "cmath"
#include "Stntuple/obj/TStnPhoton.hh"
#include "Stntuple/obj/TStnJet.hh"

#endif

class TMathModule {

	public:
		TMathModule();
		~TMathModule();
		
		void Init();
		void SetMEtaPhi(double,double,double,double,double);
		void SetEta(double,double);
		void SetPhi(double,double);
		void SetMargin(double);
		void SetDelEta();
		void SetDelPhi();
		void SetDelR();
		double GetDelPhi() { return DelPhi; }
		double GetDelPhi(double, double);
		double GetDelEta() { return DelEta; }
		double GetDelEta(double, double);
		double GetDelR(double delPhi, double delEta);
		double GetDelR() { return DelR; }
		double GetDelR(double phi_1, double eta_1, double phi_2, double eta_2);
		double GetDelR(TStnPhoton*, TStnJet*);
		double GetDelR(TStnJet*, TStnJet*);
		double GetDelEta(TStnJet*, TStnJet*);
		double GetDelPhi(TStnJet*, TStnJet*);
		bool MatchEtaPhi(double, double,double,double,double);
		bool MatchEtaPhi(TStnPhoton*,TStnJet*,double margin);
		void Print();
		void binary(int);

	private:

		//class variables
		double phi_1,phi_2, eta_1,eta_2;
		double DelPhi, DelEta, DelR;
		double margin;

	ClassDef(TMathModule,1);
};

#endif
