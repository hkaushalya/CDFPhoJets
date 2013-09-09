#ifndef EMTIMECORRECTION_HH
#define EMTIMECORRECTION_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"

#include <map>
#include <iostream>
#include <fstream>

#include <Stntuple/loop/TStnModule.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/obj/TStnEmTimingBlock.hh>
#include <Stntuple/obj/TCalDataBlock.hh>
#include <CalorGeometry/CalConstants.hh>
#include "samantha/utils/FreeFunctions.hh"
#endif

class EmTimeCorrection {
	public:
		
		typedef enum{READM=0, WRITEM=1, DUMMY=2} OperMode;
		typedef enum{EXACT, NEAREST} RunFinder;
		typedef std::map<int,TH1*> RunHistMap;
		typedef RunHistMap::iterator RHMI;
		typedef std::pair<int,TH1*> RunHistPair;
		EmTimeCorrection(OperMode, std::string);
		typedef std::map<int,float> RunCorValue;
		typedef RunCorValue::iterator RCVI;
		typedef std::pair<int,float> RunCorPair;

		void setHisto(int, float,float);
		void setMode(OperMode mode){ _opmode = mode;}
		void addEvent(int runNumber, float val);
		void endJob();

		inline OperMode mode() {return _opmode;}
		inline std::string varName() {return _varName;}
											    
		bool readInFile();
		bool hasRunCorrection(int run);
		float getRunCorrection(int run);

	protected:
		OperMode _opmode;
		int _curentRun;
		TH1* _curentHist;

		int _printCount;
		int _histNbins;
		float _histXmin, _histXmax;
		RunHistMap _runHists;
		std::string _varName;

		RunCorValue _runCors;
		RunFinder _rfinder;
		float _curentValue;
		bool _hasCor;
		EmTimeCorrection() {};
};
#endif
