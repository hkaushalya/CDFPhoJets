#ifndef COMMON_TOOLS_HH
#define COMMON_TOOLS_HH

/*******************************************************************
*  This has a bunch of utilities that is used under Stntuple
*  framework and when I make the final hitsts.
*  Author: Sam Hewamanage <samantha@fnal.gov>
********************************************************************
*	$Id: CommonTools.hh,v 1.4 2012/10/29 19:41:08 samantha Exp $
*	$Log: CommonTools.hh,v $
*	Revision 1.4  2012/10/29 19:41:08  samantha
*	Final commit.
*	
*	Revision 1.3  2011/05/31 18:27:42  samantha
*	MOVED: The output color commnads to a separated module as this is too bulky if
*	one I want to use only color stuff.
*	
*	Revision 1.2  2010/03/18 15:34:56  samantha
*	ADDED: 1. Two overloaded GetMaxBinContent() method to find the contents of the
*	maximum bin in a hist or a set of hists.
*		2. class StringVector, which can hold a set of strings to which can be
*		printed at anytime during a job. Same info can be written to a canvas
*		and saved to the root file at the end of the job. I created this so  do
*		not have to search for the log file for a job. The root file will have
*		all the info I need about the job.
*	
*	Revision 1.1  2009/12/15 23:24:02  samantha
*	A collection of utitlity functions, color couts, Draw hists, rebin, etc.
*	
*
********************************************************************/

#include "TH1.h"
#include <string>
#include <memory>
#include <ostream>
//#include <time.h>
#include <ctime>
#include <iostream>
#include <vector>
#include "TCanvas.h"
//#include "../samantha/obj/Stuple.hh"
#include "TTree.h"
#include "TChain.h"

void NormalizeBinWidth(TH1* hist);		//normalize each bin by dividing by it.
void DrawClone(const TH1* hist, const bool log=true);		//draws clone of the hist
void DumpBin(const TH1* hist, const int iBin, const std::string sMsg="");
void DumpHist(const TH1* hist, const std::string sMsg="");
int GetBin(const TH1* hist, const double val);
void FindCommonNonZeroXrange(const TH1* h1, const TH1* h2, float& xmin, float& xmax);
void FindNonZeroXrange(const TH1* h1, int& xminbin, int& xmaxbin);
std::vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4);
TH1F* MakeVariableBinHist (const std::string name, const std::string title,
							const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4);

std::auto_ptr<TH1> FillVarBinHist (const std::vector<float>& bins, TH1 *input);
TH1* MakeVariableBinHist (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4);
void WilsonInterval(float k, float n, float &center, float &error);
void WilsonInterval(TH1* &nom, TH1* denom);

/* move all this to IOColors.hh 04-07-2010
//for color outputs

//foreground colors
inline std::ostream& red      (std::ostream &s) { s << "\033[31m"; return s; }
inline std::ostream& green    (std::ostream &s) { s << "\033[32m"; return s; }
inline std::ostream& yellow   (std::ostream &s) { s << "\033[33m"; return s; }
inline std::ostream& blue     (std::ostream &s) { s << "\033[34m"; return s; }
inline std::ostream& magenta  (std::ostream &s) { s << "\033[35m"; return s; }
inline std::ostream& cyan     (std::ostream &s) { s << "\033[36m"; return s; }
inline std::ostream& white    (std::ostream &s) { s << "\033[37m"; return s; }
inline std::ostream& black    (std::ostream &s) { s << "\033[0m" ; return s; }
inline std::ostream& clearatt (std::ostream &s) { s << "\033[0m" ; return s; }
//special
inline std::ostream& bold (std::ostream &s) { s << "\033[1m" ; return s; }   //does not seem to  work
inline std::ostream& undeline (std::ostream &s) { s << "\033[4m" ; return s; }
inline std::ostream& blink (std::ostream &s) { s << "\033[5m" ; return s; }   //does not seem to  work
//background colors
inline std::ostream& redbg      (std::ostream &s) { s << "\033[41m"; return s; }
inline std::ostream& greenbg    (std::ostream &s) { s << "\033[42m"; return s; }
inline std::ostream& yellowbg   (std::ostream &s) { s << "\033[43m"; return s; }
inline std::ostream& bluebg     (std::ostream &s) { s << "\033[44m"; return s; }
inline std::ostream& magentabg  (std::ostream &s) { s << "\033[45m"; return s; }
inline std::ostream& cyanbg     (std::ostream &s) { s << "\033[46m"; return s; }
inline std::ostream& whitebg    (std::ostream &s) { s << "\033[47m"; return s; }
inline std::ostream& blackwhite    (std::ostream &s) { s << "\033[7m" ; return s; }; //blk text on white background
*/
double GetMaxBinContent(const TH1* hist);
double GetMaxBinContent(const std::vector<TH1*> hist);
//TTree* GetStupleTree(const std::string sTreeName, Stuple& stuple);
//void PrepStupleTreeForReading(TChain *tree, Stuple *st);


class StringVector
{
	public:
		StringVector();
		~StringVector();
		inline void clear() { strVec.clear(); }
		inline unsigned size() { return strVec.size(); }
		inline void print() 
		{
			for (std::vector<std::string>::iterator it = strVec.begin(); it != strVec.end(); ++it)
				std::cout << *it << std::endl;
		};
		inline std::string at(const unsigned i) 
		{
			if (i < strVec.size()) return strVec.at(i);
			else 
			{
				std::cout << "ERROR: Out of bound request! returning empty string!" << std::endl;
				return "";
			}
		};
		TCanvas* writeToCanvas();
		std::vector<std::string>& vec() { return strVec; }
		inline void push_back(const std::string s) { strVec.push_back(s); }
		TCanvas* canvas() { return c; }

	private:
		std::vector<std::string> strVec;
		TCanvas* c;
};

class Clock
{
	public:
		Clock();
		Clock(const std::string n);
		void Start();
		void Stop();
		//void Start() {start = time_t(NULL); std::cout << "start time = " << start << std::endl;}
		//void Stop() { end = time_t(NULL);  std::cout << "end  time = " << end << std::endl;}

		void ShowElapsedSeconds() const 
		{ 
			std::cout << name << ": Elapsed " <<  time_t(NULL) - start << " (sec)" << std::endl; }
		void ShowElapsedMinutes() const 
		{ 
			std::cout << name << ": Elapased " <<  (time_t(NULL) - start)/60. << " (min)" << std::endl; 
		}
		void ShowElapsedHours() const 
		{ 
			std::cout << name << ": Elapased " <<  (time_t(NULL) - start)/3600. << " (hrs)" << std::endl; 
		}
		void ShowEndTimeHours() const; 
/*		{ 
			if (end > 0)
				std::cout << name << ": Total Time " <<  (end - start)/3600. << " (hrs)" << std::endl; 
			else 
				std::cout << "Clock has not stopped!" << std::endl;
			
			std::cout << name << ": Total Time " <<  (end - start)/3600. << " (hrs)" << std::endl; 
		}
*/
		double GetElaspsedSeconds() const { return time_t(NULL) - start; }
		double GetElaspsedMinutes() const { return (time_t(NULL) - start)/(time_t)60.; }
		double GetElaspsedHours() const { return (double)(time_t(NULL) - start)/3600.; }
		
		double EndTimeSeconds() const { return end - start; }
		double EndTimeMinutes() const { return (end - start)/60.; }
		double EndTime() const { return (end - start)/3600.; }
		
	private: 
		std::string name;
		time_t start, end;		
	

};

#endif
