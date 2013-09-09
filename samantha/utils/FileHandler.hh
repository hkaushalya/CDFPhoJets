///////////////////////////////////////////////////////////
// This will be the common file handler for all mods.    //
// Once a intance is added to the job list, there will   //
// be a log file                                         //
// which everybody can write. this MUST be listed in top //
// before all mods.                                      //
// first have an std::string vector with all outputs,    //
// once an event is processes, write it to the file, and //
// flush before next event to avoid possible memory      //
// problems.
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#ifndef FILEHANDLER_HH
#define FILEHANDLER_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "samantha/utils/FreeFunctions.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#endif

class FileHandler: public TStnModule {

	protected:
		TStnHeaderBlock *fHeaderBlock;


	public:
		FileHandler(const char* name="FileHandler", const char* title="FileHandler");
		~FileHandler();

		// ****** accessors
		TStnHeaderBlock* GetHeaderBlock() { return fHeaderBlock; }

		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void Cleanup();
		
		void SetFile(std::string name_) { sOutFileName = name_; }
		std::string GetFileName() const { return sOutFileName; }
		void Write(std::string text_) { vText2File.push_back(text_); }
		void Write(double text_) 		{ vText2File.push_back(ToStr(text_)); }
		void WriteDirect(std::string text_);	//only to write end jobs summaries
		void WriteDirect(double num_);			//only to write end jobs summaries
		void WriteToFile(); 							// dump the vector content to file.
		void WritePass(TStnModule*);
		void Write(TStnModule* ,std::string text_);
		void Write(TStnModule*, double text_);
		void WriteDirect(TStnModule*, std::string text_);	//only to write end jobs summaries
		void WriteDirect(TStnModule*, double num_);			//only to write end jobs summaries
		void CloseFile() { OutFile.close(); }	

	private:
		ofstream OutFile;					// the log file out stream
		std::string sOutFileName;		//log file name
		std::vector<std::string> vText2File; // text to be written to file
														 // all mods out puts first comes here
		bool bWritePermit;				// makes sure file/s is open to write
												// also I can stop writing without much trouble.
		
	ClassDef(FileHandler,1)
};

#endif
