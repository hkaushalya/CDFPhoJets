///////////////////////////////////////////////////////////
//  //
//  //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

// require "boost" library installed. 
#include "samantha/utils/FileHandler.hh"

ClassImp(FileHandler)

//_____________________________________________________________________________
FileHandler::FileHandler(const char* name, const char* title):
  TStnModule(name,title),
  sOutFileName("NO_NAME.sxt"),
  bWritePermit(false)
{
	std::cout << "Hello I am FileHandler module" << std::endl;
}

//_____________________________________________________________________________
FileHandler::~FileHandler() {
}
//_____________________________________________________________________________
int FileHandler::BeginJob()
{
	OutFile.open(sOutFileName.c_str());
	if (OutFile.is_open()) {
		StdOut(__FILE__,__LINE__,0,"Log File "+sOutFileName+ " opened successfully.");
		bWritePermit = true;
	} else {
		StdOut(__FILE__,__LINE__,3,"Log File "+sOutFileName+ " is not open. pl. check.");
		bWritePermit = false;
	}

	fHeaderBlock   = (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	if (! fHeaderBlock) { std::cout << "************* header error" << std::endl; }
	
}

//_____________________________________________________________________________
int FileHandler::BeginRun()
{

//	Write(">>>>>>>>>> Begin Run " fHeaderBlock->RunNumber() );

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int FileHandler::Event(int ientry)
{
	SetPassed(1);
	// now write the stuff of the previous event
	// output is always delayed by 1 event. bcos of the way the mods are run
	// in this frame work. so need to add this same steps at the end of job
	// still if this is at the top of the mods list, we can't write end job
	// summaries as this will be closed first.

	WriteToFile();
	Cleanup();
	return 0;
} // Event

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void  FileHandler::WriteToFile()
{
	for (std::vector<std::string>::iterator it = vText2File.begin(); 
			it != vText2File.end(); it++) {
		OutFile << (*it).c_str() << "\n";
	}
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void  FileHandler::WriteDirect(double num_)
{
	WriteDirect(ToStr(num_));
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void  FileHandler::WriteDirect(std::string text_)
{
	// only to write 'end job' summaries to the file
	// after the Event loop is over
	if (OutFile.is_open()) {
		OutFile << text_.c_str() << "\n";
	} else {
		StdOut(__FILE__,__LINE__,3,"Log File "+sOutFileName+ " is not open. pl. check.");
	}
}
/*-------------------------------------------------------------------*/
void  FileHandler::WriteDirect(TStnModule* mod, std::string text_)
{
	// only to write 'end job' summaries to the file
	// after the Event loop is over
	std::string modname = mod->GetName();
	WriteDirect(modname + ":: " + text_);
}
/*-------------------------------------------------------------------*/
void  FileHandler::WriteDirect(TStnModule* mod, double num_)
{
	// only to write 'end job' summaries to the file
	// after the Event loop is over
	std::string modname = mod->GetName();
	WriteDirect(modname + ":: " + ToStr(num_));
}

/*-------------------------------------------------------------------*/
void  FileHandler::Write(TStnModule* mod, std::string text_)
{
	std::string modname = mod->GetName();
	Write(modname + ":: " + text_);
}
/*-------------------------------------------------------------------*/
void  FileHandler::Write(TStnModule* mod, double num_)
{
	std::string modname = mod->GetName();
	Write(modname + ":: " + ToStr(num_));
}


/*-------------------------------------------------------------------*/
void  FileHandler::WritePass(TStnModule* mod)
{
	// call this to put a record in the log for events passed
	std::string modname = mod->GetName();
	std::string str = modname +"::Pass Run,Evt: " + ToStr(fHeaderBlock->RunNumber()) 
							+ "\t" + ToStr(fHeaderBlock->EventNumber());
	WriteDirect(str);		//write events pass immedietly. if the next module crashed,
								// i can know when/at what event it crached.
}
/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void FileHandler::Cleanup()
{
	// clear the string holder
	vText2File.clear();
} 


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int FileHandler::EndJob() {

	// write the last events stuff to file but do not close, let others write too.
	Write(">>>>>>>>>>>>>>>>>>>>>> END JOB SUMMARIES <<<<<<<<<<<<<<<<<<<<<");
	WriteToFile();
	printf("[FHD:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[FHD:01:] Log file " << GetFileName() << " was created." << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
