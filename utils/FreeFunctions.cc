/*
keep it simple and readable. forget about performance!
*/
/*
all free functions
*/

#include "samantha/utils/FreeFunctions.hh"
//#include <boost/lexical_cast.hpp>
#include "Stntuple/obj/TGenParticle.hh"
#include<assert.h>

enum Severity_t {
	INFO = 0,
	WARN = 1,
	ERRO = 2,
	FERR = 3,
};

std::vector<std::string> vsSeverity;
//vsSeverity.push_back("INFO");
//vsSeverity.push_back("WARN");
//vsSeverity.push_back("ERRO");
//vsSeverity.push_back("FERR");



/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2)
{
	return (c1->Energy() > c2->Energy()) ? true : false;
}


/*-------------------------------------------------------------------*/
void StdOut(const std::string& msg, const double num, const std::string& ModName)
{
	std::cout << ModName << "::"<< msg << "=" << num << std::endl;
}
/*-------------------------------------------------------------------*/
void StdOut(const std::string& msg, const std::string& ModName)
{
	std::cout << ModName << "::"<< msg << std::endl;
}


/*-------------------------------------------------------------------*/
TFolder* GetHistFolder(TStnModule* const mod,const std::string& folderName, const std::string& folderTitle)
{
	if (mod) {
		TFolder *histFolder = (TFolder*) mod->GetFolder()->FindObject("Hist");
		if (histFolder) {
			TFolder *newFolder = (TFolder*) histFolder->FindObject(folderName.c_str());
			if (!newFolder) newFolder = histFolder->AddFolder(folderName.c_str(),folderTitle.c_str());
			return newFolder;
		} else {
			std::cout << "ERROR::" << __FILE__ << "::" << __LINE__ << ":: Cannot find parent folder'Hist'. Returning NULL." << std::endl;
		}
	} else std::cout<< "ERROR::" << __FILE__ << "::" << __LINE__ << ":: Passed in TStnModule is NULL. Returning NULL." << std::endl;

return NULL;
}


/*-------------------------------------------------------------------*/
void StdErr(const std::string& fileName, const unsigned int lineNumber,
				const int severity, const std::string& message)
{
	/*
		0 = info
		1 = warning
		2 = error
		3 = fatal error
	*/
	std::string sev;
	switch (severity) {
		case 0: {
			sev = "INFO";
			break;
		}
		case 1: {
			sev = "WARN";
			break;
		}
		case 2: {
			sev = "ERRR";
			break;
		}

		case 3: {
			sev = "FERR";
			break;
		}
	}

	std::cerr << sev << "::" << fileName << "::" << lineNumber << ":: " << message << std::endl;
	//std::cout << sev << "::" << fileName << "::" << lineNumber << ":: " << message << std::endl;
	
}

/*-------------------------------------------------------------------*/
void StdOut(const std::string& fileName, const unsigned int lineNumber,
				const int severity, const std::string& message)
{
	/* 0 = info
		1 = warning
		2 = error
		3 = fatal error
	*/
		std::string sev;
	switch (severity) {
		case 0: {
			sev = "INFO";
			break;
		}
		case 1: {
			sev = "WARN";
			break;
		}
		case 2: {
			sev = "ERRR";
			break;
		}

		case 3: {
			sev = "FERR";
			break;
		}
	}
	std::cout << sev << "::" << fileName << "::" << lineNumber << ":: " << message << std::endl;
	
}


/*------good  for error reporting------------------------------------*/
void ToFile(ofstream* ofFile, const std::string& fileName, const unsigned int lineNumber,
				const int severity, const std::string& message)
{
	std::string sev;
	if (severity >=0 && severity < vsSeverity.size()) {		// if severity is out of the std values, consider
		sev = vsSeverity[severity];											// as info.
	} else sev = vsSeverity[0];

	if (ofFile) {
		*ofFile << "-----------------------------------------------------------------" << std::endl;
		*ofFile << fileName << "::" << lineNumber << std::endl;
		*ofFile << sev << "::" << message << std::endl;
	} else {
		std::cout << "-----------------------------------------------------------------" << std::endl;
		std::cout << __FILE__ << "::" << __LINE__ << ":: I can't write to a bad file pointer.!" << std::endl;
		std::cout << fileName << "::" << lineNumber << std::endl;
		std::cout << sev << "::" << message << std::endl;
	}
}

	

/*------good  just printing stuff------------------------------------*/
void ToFile(ofstream* ofFile, const std::string& message)
{
	if (ofFile) {
		*ofFile << message;
	} else {
		std::cout << message;			// i don't want to loose the summary, do i?
	}
}
void ToFile(ofstream* ofFile, const char* message)
{
	if (ofFile) {
		*ofFile << message;
	} else {
		std::cout << message;			// i don't want to loose the summary, do i?
	}
}


/*-------------------------------------------------------------------*/
void PrintInBinaryForm(const int& number)
{
	if (number < 0) {
		std::cout <<  "ERROR::" << __FILE__ << "::" << __LINE__ << "::Number given is Negative!" << std::endl;
		return;
	}
	int num = number;
	binary(num);
	std::cout << std::endl;
}


/*-------------------------------------------------------------------*/
void binary(const int& number)
{
	int remainder;
	if (number <= 1 ) {
		std::cout << number;
		return;
	}
	remainder = number%2;
	binary(number >> 1);    
	std::cout << remainder;
}


//==========================================================
double DelPhi(double p1, double p2)
{
	//i assume p1, p2 to from 0 to 2pi. but should this be a must??
	//return Phi separation from -pi to pi range
	double kPI = TMath::Pi();
	double kTWOPI = TMath::TwoPi();

	double dDelPhi = p1 - p2;

	if (dDelPhi >kPI) {
		dDelPhi -= kTWOPI;
	}
	if (dDelPhi < -kPI) {
		dDelPhi += kTWOPI;
	}
									
	return dDelPhi;
}

/*-------------------------------------------------------------------*/
std::string ToStr(float num)
{
	// coverts a number to string
/*	std::string str;
	try {
		str = boost::lexical_cast<std::string>(num);
	} 
	catch (boost::bad_lexical_cast &) {
		std::cout<<__FILE__<<"::"<<__LINE__<<":: double->str casting FAILED for value:"<<num<<std::endl;
		str = "cast failed!";
	}
	*/
	std::stringstream s;
	s << num;
	return s.str();
}

/*-------------------------------------------------------------------*/
TLorentzVector FindMatchingHEPGPar(TGenpBlock* fGenpBlock, const TLorentzVector tlObj,
									 const float fDelR, const int iPDGcode, const int iStatus,
									 const int iLoopTo, const bool bDebug)
{
	// For a given detector obj (4-vec) this will find the matching HEPG particle
	// and returns its 4-vec
	// iStatus must be 1,2, or 3. see elog #1201

	if (bDebug)
	{
		std::cout << "::::" << __FILE__<<"::"<<__FUNCTION__ << std::endl;
	}

	assert (fGenpBlock != NULL && "FreeFunctions::FindMatchingHEPGPar():: Passed fGenpBlock is null!");
	assert (fDelR>0  && "FreeFunctions::FindMatchingHEPGPar():: DelR separation must be >0!");

	if (iStatus<0 || iStatus>3)
	{
		StdOut(__FILE__,__LINE__,3,"Requiring a unknown HEPG particle status!");
		exit (1);
	}

	int Nparticles = fGenpBlock->NParticles();

	if (iLoopTo > -1 && iLoopTo < Nparticles)
	{
		Nparticles = iLoopTo;
	}

	if (bDebug)
	{
		std::cout << __FUNCTION__ << "::Values of fDelR("<< fDelR << "), PdgCode(" << iPDGcode 
						<< "), iStatus(" << iStatus << "), loopTo("<< iLoopTo << ")" << std::endl;
	}

	TLorentzVector hepgVec(0,0,0,0);
	bool bNotFound = true;

	TGenParticle *par, *mom;

	for (int i = 0 ; i < Nparticles ; i++)	
	{
		par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
			mom = fGenpBlock->Particle(im);

			if (mom != 0) {
				int par_id   = par->GetPdgCode();
				int par_stat = par->GetStatusCode();
				TLorentzVector parvec;
				par->Momentum(parvec);

				if (iStatus>0)
				{
					if (par_stat != iStatus) continue; 			// find a stable particle first. see elog#1201
				}
				if (par_id != iPDGcode) continue;

				if (tlObj.DeltaR(parvec) < fDelR) {
					if (bDebug)
					{
						std::cout << "FOUND MATCH "<< std::endl;	
						std::cout << "i,par_id,pat_stat,delR,PtRat = " << i << "\t" << par_id << "\t" << par_stat;
						std::cout << "\t" << tlObj.DeltaR(parvec) << "\t" << tlObj.Pt()/parvec.Pt() << std::endl;
					}

					hepgVec = parvec;
					bNotFound = false;
					break;
				}

			}
		} 
	}


	if (bNotFound)
	{
		StdOut(__FILE__, __LINE__,2, ":; No matching HEPG particle found!");
	}

	return hepgVec;
}

/*-------------------------------------------------------------------*/
TLorentzVector FindHEPGPar(TGenpBlock* fGenpBlock, const int iPDGcode, 
									const int iStatus, const int iLoopTo, 
									const bool bDebug)
{
	// returns its 4-vec of the first particle with the given PDG code
	// and status
	// iStatus must be 1,2, or 3. see elog #1201
	// iLoopTo can be set to search first 20 particles or so to speed 
	// up the search. Default is loop over all particles

	if (bDebug)
	{
		std::cout << "::::" << __FILE__<<"::"<<__FUNCTION__ << std::endl;
	}

	assert (fGenpBlock != NULL && "FreeFunctions::FindMatchingHEPGPar():: Passed fGenpBlock is null!");

	if (iStatus<0 || iStatus>3)
	{
		StdOut(__FILE__,__LINE__,3,"Requiring a unknown HEPG particle status!");
		exit (1);
	}

	int Nparticles = fGenpBlock->NParticles();

	if (iLoopTo > -1 && iLoopTo < Nparticles)
	{
		Nparticles = iLoopTo;
	}

	if (bDebug)
	{
		std::cout << __FUNCTION__ << "::Values of  PdgCode(" << iPDGcode 
						<< "), iStatus(" << iStatus << "), loopTo("<< iLoopTo << ")" << std::endl;
	}

	TLorentzVector hepgVec(0,0,0,0);
	bool bNotFound = true;

	TGenParticle *par, *mom;

	for (int i = 0 ; i < Nparticles ; i++)	
	{
		par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if (im >= 0) 	//_________________________________im=-1 means no mother for imcoming particles
		{
			mom = fGenpBlock->Particle(im);

			if (mom != 0)
			{
				int par_id   = par->GetPdgCode();
				int par_stat = par->GetStatusCode();
				TLorentzVector parvec;
				par->Momentum(parvec);

				if (bDebug)
				{
					std::cout << "i,par_id,par_stat = " << i << "\t" << par_id << "\t" << par_stat << std::endl;
				}
				
				if (iStatus>0)
				{
					if (par_stat != iStatus) continue; 			// find a stable particle first. see elog#1201
				}
				if (par_id != iPDGcode) continue;

				if (bDebug)
				{
					std::cout << "FOUND MATCH "<< std::endl;	
					//std::cout << "i,par_id,par_stat = " << i << "\t" << par_id << "\t" << par_stat << std::endl;
				}

				hepgVec = parvec;
				bNotFound = false;
				break;

			}
		} 
	}

	if (bNotFound)
	{
		StdOut(__FILE__, __LINE__,2, ":; No matching HEPG particle found!");
	}

	return hepgVec;
}


//to sort a TLorentzVector in Pt descending order using std::sort
bool CompareTLorentzVectorPt(TLorentzVector tl1, TLorentzVector tl2) { return (tl1.Pt()>tl2.Pt()); }
