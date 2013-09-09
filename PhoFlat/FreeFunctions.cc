#include "PhoFlat/FreeFunctions.hh"
#include <sstream>
#include "TLorentzVector.h"
#include "TMath.h"

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


std::string ToStr(float num)
{
	// coverts a number to string
	std::ostringstream os;
	os << num;
	return os.str();
}

void CorrectPhotonEnergy(Stuple& stuple, int JES)
{
	float EmEnergyCorr = 0.01;		// 1%

	for (int i=0; i < stuple.pho_num ; i++) {
		TLorentzVector vec(stuple.pho_Px[i], stuple.pho_Py[i], stuple.pho_Pz[i], stuple.pho_E[i]);
		
		if (JES>0) 	vec += vec * EmEnergyCorr;
		if (JES<0)	vec -= vec * EmEnergyCorr;
		stuple.pho_Px[i] = vec.Px();
		stuple.pho_Py[i] = vec.Py();
		stuple.pho_Pz[i] = vec.Pz();
		stuple.pho_E[i] = vec.E();
		stuple.pho_Etc[i] = vec.E() * TMath::Sin(vec.Theta());
	}
}


