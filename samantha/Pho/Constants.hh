///////////////////////////////////////////////////////////
//  //
//  //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#ifndef CONSTANTS_HH
#define CONSTANTS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#endif

class Constants: public TStnModule {

	protected:

	public:
		int iLeadJetIndex;
		int iSecondLeadIndex;
		int iTightPhotonIDvars;
		int iLoosePhotonIDvars;
		int iTightElectronIDvars;
		int iLooseElectronIDvars;
	
	private:
	
	ClassDef(Constants,1)
};

#endif
