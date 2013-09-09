/*{{{*/
/*  $Id: ExoticEvts.C,v 1.1 2013/02/28 03:36:46 samantha Exp $
 *  $Log: ExoticEvts.C,v $
 *  Revision 1.1  2013/02/28 03:36:46  samantha
 *  Final commit. no checks.! these were never commited to cvs before!
 *
 */

/*}}}*/

#include <iostream>
#include <sstream>
#include "TBenchmark.h"
#include "TCanvas.h"
#include <iomanip>
#include "TROOT.h"
//for debugging
#include <exception>
#include <stdexcept> // out_of_range exception
#include "PhoFlat/Stuple.hh"
#include "PhoFlat/ReadInAll.hh"
#include "../RootTools/IOColors.hh"
#include "TFile.h"

bool ExoticEvent(const Stuple& st)
{
	if (st.met_Met<900) return false;
	if (st.jet_NJet15<1) return false;
	//if (st.met_Met<150) return false;
	//if (st.pho_Ntight > 0 && st.ele_Ntight > 0 > st.met_Ht>400)
	if (st.pho_Ntight> 0)
	{
		return true;
		for (int i=0; i< st.pho_num ; ++i)
		{
			if (st.pho_TightId[i] != 0) continue;
			if (fabs(st.pho_EmTime[i])> 4.8) continue;
			if (st.pho_Etc[i] < 30.) continue;
			std::cout << "pho et = " << st.pho_Etc[i] << std::endl;
			std::cout << " *********************************************************** " << std::endl;
			std::cout << "Run, Evt = " << st.evt_RunNumber << ", " << st.evt_EventNumber << std::endl;
			return true;
		}
	}
	return false;
}
 
//------ main -------------------------------------------------------
void ExoticEvts(TChain *ch, int iRunEvents)
{
	TBenchmark timer;
	timer.Start("phoana_time");
	
	
	if (ch == NULL) {
		std::cout << "NULL chain!" << std::endl;
		assert(false);
	}
	if (ch->GetEntries()<1)
	{
		std::cout << "No entries found!" << std::endl;
		return;
	}

	gROOT->Reset();
	Stuple stuple;
	ReadInAll *readIn = new ReadInAll(ch, &stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);
	readIn->Init();
	
	unsigned int iENTRIES = (unsigned int) ch->GetEntries();
	unsigned int iEvt2Process = 0; 	//_____________________________________________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;			//______________________________________________ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//______________________________________________ default

	std::cout << " ==> Entries found :: " << iENTRIES << std::endl;

	//do a check for the MC. this is a common mistake I make
	
	
	  if (stuple.evt_McFlag== 1) 
	  {
		  std::cout << red  << "THIS IS NOT DATA!" 
			  << " returning!" << clearatt << std::endl;
		  assert (false);
	  }

	unsigned int iEvtProc = 0;			//______________________________________________ number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;
	std::cout << "Init done. Ready, Set , GO... " << std::endl; 
	
	timer.Start("looptimer");
	
	int iCurrTree = -1;
	unsigned iProgressBy = (unsigned) (iEvt2Process/10.);
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
		ch->GetEntry(iEvtProc);
		if (iCurrTree != ch->GetTreeNumber())
		{
			std::cout << green << "Opening file " << ch->GetFile()->GetName() << clearatt << std::endl;
			iCurrTree = ch->GetTreeNumber();
		}

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed [" << (int)( iEvtProc/(double)iEvt2Process * 100)<< "%] ";
			timer.Show("looptimer");
			timer.Start("looptimer");
		}

		  //___________________________________________________________________________ drop the first 400pb-1
		if (stuple.evt_McFlag == 0) {  //______________________________________________ if data
			if (stuple.evt_RunNumber < 190851) {
				iFirst400++;
				continue;
			}
		}
	
		if (ExoticEvent(stuple)) 
		{
			std::cout << " Entry # = " << iEvtProc << std::endl;
			std::cout << green << " ----------  EXOTIC EVENT: Run, Event: " << stuple.evt_RunNumber 
						<< ", " << stuple.evt_RunNumber << " ------------------ " << std::endl;
			std::cout << std::setw(8) << "npho[T/L]" << std::setw(8) <<  "nele[T/L]" << std::setw(6) <<  "Njet" << std::setw(10) <<  "Ht" <<std::setw(10) <<  "Met" << std::endl;
			std::cout << std::setw(1) << stuple.pho_num << "[" << stuple.pho_Ntight << "/" << stuple.pho_num - stuple.pho_Ntight << "]"
						 << std::setw(8) << stuple.ele_num << "[" << stuple.ele_Ntight << "/" << stuple.ele_num - stuple.ele_Ntight << "]"
				       << std::setw(4) << stuple.jet_NJet15 << std::setw(10) << stuple.met_Ht << std::setw(10) << stuple.met_Met << std::endl; 
			std::cout << blue << std::setw(3) << "obj" << std::setw(3) << "ind"  << std::setw(8) << "TightId" <<  std::setw(8) << "Et" <<  std::setw(10) << "eta" << std::endl;
			for (unsigned int i = 0; i < stuple.pho_num; ++i)
			{
				std::cout << cyan << std::setw(3) << "g" <<std::setw(3) << i <<  std::setw(8) << stuple.pho_TightId[i] <<  std::setw(8) << stuple.pho_Etc[i] <<  std::setw(10) << stuple.pho_DetEta[i] << std::endl;
			}
			for (unsigned int i = 0; i < stuple.ele_num; ++i)
			{
				std::cout << std::setw(3) << "e" <<std::setw(3) << i <<  std::setw(8) << stuple.ele_TightId[i] <<  std::setw(8) << stuple.ele_Etc[i] <<  std::setw(10) << stuple.ele_DetEta[i] << std::endl;
			}
			for (unsigned int i = 0; i < stuple.jet_num; ++i)
			{
				std::cout << red << std::setw(3) << "jet" <<std::setw(3) << i <<  std::setw(8) <<  "" <<  std::setw(8) << stuple.jet_Pt[i] <<  std::setw(10) << stuple.jet_DetEta[i] << std::endl;
			}


		} else continue;
		
		std::cout << clearatt << std::endl;
		
	}
	
	//print out the job summary
	std::cout << magenta << "======== PhotonJets =======================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Entries request  = " << iRunEvents << " (" << iRunEvents/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Events Processed = " << iEvtProc << " (" << iEvtProc/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << red << "__________________________________________" << std::endl;

	std::cout << "======= END PhotonJets ====================" << std::endl;

	timer.Show("phoana_time");
} // Main



