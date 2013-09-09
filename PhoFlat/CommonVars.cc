#include "PhoFlat/CommonVars.hh"
#include <iomanip>
#include <iostream>

CommonVars::CommonVars()
{
	evt_McFlag 		= -1;
	evt_RunNumber 	= -1;
	evt_EventNumber 	= -1;
	tri_pho25iso 		= -1;
	tri_pho50 			= -1;
	tri_pho70 			= -1;

	vtx_N 		= -1;
	vtx_NClass12= -1;
	vtx_z 		= -1.0;
	vtx_Ntracks = -1;
	vtx_SumPt 	= -1.0;

	met_Met 		= -1.0;	
	met_RawMet 	= -1.0;	
	met_MetX		= -1.0;	
	met_MetX		= -1.0;	
	met_SumEt 	= -1.0;
	met_Ht 		= -1.0;	
	met_MetPhi 	= -1.0;
	met_Gen_d 	= -1.0;	
	met_Gen_m 	= -1.0;
	met_Gen_p 	= -1.0;
	met_Gen_mUn = -1.0;
	met_Gen_pUn = -1.0;

}

void CommonVars::Dump()
{
	std::cout << "========= Event Variable Info =============" << std::endl;
	std::cout << 
				"MC Flag     = " << std::setw(8) << evt_McFlag << "\n" <<
				"Run #       = " << std::setw(8) << evt_RunNumber << "\n" <<
				"Evt #       = " << std::setw(8) << evt_EventNumber << "\n" <<
				"Trig. 25iso/50/70 = " << std::setw(5) << tri_pho25iso
								<< "/" << tri_pho50 << "/" << tri_pho70 << "\n" <<
				"Nvtx/NvtxC12= " << std::setw(7) << vtx_N << "/" << vtx_NClass12 << "\n" <<
				"Vtx z       = " << std::setw(8) << vtx_z << "\n" <<
				"SumEt       = " << std::setw(8) << met_SumEt << "\n"	<<
				"Ht          = " << std::setw(8) << met_Ht << "\n"	<<
				"MEt         = " << std::setw(8) << met_Met << "\n" <<
				"MEtRaw      = " << std::setw(8) << met_RawMet << "\n" <<
				"MEtX        = " << std::setw(8) << met_MetX << "\n" <<
				"MEtY        = " << std::setw(8) << met_MetY << "\n" <<
				"MEt phi     = " << std::setw(8) << met_MetPhi << std::endl;

}
void CommonVars::PrintHeader() const
{
	std::cout << " ============== Run, Event = " << evt_RunNumber << ", " << evt_EventNumber << std::endl;
}
