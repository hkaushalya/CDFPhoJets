#include "PhoFlat/PhotonList.hh"
//#include "../RootTools/CommonTools.hh"
#include "../RootTools/IOColors.hh"
#include <iostream>
#include <iomanip>
#include <memory>

PhotonList::PhotonList()
{
	pho_num = 0;
	pho_Ntight = 0;
	pho_Nloose = 0;
}

void PhotonList::Clear()
{
	pho_num = 0;
	pho_Ntight = 0;
	pho_Nloose = 0;
	pho_Index.clear();
	pho_PhoBlockIndex.clear();
	pho_Etc.clear();
	pho_E.clear();
	pho_Px.clear();
	pho_Py.clear();
	pho_Pz.clear();
	pho_Detector.clear();
	pho_DetEta.clear();
	pho_DetPhi.clear();
	pho_XCes.clear();
	pho_ZCes.clear();
	pho_HadEm.clear();
	pho_Chi2Mean.clear();
	pho_N3d.clear();
	pho_Iso4.clear();
	pho_TrkPt.clear();
	pho_TrkIso.clear();
	pho_CesWireE2.clear();
	pho_CesStripE2.clear();
	pho_PhiWedge.clear();
	pho_NMuonStubs.clear();
	pho_EmTime.clear();
	pho_TightId.clear();
	pho_LooseId.clear();
	pho_PhoenixId.clear();
	pho_Halo_seedWedge.clear();
	pho_Halo_eastNhad.clear();
	pho_Halo_westNhad.clear();
	pho_matchJetIndex.clear();
	pho_CprWgt.clear();
	pho_CprSys1.clear();			// each vector will hold 8 syst values
	pho_CprSys2.clear();
	pho_CprSys3.clear();
	pho_CprSys4.clear();
	pho_CprSys5.clear();
	pho_CprSys6.clear();
	pho_CprSys7.clear();
	pho_CprSys8.clear();
}

void PhotonList::DumpBasic() const
{
	std::cout << "========= PhotonList Basic Info ===========" << std::endl;
	std::cout << "pho_num    = " << pho_num << std::endl;
	std::cout << "pho_Ntight = " << pho_Ntight << std::endl;
	std::cout << "pho_Nloose = " << pho_Nloose << std::endl;
	std::cout << "===========================================" << std::endl;
};

void PhotonList::Dump() const
{
	DumpBasic();

	std::cout << blue 
				<< std::setw(7) << "Index"
				<< std::setw(8) << "PhoEt"
				<< std::setw(8) << "PhoE"
				<< std::setw(8) << "Eta"
				<< std::setw(8) << "TightId"
				<< std::setw(8) << "LooseId"
				<< std::setw(12) << "MatchJetInd"
				<< std::setw(10) << "PhoBlkInd"
				<< clearatt << std::endl;

	for (unsigned int i = 0 ; i<pho_num; ++i)
	{
		std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) 
			<< std::setw(7) << i
			<< std::setw(8) << pho_Etc.at(i)
			<< std::setw(8) << pho_E.at(i)
			<< std::setw(8) << pho_DetEta.at(i)
			<< std::setw(8) << pho_TightId.at(i)
			<< std::setw(8) << pho_LooseId.at(i)
			<< std::setw(12) << pho_matchJetIndex.at(i)
			<< std::setw(10) << pho_PhoBlockIndex.at(i)
			<< std::endl;
	}
	std::cout << "===========================================" << std::endl;
}
