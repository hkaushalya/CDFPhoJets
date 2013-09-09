#include "PhoFlat/ElectronList.hh"
#include <sstream>
#include <iomanip>
#include <iostream>
//#include "../RootTools/CommonTools.hh"
#include "../RootTools/IOColors.hh"

ElectronList::ElectronList()
{
	ele_num = 0;
	ele_Ntight = 0;
	ele_Nloose = 0;
}

void ElectronList::Clear()
{
	ele_num = 0;
	ele_Ntight = 0;
	ele_Nloose = 0;
	ele_Index.clear();
	ele_PhoBlockIndex.clear();
	ele_EleBlockIndex.clear();
	ele_Etc.clear();
	ele_E.clear();
	ele_Px.clear();
	ele_Py.clear();
	ele_Pz.clear();
	ele_Detector.clear();
	ele_DetEta.clear();
	ele_DetPhi.clear();
	ele_XCes.clear();
	ele_ZCes.clear();
	ele_HadEm.clear();
	ele_Chi2Mean.clear();
	ele_N3d.clear();
	ele_Iso4.clear();
	ele_TrkIso.clear();
	ele_CesWireE2.clear();
	ele_CesStripE2.clear();
	ele_PhiWedge.clear();
	ele_NMuonStubs.clear();
	ele_EmTime.clear();
	ele_PhoenixId.clear();
	ele_Halo_seedWedge.clear();
	ele_Halo_eastNhad.clear();
	ele_Halo_westNhad.clear();
	ele_matchJetIndex.clear();
	ele_Ntracks.clear();
	ele_Emfr.clear();
	ele_EoverP.clear();
	ele_TrackPt.clear();
	ele_TrackBcPt.clear();
	ele_TrackPhi.clear();
	ele_Nssl.clear();
	ele_Nasl.clear();
	ele_TightId.clear();
	ele_LooseId.clear();
	ele_ConversionId.clear();
}

void ElectronList::Delete(const int i)
{
//deletes an element
std::stringstream msg;
msg << __FILE__ <<":" << __FUNCTION__ << ": Element index given is out of range. Nothing is deleted!"; 
assert ((i>=0 && i< (int)ele_num) && msg.str().c_str()); 


	--ele_num;
	if (ele_TightId.at(i) == 0) --ele_Ntight;
	if (ele_LooseId.at(i) == 0) --ele_Nloose;
	
	
	ele_Index.erase(ele_Index.begin()+i);
	ele_PhoBlockIndex.erase(ele_PhoBlockIndex.begin()+i);
	ele_EleBlockIndex.erase(ele_EleBlockIndex.begin()+i);
	ele_Etc.erase(ele_Etc.begin()+i);
	ele_E.erase(ele_E.begin()+i);
	ele_Px.erase(ele_Px.begin()+i);
	ele_Py.erase(ele_Py.begin()+i);
	ele_Pz.erase(ele_Pz.begin()+i);
	ele_Detector.erase(ele_Detector.begin()+i);
	ele_DetEta.erase(ele_DetEta.begin()+i);
	ele_DetPhi.erase(ele_DetPhi.begin()+i);
	ele_XCes.erase(ele_XCes.begin()+i);
	ele_ZCes.erase(ele_ZCes.begin()+i);
	ele_HadEm.erase(ele_HadEm.begin()+i);
	ele_Chi2Mean.erase(ele_Chi2Mean.begin()+i);
	ele_N3d.erase(ele_N3d.begin()+i);
	ele_Iso4.erase(ele_Iso4.begin()+i);
	ele_TrkIso.erase(ele_TrkIso.begin()+i);
	ele_CesWireE2.erase(ele_CesWireE2.begin()+i);
	ele_CesStripE2.erase(ele_CesStripE2.begin()+i);
	ele_PhiWedge.erase(ele_PhiWedge.begin()+i);
	ele_NMuonStubs.erase(ele_NMuonStubs.begin()+i);
	ele_EmTime.erase(ele_EmTime.begin()+i);
	ele_PhoenixId.erase(ele_PhoenixId.begin()+i);
	ele_Halo_seedWedge.erase(ele_Halo_seedWedge.begin()+i);
	ele_Halo_eastNhad.erase(ele_Halo_eastNhad.begin()+i);
	ele_Halo_westNhad.erase(ele_Halo_westNhad.begin()+i);
	ele_matchJetIndex.erase(ele_matchJetIndex.begin()+i);
	ele_Ntracks.erase(ele_Ntracks.begin()+i);
	ele_Emfr.erase(ele_Emfr.begin()+i);
	ele_EoverP.erase(ele_EoverP.begin()+i);
	ele_TrackPt.erase(ele_TrackPt.begin()+i);
	ele_TrackBcPt.erase(ele_TrackBcPt.begin()+i);
	ele_TrackPhi.erase(ele_TrackPhi.begin()+i);
	ele_Nssl.erase(ele_Nssl.begin()+i);
	ele_Nasl.erase(ele_Nasl.begin()+i);
	ele_TightId.erase(ele_TightId.begin()+i);
	ele_LooseId.erase(ele_LooseId.begin()+i);
	ele_ConversionId.erase(ele_ConversionId.begin()+i);


}


void ElectronList::DumpBasic() const
{
	std::cout << "========= ElectronList Basic Info =========" << std::endl;
	std::cout << "ele_num    = " << ele_num << std::endl;
	std::cout << "ele_Ntight = " << ele_Ntight << std::endl;
	std::cout << "ele_Nloose = " << ele_Nloose << std::endl;
	std::cout << "===========================================" << std::endl;
};

void ElectronList::Dump() const
{
	DumpBasic();

	std::cout << blue << std::setw(7) << "Index"
				<< std::setw(8) << "EleEt"
				<< std::setw(8) << "Eta"
				<< std::setw(8) << "TightId"
				<< std::setw(8) << "LooseId"
				<< std::setw(12) << "MatchJetInd"
				<< std::setw(10) << "EleBlkInd"
				<< clearatt << std::endl;

	for (int i = 0 ; i<ele_num; ++i)
	{
		std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) 
			<< std::setw(8) << ele_Etc.at(i)
			<< std::setw(8) << ele_DetEta.at(i)
			<< std::setw(8) << ele_TightId.at(i)
			<< std::setw(8) << ele_LooseId.at(i)
			<< std::setw(12) << ele_matchJetIndex.at(i)
			<< std::setw(10) << ele_EleBlockIndex.at(i)
			<< std::endl;
	}
	std::cout << "===========================================" << std::endl;
}
