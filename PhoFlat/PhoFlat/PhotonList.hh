#ifndef PhotonList_HH
#define PhotonList_HH

#include <vector>

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>
/*

*/
class PhotonList {

	protected:

	public:
			PhotonList();
			
			unsigned int  pho_num;							// number of electrons so we know how much to read from arrays
			unsigned int  pho_Ntight;						// # of tight photons
			unsigned int  pho_Nloose;						// # of loose photons
			std::vector<int>   pho_Index;				//this is the index(order) in the Photon block
			std::vector<int>   pho_PhoBlockIndex;  	//this is the index(order) in the Photon block
			std::vector<float> pho_Etc;
			std::vector<float> pho_E;
			std::vector<float> pho_Px;
			std::vector<float> pho_Py;
			std::vector<float> pho_Pz;
			std::vector<int>   pho_Detector;
			std::vector<float> pho_DetEta;
			std::vector<float> pho_DetPhi;
			std::vector<float> pho_XCes;
			std::vector<float> pho_ZCes;
			std::vector<float> pho_HadEm;
			std::vector<float> pho_Chi2Mean;
			std::vector<int>   pho_N3d;
			std::vector<float> pho_Iso4;
			std::vector<float> pho_TrkPt;				// from photon block
			std::vector<float> pho_TrkIso;				// from photon block
			std::vector<float> pho_CesWireE2;
			std::vector<float> pho_CesStripE2;
			std::vector<int>   pho_PhiWedge;
			std::vector<int>   pho_NMuonStubs;			//trkless muon stubs around 30 degree cone of the photon
			std::vector<float> pho_EmTime;
			std::vector<int>   pho_TightId;
			std::vector<int>   pho_LooseId;
			std::vector<int>   pho_PhoenixId;
			std::vector<int>   pho_Halo_seedWedge;
			std::vector<int>   pho_Halo_eastNhad;
			std::vector<int>   pho_Halo_westNhad;
			std::vector<int>   pho_matchJetIndex; 	//if a matching jet is found and removed from jet list
			std::vector<float> pho_CprWgt;
			std::vector<float> pho_CprSys1;			// each vector will hold 8 syst values
			std::vector<float> pho_CprSys2;
			std::vector<float> pho_CprSys3;
			std::vector<float> pho_CprSys4;
			std::vector<float> pho_CprSys5;
			std::vector<float> pho_CprSys6;
			std::vector<float> pho_CprSys7;
			std::vector<float> pho_CprSys8;

		void Clear();

		void DumpBasic() const;
		void Dump() const;

	private:

};

#endif
