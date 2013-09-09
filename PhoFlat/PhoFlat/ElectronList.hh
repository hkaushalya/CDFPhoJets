#ifndef ElectronList_HH
#define ElectronList_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>
/*

*/
class ElectronList {

	protected:

	public:
			ElectronList();

			//electrons collections
			unsigned int  ele_num;							// number of electrons so we know how much to read from arrays
			unsigned int  ele_Ntight;
			unsigned int  ele_Nloose;
			std::vector<int>   ele_Index;				//this is the index(order) in the Photon block
			std::vector<int>   ele_PhoBlockIndex;		//this is the index(order) in the Photon block
			std::vector<int>   ele_EleBlockIndex;		//this is the index(order) in the Electron block
			std::vector<float> ele_Etc;
			std::vector<float> ele_E;
			std::vector<float> ele_Px;
			std::vector<float> ele_Py;
			std::vector<float> ele_Pz;
			std::vector<int>   ele_Detector;
			std::vector<float> ele_DetEta;
			std::vector<float> ele_DetPhi;
			std::vector<float> ele_XCes;
			std::vector<float> ele_ZCes;
			std::vector<float> ele_HadEm;
			std::vector<float> ele_Chi2Mean;
			std::vector<int>   ele_N3d;
			std::vector<float> ele_Iso4;
			std::vector<float> ele_TrkIso;				//from photon block
			std::vector<float> ele_CesWireE2;
			std::vector<float> ele_CesStripE2;
			std::vector<int>   ele_PhiWedge;
			std::vector<int>   ele_NMuonStubs;			//trkless muon stubs around 30 degree cone of the photon
			std::vector<float> ele_EmTime;
			std::vector<int>   ele_PhoenixId;
			std::vector<int>   ele_Halo_seedWedge;
			std::vector<int>   ele_Halo_eastNhad;
			std::vector<int>   ele_Halo_westNhad;
			std::vector<int>   ele_matchJetIndex;		//if a matching jet is found and removed from jet list
			std::vector<int>   ele_Ntracks;				//these are from TStnElectron
			std::vector<float> ele_Emfr;
			std::vector<float> ele_EoverP;
			std::vector<float> ele_TrackPt;
			std::vector<float> ele_TrackBcPt;
			std::vector<float> ele_TrackPhi;
			std::vector<int>   ele_Nssl;
			std::vector<int>   ele_Nasl;
			std::vector<int>   ele_TightId;
			std::vector<int>   ele_LooseId;
			std::vector<int>   ele_ConversionId;

		void Clear();
		void Delete(const int i); //to delete an item from the list
		void DumpBasic() const;
		void Dump() const;

	private:

};

#endif
