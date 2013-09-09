#ifndef JetList_HH
#define JetList_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>
/*

*/
class JetList {

	protected:

	public:
			JetList();

			unsigned int  jet_num;							// number of jets so we know how much to read from arrays
			unsigned int  jet_NJet15;
			std::vector<int>   jet_Index;				//original index in jet block
			std::vector<float> jet_Pt;					//easy to get from JetFilter
			std::vector<float> jet_E;
			std::vector<float> jet_Px;
			std::vector<float> jet_Py;
			std::vector<float> jet_Pz;
			std::vector<float> jet_DetEta;
			std::vector<float> jet_DetPhi;
			std::vector<float> jet_HadEm;
			std::vector<float> jet_Emfr;
			std::vector<int>   jet_Ntowers;
			std::vector<int>   jet_Ntracks;
			std::vector<int>   jet_SeedIPhi;
			std::vector<int>   jet_SeedIEta;
			std::vector<int>   jet_EmTime;
			std::vector<int>	 jet_SecVtxTag;
			std::vector<float> jet_SecVtxppb;
			std::vector<float> jet_SecVtxnpb;
			std::vector<float> jet_SecVtxTrkmass;

			unsigned int  jet_raw_num;
			std::vector<int>   jet_raw_Index;
			std::vector<float> jet_raw_Pt;
			std::vector<float> jet_raw_E;
			std::vector<float> jet_raw_Px;
			std::vector<float> jet_raw_Py;
			std::vector<float> jet_raw_Pz;


			void Dump() const;
			void DumpAll() const;
			void DumpRawJets() const;
			void Clear();
			void ApplyEtEtaCuts(const float fMinPt, const float fMaxPt,
									const float fMinEta, const float fMaxEta);
	private:

};

#endif
