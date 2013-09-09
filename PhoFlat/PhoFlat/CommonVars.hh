#ifndef CommonVars_HH
#define CommonVars_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>
/*
*/
class CommonVars {

	protected:

	public:
			CommonVars();

			int evt_McFlag;		// =1 if MC , else 0
			int evt_RunNumber;
			int evt_EventNumber;
			int tri_pho25iso;
			int tri_pho50;
			int tri_pho70;

			int   vtx_N;
			int   vtx_NClass12;
			float vtx_z;				// highest sum pt vertex Z
			int   vtx_Ntracks;		// in highest sum pt vertex
			float vtx_SumPt;			// highest

			float met_Met;				// corrected MET
			float met_RawMet;				// raw MET
			float met_MetX;				// corrected MET
			float met_MetY;				// corrected MET
			float met_SumEt;			// corrected
			float met_Ht;				// corrected
			float met_MetPhi;			// corrected
			float met_Gen_d;			//defualt generated met 
			float met_Gen_m;			//defualt generated met  - dev
			float met_Gen_p;			//defualt generated met  + dev
			float met_Gen_mUn;			//defualt generated met  - devUn
			float met_Gen_pUn;			//defualt generated met  + devUn

			void Dump();
			void PrintHeader() const;

	private:

};

#endif
