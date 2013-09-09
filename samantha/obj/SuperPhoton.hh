#ifndef SUPERPHOTON_HH
#define SUPERPHOTON_HH


// Author: Samantha Hewamanage <samantha@fnal.gov>
/*
10-02-2007 - sam
This will hold all the information about the photon.
Almost all the information in the TStnPhoton class
and some more stuff like Tight Photon, Cosmic, Beam Halo, ... tags.

I don't need to copy most of these from TStnPhoton. Just keeping a
pointer to a TStnPhoton and use that to get the information as
needed would be fine. If I do that, I'll have to put additional
'if' to check for NULL pointer usage, just be safe and for easy 
debugging. In a way, even if there is little bit of copying here,
since all my later modules will be talking to this object constantly,
it may actually be little more efficient as it will be one step intead
of two to get the required info.

So let keep this as it is. May be later I can reduce the copying.

*/


#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include <string>
#include "samantha/utils/FreeFunctions.hh"

class SuperPhoton {

	protected:


	public:
		SuperPhoton();
		~SuperPhoton();

		// time window widths are set by user. so they are not exact 
		// as in the comments below
		enum TimeWindowBits_t {  // bit i will be set if it is in that time window
			kEarlyTimeBit     = 0x1 << 0,		// time < -30ns
			kBeamHaloTimeBit  = 0x1 << 1,		// time around Beam Halo >-12 & < -6ns
			kInTimeBit    		= 0x1 << 2,		// around >-5.2 & <+5.2 ns
			kLateTimeBit		= 0x1 << 3		// time > +30 ns
		};

		void Print(Option_t *opt="a", const int iPrtLevel=10);		// print its info

		//accessors

		// needed these for JetFilterModule to remove electrons properly.
		// It uses the index from the Electron Block. not from photon Block.
		// -- sam, 05-13-2008
		int GetPhotonBlockIndex() 	const { return iPhoBlockIndex;}
		int GetElectronBlockIndex() 	const { return iEleBlockIndex;}
		// So I'll make this obsolete.
		int GetIndex() 				const { 
			std::cout << __FILE__ <<"::" << __LINE__ << ":: This function is obsolete. See header file for the new accessors." << std::endl;
			return -1;
		}
		
		TStnPhoton* GetPhoton() 	const { return aPhoton    ;}				
		TLorentzVector GetCorVec() const { return tlCorVec   ;}
		TLorentzVector GetRawVec() const { return tlRawVec	  ;}
		int 	GetDetector() 			const { return iDetector  ;}
		float GetDetEta() 			const { return fDetEta    ;}
		float GetDetPhi() 			const { return fDetPhi    ;}
		float GetEtRaw() 				const { return fEtRaw     ;}			
		float GetECorr() 				const { return fECorr     ;}	
		float GetEtCorr() 			const { return fEtCorr	  ;}
		float GetXCes() 				const { return fXCes		  ;}		
		float GetZCes() 				const { return fZCes		  ;}
		float GetHadEm() 				const { return fHadEm     ;}
		float GetChi2Mean() 			const { return fChi2Mean  ;}
		int 	GetN3d() 				const { return iN3d       ;}
		float GetIso4() 				const { return fIsoEtCorr ;}
		float GetTrkPt() 				const { return fTrkPt     ;}	
		float GetTrkPt2() 			const { return fTrkPt2     ;}	
		float GetTrkIso() 			const { return fTrkIso    ;}
		float GetCesWireE2() 		const { return fCesWireE2 ;}
		float GetCesStripE2() 		const { return fCesStripE2;}
		int   GetLoosePhotonId()		const { return iLoosePhotonId   ;}	
		int   GetTightPhotonId() 		const { return iTightPhotonId   ;}
		int   GetLooseElectronId()	const { return iLooseElectronId;}			// this is photon-like electron id cuts
		int   GetStdLooseElectronId()	const { return iStdLooseElectronId;}	// this is std. electron id cuts
		int   GetTightElectronId()	const { return iTightElectronId;}
		int   GetStdTightElectronId()	const { return iStdTightElectronId;}
		int   GetBeamHaloId() 			const { return iIsBeamHalo;}
		int   GetCosmicId() 			const { return iIsCosmic  ;}
		int   GetConversionId() 		const { return iIsConversion;}
		int   GetPhoenixId() 			const { return iIsPhoenix ;}
		int   GetPMTSpikeId() 	      const { return iIsPMTSpike;}	
		float GetHadTdcTime() 		const { return fHadTdcTime;}
		float GetSinTheta() 			const { return fSinTheta  ;}
		float GetEmTime()          const { return fEmTime    ;}	
		int   GetTimeWindowBit()     const { return iTimeWindow;}
		bool  IsEarly()        const { return (iTimeWindow == kEarlyTimeBit) ? true : false ;}
		bool  IsTimeBeamHalo() const { return (iTimeWindow == kBeamHaloTimeBit) ? true : false ;}
		bool  IsInTime()       const { return (iTimeWindow == kInTimeBit) ? true : false ;}
		bool  IsLate()         const { return (iTimeWindow == kInTimeBit) ? true : false ;}
		int   GetNMuonStubs()  const { return iMuonStubs ;}
		bool  HasMuonStubs()   const { return (iMuonStubs>0) ? true : false ;}
				// here the magic number -100, is some time before the EM time gate starts
				// this wasy I can know EM timing info is available or not.
				// this will be some huge negative value in the later case.
				// fisrt 400pb-1 does not have EM timinng info.
		bool HasEmTime()      const { return (fEmTime > -100) ? true : false ;}
												
		//synonyms
		bool IsTightPhoton() const { 
			int id = GetTightPhotonId();
			if (id == 0) return true;
			else if (id > 0) return false;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Tight Photon ID is not set. returning 'false'.");
				return false;
			}
		}
		
		bool IsLoosePhoton() const { 
			int id = GetLoosePhotonId();
			if (id == 0) return true;
			else if (id > 0) return false;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Loose Photon ID is not set. returning 'false'.");
				return false;
			}
		}
		
		bool IsSidebandPhoton() const {
			int tid = GetTightPhotonId();
			int lid = GetLoosePhotonId();
			if (tid != 0 && lid == 0) return true;
			else return false;
		}



		bool IsTightElectron() const { 
			int id = GetTightElectronId();
			if (id == 0) return true;
			else if (id > 0) return false;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Tight Electron ID is not set. returning 'false'.");
				return false;
			}
		}
		
		bool IsStdTightElectron() const { 
			int id = GetStdTightElectronId();
			if (id == 0) return true;
			else if (id > 0) return false;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Standard Tight Electron ID is not set. returning 'false'.");
				return false;
			}
		}

		
		bool IsLooseElectron() const { 
			int id = GetLooseElectronId();
			if (id == 0) return true;
			else if (id > 0) return false;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Loose Electron ID is not set. returning 'false'.");
				return false;
			}
		}
		
		bool IsStdLooseElectron() const { 
			int id = GetStdLooseElectronId();
			if (id == 0) return true;
			else if (id > 0) return false;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Standard Loose Electron ID is not set. returning 'false'.");
				return false;
			}
		}

		
		bool IsConversion() const { 
			int id = GetConversionId();
			if (id == 0) return false;
			else if (id > 0) return true;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Conversion ID is not set. returning 'false'.");
				return false;
			}
		}
		
		bool IsBeamHalo() const { 
			int id = GetBeamHaloId();
			if (id == 0) return false;
			else if (id > 0) return true;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Beam Halo ID is not set. returning 'false'.");
				return false;
			}
		}
		
		bool IsCosmic() const {
			int id = GetCosmicId();
			if (id == 0) return false;
			else if (id > 0) return true;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Cosmic ID is not set. returning 'false'.");
				return false;
			}
		}

		bool IsPhoenix() const { 
			int id = GetPhoenixId();
			if (id == 0) return false;
			else if (id > 0) return true;
			else {
				if (iPrintLevel>5)
				StdOut(__FILE__,__LINE__,2,"Phoenix ID is not set. returning 'false'.");
				return false;
			}
		}

		int GetPrintLevel() const { return iPrintLevel;}

		//setters
		
		// needed these for JetFilterModule to remove electrons properly.
		// It uses the index from the Electron Block. not from photon Block.
		// -- sam, 05-13-2008
		void SetPhotonBlockIndex(const int ind_) { iPhoBlockIndex = ind_;}
		void SetElectronBlockIndex(const int ind_) { iEleBlockIndex = ind_;}
		// So I'll make this obsolete.
		void SetIndex(int ind_) { 
			std::cout << __FILE__ <<"::" << __LINE__ << ":: This function is obsolete. See header file for the new accessors." << std::endl;
		}
		
		void SetPhoton(TStnPhoton* pho_) { aPhoton = pho_ ;}
		void SetCorVec(TLorentzVector corvec_) { tlCorVec = corvec_ ;}
		void SetRawVec(TLorentzVector rawvec_) { tlRawVec = rawvec_ ;}
		void SetDetector(int det_) 	{ iDetector   = det_  ;}
		void SetDetEta(float eta_) 	{ fDetEta 	  = eta_  ;}
		void SetDetPhi(float phi_) 	{ fDetPhi 	  = phi_  ;}
		void SetEtRaw(float etr_) 		{ fEtRaw  	  = etr_  ;}			
		void SetECorr(float ec_) 		{ fECorr  	  = ec_   ;}		
		void SetEtCorr(float etc_) 	{ fEtCorr 	  = etc_  ;}		
		void SetXCes(float xces_) 		{ fXCes   	  = xces_ ;}			
		void SetZCes(float zces_) 		{ fZCes 		  = zces_ ;}
		void SetHadEm(float hadem_) 	{ fHadEm      = hadem_;}
		void SetChi2Mean(float cm_) 	{ fChi2Mean   = cm_   ;}
		void SetN3d(int n3d_) 			{ iN3d        = n3d_  ;}
		void SetIso4(float iso_) 		{ fIsoEtCorr  = iso_  ;}	
		void SetTrkPt(float tp_) 		{ fTrkPt      = tp_   ;}		
		void SetTrkPt2(float tp_) 		{ fTrkPt2     = tp_   ;}		
		void SetTrkIso(float ti_) 		{ fTrkIso     = ti_   ;}
		void SetCesWireE2(float c2_) 	{ fCesWireE2  = c2_   ;}
		void SetCesStripE2(float s2_) { fCesStripE2 = s2_   ;}
		void SetLoosePhotonId(int lid_) { iLoosePhotonId    = lid_  ;} 
		void SetTightPhotonId(int tid_) { iTightPhotonId    = tid_  ;}
		void SetLooseElectronId(int lid_) { iLooseElectronId = lid_  ;} 
		void SetStdLooseElectronId(int lid_) { iStdLooseElectronId = lid_ ;} 
		void SetTightElectronId(int eid_) { iTightElectronId = eid_  ;}
		void SetStdTightElectronId(int eid_) { iStdTightElectronId = eid_ ;}
		void SetBeamHaloId(int bhd_) 	{ iIsBeamHalo = bhd_  ;}	
		void SetCosmicId(int bc_) 		{ iIsCosmic   = bc_   ;}
		void SetConversionId(int bc_) { iIsConversion = bc_ ;}
		void SetPhoenixId(int bp_) 	{ iIsPhoenix  = bp_   ;}	
		void SetPMTSpikeId(int bp_) 	{ iIsPMTSpike = bp_   ;}	
		void SetHadTdcTime(float hd_) { fHadTdcTime = hd_   ;}
		void SetSinTheta(float st_) 	{ fSinTheta   = st_   ;}
		void SetEmTime(float time_)  	{ fEmTime     = time_ ;}	
		void SetEarlyTimeBit() 			{ iTimeWindow = kEarlyTimeBit 	;}
		void SetBeamHaloTimeBit() 		{ iTimeWindow = kBeamHaloTimeBit ;}
		void SetInTimeBit() 				{ iTimeWindow = kInTimeBit 		;}
		void SetLateTimeBit() 			{ iTimeWindow = kLateTimeBit 		;}
		void SetNMuonStubs(int mu_)    { iMuonStubs  = mu_   				;}
		void SetPrintLevel(const int pp_) { iPrintLevel = pp_;}


		
	private:
		int iIndex;     			//orginal order in the photon Block		-- no longer used. see comments above on 05-13-2008, sam
		int iPhoBlockIndex;    			//orginal order in the photon Block, replacing iIndex - 05-13-2008
		int iEleBlockIndex;    			//orginal order in the electron Block, replacing iIndex - 05-13-2008
		TStnPhoton*  aPhoton;				
		TLorentzVector tlCorVec;		//corrected 4-vec
		TLorentzVector tlRawVec;		//raw 4-vec
		int 	iDetector;	
		float fDetEta;
		float fDetPhi;
		float fEtRaw;				// Et()
		float fECorr;				// ECorr()
		float fEtCorr;				// ECorr() * SinTheta()
		float fXCes;			
		float fZCes;
		float fHadEm;
		float fChi2Mean;
		int   iN3d;
		float fIsoEtCorr;			// EIso4(2)
		float fTrkPt;				// Pt()
		float fTrkPt2;				// Pt2()  //2ns highest pt track
		float fTrkIso;				// SumPt4()
		float fCesWireE2;
		float fCesStripE2;
		int 	iLoosePhotonId;			// if passed =0, else > 0, not set < 0 
		int 	iTightPhotonId;			// if passed =0, else > 0, not set < 0
		int 	iLooseElectronId;			//    Photon-like electron ID , if passed =0, else > 0 , not set < 0 ,loose E/P compare to tight
		int 	iStdLooseElectronId;		// standard loose electron ID , if passed =0, else > 0 , not set < 0 ,loose E/P compare to tight
		int 	iTightElectronId;			//    Photon-like electron ID , if passed =0, else > 0 , not set < 0
		int 	iStdTightElectronId; 	// standard tight electron ID , if passed =0, else > 0 , not set < 0
		int 	iIsBeamHalo;				// -1 = not set, 0= no 1=yes
		int 	iIsCosmic;					// -1 - not set, 0= no 1=yes
		int 	iIsConversion;				// -1 - not set, 0= no 1=yes
		int 	iIsPhoenix;					// -1 - not set, 0= no 1=yes
		int 	iIsPMTSpike;				// -1 - not set, 0= no 1=yes
		float fHadTdcTime;				// timing of Had TDCs // not usefull. exist only
												// if there is really large shower extending to HAD cal.
		float fSinTheta;
		float fEmTime;       	// EM timing. 
		int iTimeWindow;     	// in which time window is this photon in.
										// see elog#250 and 2nd log book pp.2
										// -1 = not set, 0= not defined
		int iMuonStubs;			// Trkless muon stubs in the 30 degree cone

		int iPrintLevel;			// larger means more print statements
};

#endif
