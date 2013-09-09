#ifndef HALOTYPES_HH
#define HALOTYPES_HH

class HaloTypes
{
	public:
		HaloTypes(int iSw = 0, int iNh = 0);
		bool IsType(const int) const;
		void SetSeedWedge(int isw) { iSeedWedge = isw; }
		void SetNhad(int inh) { iNhad = inh; }
		int GetSeedWedge() const { return iSeedWedge; }
		int GetNhad() const { return iNhad; }
		static const int iNHALOTYPES = 6;  // total number of halo types
		
	
	private:
		int iSeedWedge;
		int iNhad;			// east+west

};
#endif
