#include "PhoFlat/HaloTypes.hh"
/********************************************************************
 * $Id: HaloTypes.cc,v 1.2 2010/02/13 03:32:09 samantha Exp $
 *
 * $Log: HaloTypes.cc,v $
 * Revision 1.2  2010/02/13 03:32:09  samantha
 * MODIFIED:defined the number of halo types as a static constant of the class.
 * And added automatic CVS log info to the file.
 *
 * *****************************************************************/

const int HaloTypes::iNHALOTYPES;  // total number of halo types

//-------------------------------------------------------------------
HaloTypes::HaloTypes(int iSw, int iNh) 
{
	//Nhad = eastNhad + westNhad
	iSeedWedge = iSw;
	iNhad = iNh;
}

//-------------------------------------------------------------------
bool HaloTypes::IsType(const int iType) const
{
	bool bYes = false;
	
	switch (iType) {
		case 0:{
				if (iSeedWedge > 8 || iNhad > 2) bYes = true;
				break;
			}
		case 1: {
				if (iSeedWedge > 4 && iNhad > 1) bYes = true;
				break;
		}
		case 2: {
				if (iSeedWedge > 4 && iNhad > 2) bYes = true;
				break;
		}
		case 3: {
				if (iSeedWedge > 7 && iNhad > 2) bYes = true;
				break;
		}
		case 4: {
				if (iSeedWedge > 8 && iNhad > 2) bYes = true;
				break;
		}
		case 5: {
				if (iSeedWedge > 8 && iNhad > 3) bYes = true;
				break;
		}
	}

	return bYes;

}
