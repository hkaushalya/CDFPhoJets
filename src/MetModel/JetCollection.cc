#include "samantha/MetModel/JetCollection.hh"

JetCollection::JetCollection()
{
	ClearAll();
}

void JetCollection::ClearAll()
{
	Njet.iNjet = 0;				// all the jets
	Njet.iNjet_5 = 0;				// number if jets above Pt>5.0GeV
	Njet.iNjet_10 = 0;			// number if jets above Pt>10.0GeV
	Njet.iNjet_15 = 0;			// number if jets above Pt>15.0GeV
	Njet.iNjet_20 = 0;			// number if jets above Pt>20.0GeV
	Njet.iNjet_25 = 0;			// number if jets above Pt>25.0GeV
	Njet.iNjet_30 = 0;			// number if jets above Pt>30.0GeV
	Njet.iNjet_35 = 0;			// number if jets above Pt>35.0GeV

	JetCorrPara.iLevel = -1;
	JetCorrPara.iNvx   = -1;
	JetCorrPara.iJetCone = -1;
	JetCorrPara.iVersion = -1;
	JetCorrPara.iSys	   = -1;
	JetCorrPara.iRunNumber = -1;
	JetCorrPara.iImode   = -1;

	//vJets.clear();
	vfSmearFactors.clear();
}

float JetCollection::SumEtJets(const float fPt) const 
{
	float fSumEt = 0;
	/*for (int i=0; i < vJets.size(); ++i)
	{
		float fJetPt = vJets.at(i).tlMomentum.Pt();
		if (vJets.at(i).tlvRawMomentum.Pt() > fPt)
		{
			fSumEt += fJetPt;
		}
	}
*/
	return fSumEt;
}
