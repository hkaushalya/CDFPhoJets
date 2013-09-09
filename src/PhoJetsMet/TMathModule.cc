#include "samantha/PhoJetsMet/TMathModule.hh"
#include <TMath.h>

double const kPI= TMath::Pi();
double const kTWOPI = 2.*kPI;

//==========================================================
TMathModule::TMathModule()
{
	Init();
}
//==========================================================
TMathModule::~TMathModule()
{
}

//==========================================================
void TMathModule::Init()
{
	phi_1  = 999.99;
	phi_2  = 999.99;
	eta_1  = 999.99;
	eta_2  = 999.99;
	DelPhi = 999.99;
	DelEta = 999.99;
	DelR   = 999.99;
	margin = 999.99;
}

//==========================================================
void TMathModule::SetPhi(double p1, double p2)
{
	phi_1 = p1;
	phi_2 = p2;
}

//==========================================================
void TMathModule::SetEta(double e1, double e2)
{
	eta_1 = e1;
	eta_2 = e2;
}
//==========================================================
void TMathModule::SetMargin(double m)
{
	margin = m;
}
//==========================================================
void TMathModule::SetMEtaPhi(double p1,double e1, double p2, double e2, double m)
{
	SetEta(e1,e2);
	SetPhi(p1,p2);
	SetMargin(m);
}
	
//==========================================================
double TMathModule::GetDelPhi(double p1, double p2)
{
	SetPhi(p1,p2);
	SetDelPhi();
	return GetDelPhi();
/*
	if ( (phi_1 > kPI && phi_2 > kPI) || (phi_1 < kPI && phi_2 < kPI)) {
		DelPhi =  fabs(phi_1 - phi_2);
	} else {
		if ( fabs(phi_1 - phi_2) < kPI) {
			DelPhi = fabs(phi_1 - phi_2);
		} else {
			if (phi_1 > kPI) {
				DelPhi = fabs(kTWOPI - phi_1 + phi_2);
			} else {
				DelPhi = fabs(kTWOPI - phi_2 + phi_1);
			}
		}
	}
*/

}

//==========================================================
void TMathModule::SetDelPhi() 
{
	DelPhi = phi_1 - phi_2;

	if (DelPhi >kPI) {
		DelPhi -= kTWOPI;
	}
	if (DelPhi < -kPI) {
		DelPhi += kTWOPI;
	}

	DelPhi = fabs(DelPhi);
}

//==========================================================
double TMathModule::GetDelEta(double e1, double e2)
{
	SetEta(e1,e2);
	SetDelEta();
	return GetDelEta();
}

//==========================================================
void TMathModule::SetDelEta()
{
	DelEta = fabs(eta_1 - eta_2);
}

//==========================================================
void TMathModule::SetDelR()
{
	DelR = sqrt((pow(DelPhi,2)+pow(DelEta,2)));
}
//==========================================================
double TMathModule::GetDelR(double delPhi, double delEta)
{
	DelR = sqrt((pow(delPhi,2)+pow(delEta,2)));
	return DelR;
}

//==========================================================
double TMathModule::GetDelR(double phi_1, double eta_1, double phi_2, double eta_2)
{
	DelPhi = GetDelPhi(phi_1,phi_2);
	DelEta = GetDelEta(eta_1, eta_2);

	return GetDelR(DelPhi,DelEta);
}


//==========================================================
double TMathModule::GetDelR(TStnPhoton* pho, TStnJet* jet)
{
	phi_1 = pho->Phi();
	phi_2 = jet->Phi();
	eta_1 = pho->DetEta();
	eta_2 = jet->DetEta();
	DelR = GetDelR(phi_1,eta_1,phi_2,eta_2);
	return DelR;
}
//=========================================================
double TMathModule::GetDelPhi(TStnJet* j1, TStnJet* j2)
{
	phi_1 = j1->Phi();
	phi_2 = j2->Phi();
	DelPhi = GetDelPhi(phi_1,phi_2);
	return DelPhi;
}
//=========================================================
double TMathModule::GetDelEta(TStnJet* j1, TStnJet* j2)
{
	eta_1 = j1->DetEta();
	eta_2 = j2->DetEta();
	DelEta = GetDelEta(eta_1,eta_2);
	return DelEta;
}
//==========================================================
double TMathModule::GetDelR(TStnJet* j1, TStnJet* j2)
{
	phi_1 = j1->Phi();
	phi_2 = j2->Phi();
	eta_1 = j1->DetEta();
	eta_2 = j2->DetEta();
	DelR = GetDelR(phi_1,eta_1,phi_2,eta_2);
	return DelR;
}
//==========================================================
bool TMathModule::MatchEtaPhi(TStnPhoton* pho, TStnJet* jet, double m)
{
	DelR = GetDelR(pho,jet);
	return (DelR < m) ? true : false;
	
}
//==========================================================
bool TMathModule::MatchEtaPhi(double p1, double e1,
						double p2, double e2, double m)		//match the
{
	double DelR = GetDelR(p1,e1,p2,e2);
	return (DelR < m) ? true : false;
}

//==========================================================
void TMathModule::Print()
{
	printf("/*---------- TMathModule Summary ------------- */\n");
	printf(" phi_1,eta_1   =  %6.3f, %6.3f\n", phi_1 , eta_1); 
	printf(" phi_2,eta_2   =  %6.3f, %6.3f\n", phi_2 , eta_2); 
	printf(" DelPhi,DelEta =  %6.3f, %6.3f\n", DelPhi, DelEta); 
	printf(" DelR          =  %6.3f\n", DelR);
	printf("/*--------- END TMathModule Summary ---------- */\n");
}

//==========================================================
void TMathModule::binary(int number)
{
	int remainder;
	if(number < 0) {
		std::cout << "TMathModule:ERROR: Number given is Negative!\n";
		return;
	}
	if (number <= 1 ) {
		std::cout << number;
		return;
	}
	remainder = number%2;
	binary(number >> 1);    
	std::cout << remainder;
}
