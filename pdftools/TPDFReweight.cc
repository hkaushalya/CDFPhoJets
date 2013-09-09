
// TPDFReweight.cc contains the routines for the module which computes the reweighting 
//                 due to the PDFs (following Stephen Miller's method).
// Written by O. Gonzalez (21-V-2004)
// *************************************************************
/* Minor modifications for pho+jets analysis by sam. 
 *
 * $Id: TPDFReweight.cc,v 1.3 2009/11/19 22:37:16 samantha Exp $
 * $Log: TPDFReweight.cc,v $
 * Revision 1.3  2009/11/19 22:37:16  samantha
 * ADDED: Status of the GetConsiderAlphasEffect() to the End Job summary.
 *
 * Revision 1.2  2009/11/07 17:28:21  samantha
 * Removed the dependence on FreeFunctions module as it seems better to keep this
 * module independent. Added a smale ToStr() method to convert int/float/double to
 * strings. Added a more to end job summary output to dump the used settings.
 *
 * Revision 1.1  2009/10/23 16:48:39  samantha
 * Slight modifications to original to fit my needs. I needed store the weights
 * for later use and this does. Added a end job summary and commented out some
 * stuff I do not need in the endjob.
 *
 *
 */

#include "TF1.h"
#include <iostream>
#include <sstream>
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnLinkBlock.hh"
#include "TPDFReweight.hh"

//using namespace std;

extern "C" void pdfload_ (int *level,int *group, int *pdf, int *ierror, int *ioverload);

extern "C" void pdfprint_ (int *i);

extern "C" void getxpdf_ (int *level,int *ierror,double *x,double *q,double *xpdf);

extern "C" double getalphas_ (int *level,int *usefunc,double *q);

//_____________________________________________________________________________
TPDFReweight::TPDFReweight(const char* name, const char* title):
	TStnModule(name,title),
	bNoSummary(true)
{
	_usingMCgen   = "nonewasset";

	_usingScale   = -99;

	_pdfref_grp   = 4;   // Default is CTEQ5L
	_pdfref_set   = 46;

	_usingPDFset      = "CTEQ6M";
	_nalternativepdfs = 41;

	_hasstandarderrors = 1;

	_pdfalt_grp[0] = 4;
	_pdfalt_set[0] = 1;

	_consideralphaseffect = 1;

	_iordalphas = -1;

	for (int i=0;i<50;++i)
	{
		tot_weight[i]=0;
		tot_weightsqr[i]=0;
	}

	_indexHistogram = 0;


	_nEvents=0;

	/*
		_fail=0;
		_fail_zcut=0;           // |z|>Z_max_cut
		_pass=0;
		ZVertexVClass=12; // this means larger equal 12
		Z_max_cut=60;
		*/
}

//_____________________________________________________________________________
TPDFReweight::~TPDFReweight() {
}
//-____________________________________________________________________________
void TPDFReweight::SetReweightedwithPDF(std::string pdfname) {

	// The following sets are allowed

	if (pdfname == "CTEQ6M" )   // It is CTEQ6-M and the 40 errors
	{
		_usingPDFset       = "CTEQ6M";
		_nalternativepdfs  = 41;
		_hasstandarderrors = 1;

		_pdfalt_grp[0]=4;
		_pdfalt_set[0]=1;

	} else if (pdfname == "MRST98LO" )   // It is MRST98LO with 4 alternatives
	{
		_usingPDFset       = "MRST98LO";
		_nalternativepdfs  = 5;
		_hasstandarderrors = 0;

		for (int i=0;i<5;++i) _pdfalt_grp[i] = 3;

		_pdfalt_set[0]=72;
		_pdfalt_set[1]=73;
		_pdfalt_set[2]=74;
		_pdfalt_set[3]=75;
		_pdfalt_set[4]=76;

	} else if (pdfname == "CTEQ5L" )   // It is CTEQ5L (intended for checks)
	{
		_usingPDFset       = "CTEQ5L";
		_nalternativepdfs  = 1;
		_hasstandarderrors = 1;

		_pdfalt_grp[0] = 4;
		_pdfalt_set[0] = 46;

	} else if (pdfname == "MRST2001" )   // It is MRST2001C with 30 errors
	{
		_usingPDFset="MRST2001";
		_nalternativepdfs=31;
		_hasstandarderrors=1;

		_pdfalt_grp[0]=3;
		_pdfalt_set[0]=6;

	} else {
		std::cout<<"WARNING!!! Requested PDF for reweighting not valid "<<pdfname<<std::endl;
		_usingMCgen="nonewasset";
	}

}

//_____________________________________________________________________________
void TPDFReweight::BookHistograms() {
	char name [200];
	char title[200];
	// book histograms
	//  GetListOfHistograms()->Delete();
	DeleteHistograms();
	//-----
	// Vertex staff here
	//-----
	/*
		sprintf(name, "%s_n_Zvertex",GetName());
		sprintf(title,"%s: Number of Zverticies (qual>=%d)",GetName(),ZVertexVClass);
		fHist.fZVertexN = new TH1F(name, title,10 ,0,10);
		AddHistogram(fHist.fZVertexN);

		sprintf(name, "%s_z_Zvertex",GetName());
		sprintf(title,"%s: Z of Zverticies (qual>=%d) highest sumpt",GetName(),ZVertexVClass);
		fHist.fZVertexZ = new TH1F(name, title,160 ,-80,80);
		AddHistogram(fHist.fZVertexZ);
		*/
	sprintf(name, "%s_etgenparticles",GetName());
	sprintf(title,"%s: Et of the generated particles",GetName());
	fHist.fEtGenParticles = new TH1F(name, title,25 ,0.0,500.0);
	AddHistogram(fHist.fEtGenParticles);
	fHist.fEtGenParticles->Sumw2();

	sprintf(name, "%s_ptgenparticles",GetName());
	sprintf(title,"%s: pt of the generated particles",GetName());
	fHist.fPtGenParticles = new TH1F(name, title,20 ,0.0,300.0);
	AddHistogram(fHist.fPtGenParticles);
	fHist.fPtGenParticles->Sumw2();

	sprintf(name, "%s_etagenparticles",GetName());
	sprintf(title,"%s: Eta of the generated particles",GetName());
	fHist.fEtaGenParticles = new TH1F(name, title,14 ,-3.5,3.5);
	AddHistogram(fHist.fEtaGenParticles);
	fHist.fEtaGenParticles->Sumw2();

	sprintf(name, "%s_etdecparticles",GetName());
	sprintf(title,"%s: Et of the decayed particles",GetName());
	fHist.fEtDecParticles = new TH1F(name, title,25 ,0.0,500.0);
	AddHistogram(fHist.fEtDecParticles);
	fHist.fEtDecParticles->Sumw2();

	sprintf(name, "%s_ptdecparticles",GetName());
	sprintf(title,"%s: pt of the generated particles",GetName());
	fHist.fPtDecParticles = new TH1F(name, title,20 ,0.0,300.0);
	AddHistogram(fHist.fPtDecParticles);
	fHist.fPtDecParticles->Sumw2();

	sprintf(name, "%s_etadecparticles",GetName());
	sprintf(title,"%s: Eta of the decayed particles",GetName());
	fHist.fEtaDecParticles = new TH1F(name, title,14 ,-3.5,3.5);
	AddHistogram(fHist.fEtaDecParticles);
	fHist.fEtaDecParticles->Sumw2();

	sprintf(name, "%s_etgenparticles_orig",GetName());
	sprintf(title,"%s: Et of the generated particles (original PDF)",GetName());
	fHist.fEtGenParticles_orig = new TH1F(name, title,25 ,0.0,500.0);
	AddHistogram(fHist.fEtGenParticles_orig);

	sprintf(name, "%s_ptgenparticles_orig",GetName());
	sprintf(title,"%s: pt of the generated particles (original PDF)",GetName());
	fHist.fPtGenParticles_orig = new TH1F(name, title,20 ,0.0,300.0);
	AddHistogram(fHist.fPtGenParticles_orig);

	sprintf(name, "%s_etagenparticles_orig",GetName());
	sprintf(title,"%s: Eta of the generated particles (original PDF)",GetName());
	fHist.fEtaGenParticles_orig = new TH1F(name, title,14 ,-3.5,3.5);
	AddHistogram(fHist.fEtaGenParticles_orig);

	sprintf(name, "%s_etdecparticles_orig",GetName());
	sprintf(title,"%s: Et of the decayed particles (original PDF)",GetName());
	fHist.fEtDecParticles_orig = new TH1F(name, title,25 ,0.0,500.0);
	AddHistogram(fHist.fEtDecParticles_orig);

	sprintf(name, "%s_ptdecparticles_orig",GetName());
	sprintf(title,"%s: pt of the generated particles (original PDF)",GetName());
	fHist.fPtDecParticles_orig = new TH1F(name, title,20 ,0.0,300.0);
	AddHistogram(fHist.fPtDecParticles_orig);

	sprintf(name, "%s_etadecparticles_orig",GetName());
	sprintf(title,"%s: Eta of the decayed particles (original PDF)",GetName());
	fHist.fEtaDecParticles_orig = new TH1F(name, title,14 ,-3.5,3.5);
	AddHistogram(fHist.fEtaDecParticles_orig);
}


//_____________________________________________________________________________
int TPDFReweight::BeginJob() {

	// Checking the values of the parameters

	if ( !( 
				_usingMCgen == "nonewasset" ||
				_usingMCgen == "Isajet"     ||
				_usingMCgen == "Pythia"  
			)) 
	{
		std::cout<<"WARNING!!! MC generator set in TPDFReweight is not correct"<<std::endl;
		_usingMCgen="nonewasset";
	}

	// Setting the PDFTools values for the reference PDF
	int level     = 0;
	int ioverload = 0;
	int ierr      = 0;
	pdfload_(&level,&_pdfref_grp,&_pdfref_set,&ierr,&ioverload);

	// The alternatives PDFs... 

	level = 1;
	if (_hasstandarderrors && _nalternativepdfs>1)
	{
		ierr = 1;
	}

	pdfload_(&level,&_pdfalt_grp[0],&_pdfalt_set[0],&ierr,&ioverload);

	// register the data block

	fGenpBlock = (TGenpBlock*) RegisterDataBlock("GenpBlock","TGenpBlock");
	//  RegisterDataBlock("GenpBlock","TStnGenpBlock",&fGenpBlock);

	if (!fGenpBlock) {
		printf("TPDFReweight >>> branch *** %s *** doesn't exist \n","fGenpBlock");
		fEnabled = 0;
	}

	fHeaderBlock = (TStnHeaderBlock*)  RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	// book histograms
	BookHistograms();

	iPassed = 0;
	iProcess = 0;

	
	return 0;
}


//_____________________________________________________________________________
int TPDFReweight::BeginRun() {
	return 0;
}

//_____________________________________________________________________________
int TPDFReweight::Event(int ientry) {

	SetPassed(1);  // This module does not filter anything.

	++iProcess;
	_nEvents++;

	vPDFweights.clear();

	// Computing Q and Bjorken's x's of the event and identifying the partons

	double x1=-1;
	double x2=-1;
	double q=-1;

	int part1=-10;
	int part2=-10;

	TLorentzVector pgen1,pgen2,pdec1,pdec2;
	fGenpBlock->GetEntry(ientry);

	bool thereisdecayedparticles=0;
	int numbergeneratedparticles=0;

	//  std::cout<<"Value: "<<_usingMCgen<<std::endl; 

	if (_usingMCgen == "Isajet") {   // Code for getting kinematic variables in ISAJET 
		TGenParticle* p1 = fGenpBlock->Particle(0);
		TGenParticle* p2 = fGenpBlock->Particle(1);

		x1 = (p1->Pz()+p2->Pz())/980.0;   // X_F = x1-x2
		q = pow(p1->Energy()+p2->Energy(),2)-pow(p1->Pz()+p2->Pz(),2);   // s-hat = s*x1*x2
		q = q/(1960*1960);  // this is x1*x2 = s-hat/s

		x2 = (-x1+sqrt(x1*x1 + 4*q))/2;

		x1 = (x1+sqrt(x1*x1 + 4*q))/2;

		// The following expression is a complete shit because is the PDG mass... for the gluino 500.
		//    q = p1->GetMass()+p2->GetMass();   // In ISAJET Q es la suma de las masas
		q = sqrt(p1->Energy()*p1->Energy()-p1->Px()*p1->Px()-p1->Py()*p1->Py()-p1->Pz()*p1->Pz()) +
			sqrt(p2->Energy()*p2->Energy()-p2->Px()*p2->Px()-p2->Py()*p2->Py()-p2->Pz()*p2->Pz());

		if (q<=0 || x1<=0 || x2<=0) {
			std::cout<<"WARNING!!! Problems identifying kinematics in Isajet "<<x1<<" "<<x2<<" "<<q<<std::endl;
		}

		// Apart from identifying kinematics, one should identify the partos to evaluate the PDFs

		// In ISAJET these are the first partons after the two produced particles

		p1 = fGenpBlock->Particle(2);   // Isajet follows the PDFLib convention.
		part1 = p1->GetPdgCode();
		if (part1==21) part1=0;                  // Except for the gluon
		p2 = fGenpBlock->Particle(3);
		part2 = p2->GetPdgCode();
		if (part2==21) part2=0;


		// Getting the information on the produced particles in the Hard interaction:

		p1 = fGenpBlock->Particle(0);  // Generated particles (well, they could be 0/1 or 4/5 though)
		p2 = fGenpBlock->Particle(1);

		numbergeneratedparticles=2;

		pgen1.SetPxPyPzE(p1->Px(),p1->Py(),p1->Pz(),p1->Energy());
		pgen2.SetPxPyPzE(p2->Px(),p2->Py(),p2->Pz(),p2->Energy());

		Int_t j=0;
		thereisdecayedparticles=1;
		for (int i=6;i<fGenpBlock->NParticles();++i) {
			TGenParticle* paux = fGenpBlock->Particle(i);

			if (paux->GetNDaughters()>1 && 
					(paux->GetPdgCode() == p1->GetPdgCode() ||   // Decayed particles
					 paux->GetPdgCode() == p2->GetPdgCode() )
				) {
				if (j==0) 
					pdec1.SetPxPyPzE(paux->Px(),paux->Py(),paux->Pz(),paux->Energy());
				else if (j==1) 
					pdec2.SetPxPyPzE(paux->Px(),paux->Py(),paux->Pz(),paux->Energy());
				else
					std::cout<<"WARNING!!! Problems identifying final particles for ISAJET "<<p1->GetPdgCode()<<" "<<p2->GetPdgCode()<<std::endl;

				++j;
			}
		}

	} else if (_usingMCgen == "Pythia")   // Code for getting kinematic variables in PYTHIA 
	{

		// For pythia the calculation is a bit trickier.

		TGenParticle* prot = fGenpBlock->Particle(0); //proton PDG 2212
		TGenParticle* pbar = fGenpBlock->Particle(1); //antiproton PDG -2212

		TGenParticle* p1 = fGenpBlock->Particle(4);
		TGenParticle* p2 = fGenpBlock->Particle(5);

		// Getting the partons and information about the total momentum.
		part1 = p1->GetPdgCode();
		part2 = p2->GetPdgCode();
		if (p1->GetMother(0) != 2 || p2->GetMother(0) != 3) {
			std::cout<<"WARNING!!! TPDFReweight has problems to find PYTHIA partons "<<part1<<" "<<part2<<" "<<p1->GetMother(0)<<" "<<p2->GetMother(0)<<std::endl;
		}
		if (part1==21) part1=0;
		if (part2==21) part2=0;

		TLorentzVector totcms(p1->Px()+p2->Px(),p1->Py()+p2->Py(),p1->Pz()+p2->Pz(),p1->Energy()+p2->Energy());

		// The x's are computed using the s-hat variable which is the mass of the totcms... 

		float_t energ = totcms.M()/(prot->Energy()+pbar->Energy());   // ^s = s*x1*x2   (energ is sqrt(s^/s))

		// And a suggested variable by Torbjorn Sjostrand:

		float eta =  sqrt( (p1->Energy()+p2->Energy()+p1->Pz()+p2->Pz())/(p1->Energy()+p2->Energy()-p1->Pz()-p2->Pz())); // This one is sqrt(x1/x2)

		// Thus:
		x1 = energ * eta;
		x2 = energ/eta;

		// The standard Q in PYTHIA (no-resonances) is a bit complicated because we need pt-hat which is
		// not a Lorentz invariant. For this reason, it is only computed if the scale is not set by hand (to save time

		if (_usingScale<0) {

			// pthat is computed in the CMS: boosting there one initial and one final parton

			float alpha = -atan2(totcms.Py(),totcms.Px());
			float beta = atan2(totcms.Pt(),totcms.Pz());
			float zeta = totcms.P()/totcms.Energy();
			float gamma = 1/sqrt(1-zeta*zeta);

			// For sure this could be simplified... 

			float pcm1[3] = {p1->Px()*cos(alpha)*cos(beta) - p1->Py()*cos(beta)*sin(alpha) - p1->Pz()*sin(beta),
				p1->Px()*sin(alpha) + p1->Py()*cos(alpha),
				gamma*(p1->Px()*cos(alpha)*sin(beta) - p1->Py()*sin(beta)*sin(alpha) + p1->Pz()*cos(beta) - p1->Energy()*zeta)};

			TGenParticle* p3 = fGenpBlock->Particle(6);
			TGenParticle* p4 = fGenpBlock->Particle(7);

			float pcm2[3] = {p3->Px()*cos(alpha)*cos(beta) - p3->Py()*cos(beta)*sin(alpha) - p3->Pz()*sin(beta),
				p3->Px()*sin(alpha) + p3->Py()*cos(alpha),
				gamma*(p3->Px()*cos(alpha)*sin(beta) - p3->Py()*sin(beta)*sin(alpha) + p3->Pz()*cos(beta) - p3->Energy()*zeta)};

			float pthat = pcm1[0]*pcm2[0] + pcm1[1]*pcm2[1] + pcm1[2]*pcm2[2];
			pthat = pthat*pthat/(pcm1[0]*pcm1[0] + pcm1[1]*pcm1[1] + pcm1[2]*pcm1[2]);
			pthat = (pcm2[0]*pcm2[0] + pcm2[1]*pcm2[1] + pcm2[2]*pcm2[2] - pthat);  // It is pthat**2

			// Q in PYTHIA:

			q = sqrt(pthat + (p1->Energy()*p1->Energy()-p1->Px()*p1->Px()-p1->Py()*p1->Py()-p1->Pz()*p1->Pz() 
						+p2->Energy()*p2->Energy()-p2->Px()*p2->Px()-p2->Py()*p2->Py()-p2->Pz()*p2->Pz() 
						+p3->Energy()*p3->Energy()-p3->Px()*p3->Px()-p3->Py()*p3->Py()-p3->Pz()*p3->Pz() 
						+p4->Energy()*p4->Energy()-p4->Px()*p4->Px()-p4->Py()*p4->Py()-p4->Pz()*p4->Pz())/2);
			

		}

		// Getting the generated particles (no decayed particles)
		// IN PYTHIA THE GENERATED PARTICLES ARE 7 and so on... marked with first mother=0

		p1 = fGenpBlock->Particle(6);
		p2 = fGenpBlock->Particle(7);
		int mothidx=p1->GetMother(0);   // It is usually 4, but just in case...
		numbergeneratedparticles++;
		pgen1.SetPxPyPzE(p1->Px(),p1->Py(),p1->Pz(),p1->Energy());

		if (p2->GetMother(0)==mothidx) {
			numbergeneratedparticles++;
			pgen2.SetPxPyPzE(p2->Px(),p2->Py(),p2->Pz(),p2->Energy());
		}

	}

	// checks...

	if (part1<-6 || part1>6 || part2<-6 || part2>6) {
		std::cout<<"WARNING!!! Strange parton in the proton/antiproton: "<<part1<<" "<<part2<<std::endl;
		fGenpBlock->Print();
		q=-1;
	}

	// If Q is set by hand, that one is used

	if (_usingScale>0) q=_usingScale;

	//this is the final q value used in calcuating weights. So I'll use it in the next modules
	Q = q;	

	
	// with the x's and q we compute the weight of the event

	double weight=1;

	if (x1>=0 && x2>=0 && q>=0) {  // Only if no problems/or MC information available

		int level=0;
		int ierr=0;
		double xpdf[13];
		for (int i=0;i<13;++i) xpdf[i]=0;

		getxpdf_(&level,&ierr,&x1,&q,xpdf);
		weight = 1/xpdf[part1+6];

		getxpdf_(&level,&ierr,&x2,&q,xpdf);

		weight /= xpdf[-part2+6];   // antiproton

		if (_consideralphaseffect) weight /= pow(getalphas_ (&level,&_iordalphas,&q),2);

		// Now the alternative PDFs:

		level=1;
		int ioverload=-1;  // Overwriting silently
		for (int i=0;i<_nalternativepdfs;++i)
		{
			double weight2;
			if (_hasstandarderrors) {  // There is one PDF + errors (or just one PDF)
				getxpdf_(&level,&i,&x1,&q,xpdf);
				weight2=xpdf[part1+6];

				//	std::cout<<"uno "<<weight2<<std::endl;

				getxpdf_(&level,&i,&x2,&q,xpdf);
				weight2 *=xpdf[-part2+6];

				//	std::cout<<"dos "<<weight2<<" "<<xpdf[-part2+6]<<std::endl;

			} else {   // There are alternatives PDFs
				pdfload_(&level,&_pdfalt_grp[i],&_pdfalt_set[i],&ierr,&ioverload);
				getxpdf_(&level,&ierr,&x1,&q,xpdf);
				weight2=xpdf[part1+6];

				getxpdf_(&level,&ierr,&x2,&q,xpdf);
				weight2 *=xpdf[-part2+6];
			}
			
			if (_consideralphaseffect) weight2 *= pow(getalphas_ (&level,&_iordalphas,&q),2);

			weight2 *= weight;  // Final weight of the event

			vPDFweights.push_back(weight2); //storing this weight for other modules to access this weight - sam

			// For the _indexHistogram PDF we fill some histograms with the event inforamtion
			// We do the same without the weight.

			if (i==_indexHistogram) {
				fHist.fEtGenParticles->Fill(pgen1.E()*pgen1.Pt()/pgen1.P(),weight2);   // Generated particles
				fHist.fPtGenParticles->Fill(pgen1.Pt(),weight2);
				fHist.fEtaGenParticles->Fill(pgen1.PseudoRapidity(),weight2);

				fHist.fEtGenParticles_orig->Fill(pgen1.E()*pgen1.Pt()/pgen1.P());
				fHist.fPtGenParticles_orig->Fill(pgen1.Pt());
				fHist.fEtaGenParticles_orig->Fill(pgen1.PseudoRapidity());

				if (numbergeneratedparticles>1) {
					fHist.fEtGenParticles->Fill(pgen2.E()*pgen2.Pt()/pgen2.P(),weight2);
					fHist.fPtGenParticles->Fill(pgen2.Pt(),weight2);
					fHist.fEtaGenParticles->Fill(pgen2.PseudoRapidity(),weight2);

					fHist.fEtGenParticles_orig->Fill(pgen2.E()*pgen2.Pt()/pgen2.P());
					fHist.fPtGenParticles_orig->Fill(pgen2.Pt());
					fHist.fEtaGenParticles_orig->Fill(pgen2.PseudoRapidity());
				}

				if (thereisdecayedparticles) {
					fHist.fEtDecParticles->Fill(pdec1.E()*pdec1.Pt()/pdec1.P(),weight2);
					fHist.fPtDecParticles->Fill(pdec1.Pt(),weight2);
					fHist.fEtaDecParticles->Fill(pdec1.PseudoRapidity(),weight2);

					fHist.fEtDecParticles->Fill(pdec2.E()*pdec2.Pt()/pdec2.P(),weight2);
					fHist.fPtDecParticles->Fill(pdec2.Pt(),weight2);
					fHist.fEtaDecParticles->Fill(pdec2.PseudoRapidity(),weight2);  


					fHist.fEtDecParticles_orig->Fill(pdec1.E()*pdec1.Pt()/pdec1.P());
					fHist.fPtDecParticles_orig->Fill(pdec1.Pt());
					fHist.fEtaDecParticles_orig->Fill(pdec1.PseudoRapidity());

					fHist.fEtDecParticles_orig->Fill(pdec2.E()*pdec2.Pt()/pdec2.P());
					fHist.fEtDecParticles_orig->Fill(pdec2.Pt());
					fHist.fEtaDecParticles_orig->Fill(pdec2.PseudoRapidity());  
				}
			}
			// For all the alternatives we keep the total weight (and the square for error calculation).

			tot_weight[i] += weight2;
			tot_weightsqr[i] += weight2*weight2;
		
		} //for _nalternativepdfs loop

	}

	/*
	//  fTrackBlock    ->GetEntry(ientry);  
	fVertexBlock ->GetEntry(ientry);  // CR test if PrimVtx is ok
	fZVertex     ->GetEntry(ientry);
	//  fHeaderBlock ->GetEntry(ientry);

	if (fPrintLevel) {
	fZVertex->Print();
//    fVertexBlock ->Print();
GetHeaderBlock()->Print();
}

int nzvtx = fZVertex->NVertices();     // # vertices from ZVertexFinder
int nzvtxcount=0;  // # vertices from ZVertexFinder > qual
int finalvertexid;
float maxsumptZVertex=0;
//    float ZVertex; 
bool ZVertex_ok=false;

// now look for highest sumpt vertex that passes qualitt cuts
//       printf("nzvtx=%d\n",nzvtx);
//     printf("%d ",nzvtx);
for (int i = 0; i < nzvtx; i++)
{
TStnVertex* v = fZVertex->Vertex(i);
	//printf(" q=%d z=%2.2f\n",v->VClass(),v->Z());
	//  printf(" %d %2.2f ",v->VClass(),v->Z());
	if(v->VClass()>(ZVertexVClass-1))
	{ 
	ZVertex_ok=true;
	nzvtxcount++;
	if ( v->SumPt() > maxsumptZVertex)
	{
	maxsumptZVertex = v->SumPt();
	finalvertexid=i; 
	}
	}
	}
//  printf("\n");

////////////  Now test PRIMVTX  //////////////

int nvtx = fVertexBlock->NVertices(); // vertices from PRIMVTX
// printf("nvtx = %d \n",nvtx);

for (int i = 0; i < nvtx; i++)
{
TStnVertex* v = fVertexBlock->Vertex(i);
	//     printf("PRIMVTX      ntrks=%d sumpt=%f z=%f \n",v->NTracks(),v->SumPt(),v->Z());
	//    printf("%d %4.2f %2.2f \n",v->NTracks(),v->SumPt(),v->Z()); 
	}

//////////////////////////////////////////////

fHist.fZVertexN->Fill(nzvtxcount);

if(ZVertex_ok)
{
TStnVertex* v = fZVertex->Vertex(finalvertexid);
if(fabs(v->Z())>Z_max_cut)
{
SetPassed(0);    
_fail_zcut++;   
}
else
{
SetPassed(1);    
_pass++;
fHist.fZVertexZ->Fill(v->Z());
	//NTracksZVertex = v->NTracks();
	}  
	}
else
{
	SetPassed(0);  
	_fail++;
}
*/ 

	if (GetPassed()) ++iPassed;

return 0;			       
}

//_____________________________________________________________________________
int TPDFReweight::EndJob() {



	//  std::cout<<"No Track Events: "<<_NoTrk<<std::endl;
	/*      std::cout<<"Events pass: "<<_pass<<std::endl;    
			  std::cout<<"Events fail (q<"<<ZVertexVClass<<"): "<<_fail<<std::endl;    
			  std::cout<<"Events fail (|z|>"<<Z_max_cut <<"): "<<_fail_zcut<<std::endl;    
			  */
	if (! GetSummaryStat())
	{
		std::string sModName(GetName());
		std::string sMsg;
		int igrp=-1, ipdfset=-1;
		GetReferencePDF(igrp, ipdfset);
		sMsg  = "[PDF:00:]----- end job: ---- " + sModName + "\n";
		sMsg += "[PDF:01:] Events Processed ----------- = " + ToStr(iProcess) + "\n";
		sMsg += "[PDF:02:] Events Passed -------------- = " + ToStr(iPassed) + "\n";
		sMsg += "[PDF:03:] MC Generator --------------- = " + GetMCGenerator() + "\n";
		sMsg += "[PDF:04:] # of alternative PDFs ------ = " + ToStr(GetNalternativePDFs()) + "\n";
		sMsg += "[PDF:05:] Reference PDF(group, pdf set)= " + ToStr(igrp) + "," + ToStr(ipdfset) + " (4,46 is CTEQ5L defaults)\n";
		sMsg += "[PDF:06:] Reweighted to PDF ---------- = " + GetReweightedwithPDF() + "\n";
		sMsg += "[PDF:07:] MC Scale (q==) ------------- = " + ToStr(GetMCScale()) +  " (<0 == not set)\n";
		sMsg += "[PDF:08:] Considering alpha_s effects  = " + ToStr(GetConsiderAlphasEffect()) +  " (1==YES,0==NO)\n";
		sMsg += "---------------------------------------------------";
		std::cout << sMsg << std::endl;
	}
	
	//I DO NOT NEED THIS STUFF  -sam 10-8-2009
	/*
	printf("Reweighting: %s %i\n",_usingPDFset.c_str(),_nalternativepdfs);
	for (int i=0;i<_nalternativepdfs;++i) {
		double weight =  GetAvrWeight(i);
		double error = GetErrorOnAvrWeight(i);

		printf("- Average Weight:  %i %f +/- %f\n (%i)\n",i,weight,error,_nEvents);
	}
	*/

	// Saving the histograms is done automatically...  
	//fHist.Write();

	return 0;
}

//_____________________________________________________________________________
double TPDFReweight::GetErrorOnAvrWeight (const unsigned ipdf) const
// Returns the error on the average weight
{
	double weight =  GetAvrWeight(ipdf);

	double error = (GetTotWeightSqr(ipdf)/_nEvents - weight*weight)/_nEvents;// We want error on averagen no RMS of weights
	if (error>0) error = sqrt(error);
	else error = -sqrt(-error);

	return error;
}

/*
//_____________________________________________________________________________
void TPDFReweight::PlotHistograms(int run_number, int slide) {
char t[10];
sprintf(t,"%d",run_number);
//  PlotHistograms(t, slide);
}
//_____________________________________________________________________________
void TPDFReweight::PlotHistograms(char *run_number, int slide) {
// plot slides

char name[120], canvas_name[120], title[120], pad_name[120];

sprintf(name,"run_%s_%i",run_number,slide);

sprintf(title,"run %s slide %i ",run_number, slide);

sprintf(canvas_name,"%s_%s",GetName(),name);
sprintf(pad_name,"%s_p1",canvas_name);

if (slide == 1) {
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
//    TPostScript ps("l1ana.ps",-111);

TCanvas* c = NewSlide(name,title,2,2);
TPad* p1 = (TPad*) c->GetPrimitive(pad_name);

p1->cd(1); 
fHist.fZVertexZ->Draw();
p1->cd(2); 
fHist.fZVertexN->Draw();

gPad->Update();
}

}
*/
/////////////////////////////////////////////////////////////////////////

std::string TPDFReweight::ToStr(const double num)
{
	std::stringstream sStr;
	sStr << num;
	return sStr.str();
}
