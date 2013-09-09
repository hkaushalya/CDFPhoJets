#include "samantha/Pho/ZIdEffiency.hh"
#include "samantha/Pho/TagConvElectrons.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TGenParticle.hh"
#include "TLorentzVector.h"
#include <iostream>

ClassImp(ZIdEffiency)

//_____________________________________________________________________________
ZIdEffiency::ZIdEffiency(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am ZIdEffiency module" << std::endl;
}

//_____________________________________________________________________________
ZIdEffiency::~ZIdEffiency() {
}

//_____________________________________________________________________________
void ZIdEffiency::SaveHistograms() {
}

//_____________________________________________________________________________
void ZIdEffiency::BookHistograms()
{
	DeleteHistograms();
	BookZHistograms();
}


//_____________________________________________________________________________
int ZIdEffiency::BeginJob()
{
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.zideff_hepgZs = 0;				//__ hepg Z events
	counter.zideff_hepgZees = 0;			//__ Z that decayed to electron+ nutrino
	counter.zideff_passGoodrun = 0;		//__ from the events with a Z, number of events pass good run
	counter.zideff_passVertex = 0;		//__ then pass 1 good vertex
	counter.zideff_passZv = 0;
	counter.zideff_passZv = 0;				//__ then pass Zv cut
	counter.zideff_detZs = 0;


	return 0;
}

//_____________________________________________________________________________
int ZIdEffiency::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
  	TagConvElectrons* tagConv = (TagConvElectrons*) ((TStnAna*) GetAna()->GetModule("TagConvElectrons"));
	if (tagConv == NULL) {
		std::cout << "WARNING! Cannot find Conversion Tagging Module required!" <<std::endl;
	}
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int ZIdEffiency::Event(int ientry)
{
	SetPassed(1);
	counter.evtsRunOver++;
	if (GetPassed()) counter.evtsPassModule++;


	if (!qMc) {
		StdOut("This is not a MC sample. nothing is done!.",GetName());
		return 0;
	}
  
  	InitSuperPhotons* initphomod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initphomod == NULL) {
		std::cout << "Cannot find InitSuperPhoton module!" <<std::endl;
		return 0;
	}
	
	
	FillDataBlocks(ientry);
	
	if (!FoundHepgZee()) return 0;

  	int NsuperPho = initphomod->GetSuperPhoSize();
	int ne = 0;
	for (int i=0; i < NsuperPho ;i++ ) {
  		if ( !(initphomod->GetSuperPhoton(i)->IsConversion()) && 
				(initphomod->GetSuperPhoton(i)->IsTightElectron()) ) {
			ne++;
		}
	}
		
	if (ne == 2) counter.zideff_detZs++;
	return 0;

} // Event

//--------------------------------------------------------------	
// measure eff of CEM PHOTON LIKE ELE ID CUTS for Zee 
/*
HEPG LEVEL
1. find hepg Z
2. apply Z window cut to remove Drell-Yan
3. find its decay e's, keep their vectors
DETECTOR LEVEL
4. apply goodrun ,vertex,Zv cut
5. look for 2 photons passing ID CUTs
*/
//--------------------------------------------------------------	
bool ZIdEffiency::FoundHepgZee()
{

	int Nparticles = fGenpBlock->NParticles();
	//int mom_id = 0;
	TLorentzVector helevec(0,0,0,0), pv(0,0,0,0);	// temp ele vec and production vertex vector for fid cut
	std::vector<TLorentzVector> elevec;					// two electron 4-vectors
	std::vector<TLorentzVector> elepv;					// two electron productions vectors
	int ndau = 0;	// number of daughter from Zee = 2

	TGenParticle *par, *mom, *dau;
	
	for (int i = 0 ; i < Nparticles ; i++) {	//_____ not much diff from running over all. my match is mostly within the fist 10 objects
		par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
			mom = fGenpBlock->Particle(im);

			if (mom != 0) {
				int par_id   = par->GetPdgCode();
				int par_stat = par->GetStatusCode();
				if (abs(par_id) == 23 && (par_stat == 2)) { 	//found a Z 
					TLorentzVector zvec;
					par->Momentum(zvec);
					ZPlot.hepgZpt->Fill(zvec.Pt());
					ZPlot.hepgZmass_b4masscut->Fill(zvec.M());
					
					if (zvec.M() < 66. || zvec.M() > 116.) return false;  // must apply this to get the correct eff. 
					
					ZPlot.hepgZmass_a4masscut->Fill(zvec.M());
					
					counter.zideff_hepgZs++;
					int p1 = par->GetFirstDaughter();
					int pl = par->GetLastDaughter();
					
					for (int i= p1 ; i <= pl || ndau == 2; i++) {
						if (ndau != 2) {
							dau = fGenpBlock->Particle(i);
							int dau_id = dau->GetPdgCode();
							int dau_stat = dau->GetStatusCode();
							
							TLorentzVector temp;
							dau->Momentum(temp);
							
							if ( (abs(dau_id) == 11) && (dau_stat == 1)) {
								ndau++;
								dau->Momentum(helevec);
								dau->ProductionVertex(pv);
								elevec.push_back(helevec);	
								elepv.push_back(pv);	
								
							} //if (electron)
							
						} else {
							break;
						}// if (ndau)
						
					} //for
					
					break;	//Z is found, but did not decay in to e+e-, so quit
				} 
			} //if Z
		} // if
	} //for
	
	if (ndau == 2) {
		counter.zideff_hepgZees++;
		return true;
	} else return false;

} // FindHepgZee

//_____________________________________________________________________________
void ZIdEffiency::BookZHistograms()
{
	ZPlot.hepgZpt  	= new TH1F("hepgZpt","HEPG Z pt",15,0,15);
	ZPlot.hepgZmass_b4masscut  	= new TH1F("hepgZmass_b4masscut","HEPG Z mass",60,0,120);
	ZPlot.hepgZmass_a4masscut  	= new TH1F("hepgZmass_a4masscut","HEPG Z mass",60,0,120);
	ZPlot.detZmass  	= new TH1F("detZmass","Z mass at Detector",60,0,120);

	TFolder* new_folder = GetHistFolder(this, "Zplots","Z histos");
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(ZPlot.hepgZpt);
	new_folder->Add(ZPlot.hepgZmass_b4masscut);
	new_folder->Add(ZPlot.hepgZmass_a4masscut);
	new_folder->Add(ZPlot.detZmass);
}

/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void ZIdEffiency::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void ZIdEffiency::FillDataBlocks(int ientry)
{
	fGenpBlock->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int ZIdEffiency::EndJob() {

	printf("----- end job: ---- %s\n",GetName());
	if (qMc)	std::cout << "RUN IS ON MC SAMPLE" << std::endl;
	else	std::cout << "RUN IS ON DATA SAMPLE" << std::endl;
	std::cout << "Events Run Over ------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "Events Pass this module ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "Hepg Z events found -------- = " << counter.zideff_hepgZs << std::endl;
	std::cout << "Hepg Zee events found ------ = " << counter.zideff_hepgZees<< std::endl;
	std::cout << "Detector Z events found ---- = " << counter.zideff_detZs << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
