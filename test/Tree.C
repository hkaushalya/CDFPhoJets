#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "../samantha/obj/Stuple.hh"
#include "../RootTools/CommonTools.hh"


void Tree()
{

	TFile f("TreeFile.root","RECREATE");
//	TTree *tree = new TTree ("PDFWeights","PDF Weights for Pho+Jets Events");
//		tree->Branch("EventNumber", iEvent,"EventNumber/I");
//		const int iNalt = iNaltenativePDFs;
//		tree->Branch("PDFWeights", &aPDFsyst.fPdfWeight,"PDFWeights[iNalt]/I");

	//TChain *ch = new TChain("PJevtsWithPDFw");
//	TChain *ch = new TChain("Ana/PhoJetsTemp/Hist/PJevtsWithPDFw/PhoJets");
	
	//TFile *f1 = new TFile("RESULTS/10222009_PDFsystFullMCDataset/11122009_ForAllPlotsWithPDFTree/PDFsyst1.root");
	//TFile *f1 = new TFile("RESULTS/10222009_PDFsystFullMCDataset/11122009_ForAllPlotsWithPDFTree/PDFsyst1.root");

//	ch->Add("~/RESULTS/10222009_PDFsystFullMCDataset/11122009_ForAllPlotsWithPDFTree/PDFsyst1.root");
//	ch->Add("~/RESULTS/10222009_PDFsystFullMCDataset/11122009_ForAllPlotsWithPDFTree/PDFsyst2.root");

	TTree *tree = 0;
	Stuple stuple;
	GetStupleTree("Stuple",tree,stuple);
	//std::cout << "entries =  " << ch->GetEntries() << std::endl;
	f.Write();
	f.Close();

}
