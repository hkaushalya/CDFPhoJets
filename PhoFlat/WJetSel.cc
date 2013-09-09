/*{{{*/
/**********************************************************
 * $Id: WJetSel.cc,v 1.2 2011/05/26 19:42:09 samantha Exp $
 * $Log: WJetSel.cc,v $
 * Revision 1.2  2011/05/26 19:42:09  samantha
 * ADDED: 1. Nvtx weights for different MC samples.
 * 2. Counting of events as cuts are applied.
 *
 * Revision 1.1  2010/08/23 03:06:54  samantha
 * W+Jet events are studied to derive additional corrections to EWK MC samples
 * to improve the MET predictions. W+jets are similar to g+jets. I use Pt and Eta
 * cuts similar to the g, on the W. I also apply a standard W transverse mass cut
 * to reduce the backgrounds.
 *
 **********************************************************
 *Samantha K. Hewamanage
 *********************************************************/
/*}}}*/

#include <iostream>
#include <sstream>
#include "PhoFlat/WJetSel.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "TBenchmark.h"
#include "TCanvas.h"
#include <iomanip>
#include "TROOT.h"
//for debugging
#include <exception>
#include <stdexcept> // out_of_range exception

void DEBUG(const std::string func, const double line, const std::string msg)
{
	std::cout << func << ":" << line << ": " << 	msg << std::endl;
}
 
//-------------------------------------------------------------------
WJetSel::WJetSel():
iProgressBy(50000),
sHistFileName("WJetSel.root"),
bMcFlag(true),
iUseNvtxWgts(0)
{
	fMinEleEt = 15.; // GeV  //do not lower this any more. See header file Set method for more details!
	fMaxEleEt = 1500.; // GeV
	fMaxEleEta = 1.1;
	iMinClass12Vtx = 1;
	iMaxClass12Vtx = 100;
	fMaxVtxz = 60.; //cm
	fMinMet = 0.0;	//minimum Missing transverse energy
	fMaxMet = 1500.0;	//maximum Missing transverse energy
	fMinWpt = 30.;
	fMinJetEt = 15.;
	fMaxJetEta = 3.;

}

//-------------------------------------------------------------
void WJetSel::Init(TChain *tree)
{
	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->CentralPhotons(1);
	readIn->Met(1);

	for (int i=0; i<4; i++) {
		vMcCount.push_back(0);
	}

	readIn->Init();

	//cannot add folder without creating a dir first!
	histFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = histFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

	vNvtxWgts.clear();
	if (UseNvtxWeights())
	{

		if (GetDataset() == 6) //Wen
		{
			//original values. nothing for nvtx 9. so I'll stuff 1 in there

			/*	                     val   +/- stat err
											bin, val, err = 1, 0.692896, 0.00590148
											bin, val, err = 2, 1.04487, 0.00946897
											bin, val, err = 3, 1.52813, 0.0193913
											bin, val, err = 4, 2.22268, 0.0455044
											bin, val, err = 5, 3.94164, 0.139866
											bin, val, err = 6, 10.6812, 0.738173
											bin, val, err = 7, 32.7236, 5.04018
											bin, val, err = 8, 185.434, 109.025
											bin, val, err = 10, 48.0754, 51.3947
											*/	
			//use central weights
			
				/* weights without MET cleanup
				vNvtxWgts.push_back(0.692896);
				vNvtxWgts.push_back(1.04487);
				vNvtxWgts.push_back(1.52813);
				vNvtxWgts.push_back(2.22268);
				vNvtxWgts.push_back(3.94164);
				vNvtxWgts.push_back(10.6812);
				vNvtxWgts.push_back(32.7236);
				vNvtxWgts.push_back(185.434);
				vNvtxWgts.push_back(1.0);
				vNvtxWgts.push_back(48.0754);
				*/

				// weights with MET cleanup
				/*vNvtxWgts.push_back(0.457523);
				vNvtxWgts.push_back(0.938594);
				vNvtxWgts.push_back(1.72987);
				vNvtxWgts.push_back(2.94801);
				vNvtxWgts.push_back(4.00134);
				vNvtxWgts.push_back(6.96233);
				vNvtxWgts.push_back(27.8493);
				*/

				//comparing w+jets data 11-04-2010
				/*vNvtxWgts.push_back(0.531763);
				vNvtxWgts.push_back(1.09089);
				vNvtxWgts.push_back(2.01056);
				vNvtxWgts.push_back(3.42637);
				vNvtxWgts.push_back(4.65062);
				vNvtxWgts.push_back(8.09207);
				vNvtxWgts.push_back(32.3683);
				*/
				//njet15==1 wgts
				/*
				vNvtxWgts.push_back(0.533008);
				vNvtxWgts.push_back(1.09741);
				vNvtxWgts.push_back(2.01858);
				vNvtxWgts.push_back(3.41263);
				vNvtxWgts.push_back(4.73342);
				vNvtxWgts.push_back(11.0714);
				vNvtxWgts.push_back(27.6785);
				*/

				//incl njet>=2
				vNvtxWgts.push_back(0.525801);
				vNvtxWgts.push_back(1.05841);
				vNvtxWgts.push_back(1.94882);
				vNvtxWgts.push_back(3.48312);
				vNvtxWgts.push_back(4.10328);
				vNvtxWgts.push_back(2.76971);


			
			if (UseNvtxWeights() == 2) // use central + stat err weights for systematics
			{
				/*
				vNvtxWgts.push_back(0.692896 + 0.00590148);
				vNvtxWgts.push_back(1.04487 + 0.00946897);
				vNvtxWgts.push_back(1.52813 + 0.0193913);
				vNvtxWgts.push_back(2.22268 + 0.0455044);
				vNvtxWgts.push_back(3.94164 + 0.139866);
				vNvtxWgts.push_back(10.6812 + 0.738173);
				vNvtxWgts.push_back(32.7236 + 5.04018);
				vNvtxWgts.push_back(185.434 + 109.025);
				vNvtxWgts.push_back(1.0 + 1.0); //100% error );
				vNvtxWgts.push_back(48.0754 + 51.3947);
				*/

				// weights with MET cleanup
				/*vNvtxWgts.push_back(0.457523 + 0.0132797);
				vNvtxWgts.push_back(0.938594 + 0.0268515);
				vNvtxWgts.push_back(1.72987 + 0.0679859);
				vNvtxWgts.push_back(2.94801 + 0.195354);
				vNvtxWgts.push_back(4.00134 + 0.53923);
				vNvtxWgts.push_back(6.96233 + 2.14862);
				vNvtxWgts.push_back(27.8493 + 28.9865);
				*/
				//comparing w+jets data 11-04-2010
				/*vNvtxWgts.push_back(0.531763+0.0154345);
				vNvtxWgts.push_back(1.09089+0.0312086);
				vNvtxWgts.push_back(2.01056+0.0790176);
				vNvtxWgts.push_back(3.42637+0.227053);
				vNvtxWgts.push_back(4.65062+0.626729);
				vNvtxWgts.push_back(8.09207+2.49727);
				vNvtxWgts.push_back(32.3683+33.69);
				*/
				//njet15==1
				/*vNvtxWgts.at(0) += 0.0167472;
				vNvtxWgts.at(1) += 0.033859;
				vNvtxWgts.at(2) += 0.0860775;
				vNvtxWgts.at(3) += 0.244785;
				vNvtxWgts.at(4) += 0.717349;
				vNvtxWgts.at(5) += 4.12607;
				vNvtxWgts.at(6) += 29.0294;
				*/

				//incl >=2 jet
				vNvtxWgts.at(0) += 0.0399206;
				vNvtxWgts.at(1) += 0.0810058;
				vNvtxWgts.at(2) += 0.197817;
				vNvtxWgts.at(3) += 0.605603;
				vNvtxWgts.at(4) += 1.20894;
				vNvtxWgts.at(5) += 1.67714;


			
			}	

		} else if (GetDataset() == 8) //wtn
		{
						
			//inclusive njet
			/*vNvtxWgts.push_back(0.541767);
			vNvtxWgts.push_back(1.09912);
			vNvtxWgts.push_back(2.26181);
			vNvtxWgts.push_back(2.75954);
			vNvtxWgts.push_back(1.20096);
			*/
			//njet15==1
			/*
			vNvtxWgts.push_back(0.515629);
			vNvtxWgts.push_back(1.02466);
			vNvtxWgts.push_back(3.65337);
			vNvtxWgts.push_back(3.33452);
			vNvtxWgts.push_back(1.81324);
			*/

			//incl >=2 jets
			vNvtxWgts.push_back(0.69517);
			vNvtxWgts.push_back(1.62533);
			vNvtxWgts.push_back(0.876305);
			vNvtxWgts.push_back(1.62533);
			vNvtxWgts.push_back(0.626632);
			/*
			vNvtxWgts.push_back(0.69517+0.287771);
			vNvtxWgts.push_back(1.62533+0.944019);
			vNvtxWgts.push_back(0.876305+0.443021);
			vNvtxWgts.push_back(1.62533+1.63509);
			vNvtxWgts.push_back(0.626632+0.636348);
			*/


			if (UseNvtxWeights() == 2) // use central + stat err weights for systematics
			{
				//inclusive njet
				/*vNvtxWgts.at(0) += 0.0848229;
				vNvtxWgts.at(1) += 0.221419;
				vNvtxWgts.at(2) += 0.802496;
				vNvtxWgts.at(3) += 1.59784;
				vNvtxWgts.at(4) += 0.85485;
				*/
				//njet15==1
				/*vNvtxWgts.at(0) += 0.0872093;
				vNvtxWgts.at(1) += 0.220091;
				vNvtxWgts.at(2) += 1.83052;
				vNvtxWgts.at(3) += 2.36328;
				vNvtxWgts.at(4) += 1.82091;
				*/
				assert (false && "syst for wtn is not implemented");
			}

		} else if (GetDataset() == 10) //ww
		{
				//inclusive njet
			/*vNvtxWgts.push_back(0.590653);
			vNvtxWgts.push_back(1.04473);
			vNvtxWgts.push_back(1.64022);
			vNvtxWgts.push_back(2.56923);
			vNvtxWgts.push_back(3.36269);
			vNvtxWgts.push_back(4.70777);
			vNvtxWgts.push_back(6.05284);
			vNvtxWgts.push_back(3.53082);
			*/

			//excl1jet
			/*vNvtxWgts.push_back(0.589689);
			vNvtxWgts.push_back(1.05429);
			vNvtxWgts.push_back(1.63497);
			vNvtxWgts.push_back(2.63068);
			vNvtxWgts.push_back(3.68145);
			vNvtxWgts.push_back(4.24303);
			vNvtxWgts.push_back(10.6076);
			vNvtxWgts.push_back(2.65189);
			*/
			//incl 2 jet
		vNvtxWgts.push_back(0.6054);
		vNvtxWgts.push_back(0.993155);
		vNvtxWgts.push_back(1.65674);
		vNvtxWgts.push_back(2.14001);
		vNvtxWgts.push_back(2.10017);
		vNvtxWgts.push_back(1.44386);
		/*
		vNvtxWgts.push_back(0.6054+0.0562463);
		vNvtxWgts.push_back(0.993155+0.0970091);
		vNvtxWgts.push_back(1.65674+0.224775);
		vNvtxWgts.push_back(2.14001+0.467692);
		vNvtxWgts.push_back(2.10017+0.734035);
		vNvtxWgts.push_back(1.44386+1.76837);
		*/

			if (UseNvtxWeights() == 2) // use central + stat err weights for systematics
			{
				//inclusive njet
				/*vNvtxWgts.at(0) += 0.0197296;
				vNvtxWgts.at(1) += 0.0359825;
				vNvtxWgts.at(2) += 0.0790681;
				vNvtxWgts.at(3) += 0.212802;
				vNvtxWgts.at(4) += 0.571547;
				vNvtxWgts.at(5) += 1.72924;
				vNvtxWgts.at(6) += 4.62293;
				vNvtxWgts.at(7) += 2.83096;
				*/

				//excl njet 1
				/*
				vNvtxWgts.at(0) += 0.0211632;
				vNvtxWgts.at(1) += 0.0388684;
				vNvtxWgts.at(2) += 0.0845184;
				vNvtxWgts.at(3) += 0.235566;
				vNvtxWgts.at(4) += 0.716573;
				vNvtxWgts.at(5) += 1.58128;
				vNvtxWgts.at(6) += 11.1253;
				vNvtxWgts.at(7) += 2.21873;
				*/
				assert (false && "syst for ww is not implemented");

			}


		} else if (GetDataset() == 11) //wz
		{
				//inclusive njet
			/*vNvtxWgts.push_back(0.581414);
			vNvtxWgts.push_back(1.04586);
			vNvtxWgts.push_back(1.72096);
			vNvtxWgts.push_back(2.57688);
			vNvtxWgts.push_back(3.52461);
			vNvtxWgts.push_back(4.25122);
			vNvtxWgts.push_back(7.89512);
			*/

			//njet15==1
			/*vNvtxWgts.push_back(0.567202);
			vNvtxWgts.push_back(1.07486);
			vNvtxWgts.push_back(1.76172);
			vNvtxWgts.push_back(2.73908);
			vNvtxWgts.push_back(3.59951);
			vNvtxWgts.push_back(16.4723);
			vNvtxWgts.push_back(4.57565);
			*/
			//incl 2 jet
			vNvtxWgts.push_back(0.593226);
			vNvtxWgts.push_back(0.965551);
			vNvtxWgts.push_back(1.70714);
			vNvtxWgts.push_back(2.46121);
			vNvtxWgts.push_back(4.34304);
			vNvtxWgts.push_back(1.92476);
			
			/*vNvtxWgts.push_back(0.593226+0.0439312);
			vNvtxWgts.push_back(0.965551+0.0690552);
			vNvtxWgts.push_back(1.70714+0.155427);
			vNvtxWgts.push_back(2.46121+0.351975);
			vNvtxWgts.push_back(4.34304+1.14669);
			vNvtxWgts.push_back(1.92476+0.976851);
			*/

			
			if (UseNvtxWeights() == 2) // use central + stat err weights for systematics
			{
				//inclusive njet
				/*vNvtxWgts.at(0) += 0.0185203;
				vNvtxWgts.at(1) += 0.0338238;
				vNvtxWgts.at(2) += 0.0777747;
				vNvtxWgts.at(3) += 0.194927;
				vNvtxWgts.at(4) += 0.551956;
				vNvtxWgts.at(5) += 1.34927;
				vNvtxWgts.at(6) += 6.02999;
				*/
				
				//njet15==1
				/*vNvtxWgts.at(0) += 0.02077;
				vNvtxWgts.at(1) += 0.0413827;
				vNvtxWgts.at(2) += 0.097705;
				vNvtxWgts.at(3) += 0.262733;
				vNvtxWgts.at(4) += 0.735991;
				vNvtxWgts.at(5) += 11.9669;
				vNvtxWgts.at(6) += 3.54428;
				*/
				assert (false && "syst for ww is not implemented");

			}

		} else if (GetDataset() == 12) //zz
		{
				//inclusive njet
			/*vNvtxWgts.push_back(0.587313);
			vNvtxWgts.push_back(1.05321);
			vNvtxWgts.push_back(1.66294);
			vNvtxWgts.push_back(2.39119);
			vNvtxWgts.push_back(3.53192);
			vNvtxWgts.push_back(5.43915);
			vNvtxWgts.push_back(9.32426);
			*/

			//njet15==1
			/*vNvtxWgts.push_back(0.573166);
			vNvtxWgts.push_back(1.05313);
			vNvtxWgts.push_back(1.78814);
			vNvtxWgts.push_back(2.6896);
			vNvtxWgts.push_back(3.77183);
			vNvtxWgts.push_back(7.28794);
			*/
		
			//incl njet2

			vNvtxWgts.push_back(0.593226);
			vNvtxWgts.push_back(0.965551);
			vNvtxWgts.push_back(1.70714);
			vNvtxWgts.push_back(2.46121);
			vNvtxWgts.push_back(4.34304);
			vNvtxWgts.push_back(1.92476);
			
			/*vNvtxWgts.push_back(0.593226+0.0439312);
			vNvtxWgts.push_back(0.965551+0.0690552);
			vNvtxWgts.push_back(1.70714+0.155427);
			vNvtxWgts.push_back(2.46121+0.351975);
			vNvtxWgts.push_back(4.34304+1.14669);
			vNvtxWgts.push_back(1.92476+0.976851);
			*/


			if (UseNvtxWeights() == 2) // use central + stat err weights for systematics
			{
				//inclusive njet
				/*vNvtxWgts.at(0) += 0.0206451;
				vNvtxWgts.at(1) += 0.0390207;
				vNvtxWgts.at(2) += 0.0876625;
				vNvtxWgts.at(3) += 0.212354;
				vNvtxWgts.at(4) += 0.679099;
				vNvtxWgts.at(5) += 2.37384;
				vNvtxWgts.at(6) += 9.705;
				*/
				//njet15==1
				/*vNvtxWgts.push_back(0.573166+0.0229915);
				vNvtxWgts.push_back(1.05313+0.0454644);
				vNvtxWgts.push_back(1.78814+0.115158);
				vNvtxWgts.push_back(2.6896+0.300805);
				vNvtxWgts.push_back(3.77183+0.932384);
				vNvtxWgts.push_back(7.28794+4.37951);
				*/
				assert (false && "syst for ww is not implemented");

			}

		} else if (GetDataset() == 13) //ttbar
		{
				//inclusive njet
			/*vNvtxWgts.push_back(0.535589);
			vNvtxWgts.push_back(1.10599);
			vNvtxWgts.push_back(2.20698);
			vNvtxWgts.push_back(2.77679);
			vNvtxWgts.push_back(1.61129);
			*/
			//njet15==1
			/*vNvtxWgts.push_back(0.573166);
			vNvtxWgts.push_back(1.05313);
			vNvtxWgts.push_back(1.78814);
			vNvtxWgts.push_back(2.6896);
			vNvtxWgts.push_back(3.77183);
			vNvtxWgts.push_back(7.28794);
			*/

			//incl njet2
			vNvtxWgts.push_back(0.529101);
			vNvtxWgts.push_back(1.03567);
			vNvtxWgts.push_back(2.28674);
			vNvtxWgts.push_back(2.96893);
			vNvtxWgts.push_back(1.90775);
			/*
			vNvtxWgts.push_back(0.529101+0.0570824);
			vNvtxWgts.push_back(1.03567+0.129538);
			vNvtxWgts.push_back(2.28674+0.464725);
			vNvtxWgts.push_back(2.96893+0.993808);
			vNvtxWgts.push_back(1.90775+0.848715);
			*/


			if (UseNvtxWeights() == 2) // use central + stat err weights for systematics
			{
				//inclusive njet
				/*vNvtxWgts.at(0) += 0.0433519;
				vNvtxWgts.at(1) += 0.113776;
				vNvtxWgts.at(2) += 0.389755;
				vNvtxWgts.at(3) += 0.810839;
				vNvtxWgts.at(4) += 0.670833;
				*/

				//njet15==1
				/*vNvtxWgts.at(0) += 0.0229915;
				vNvtxWgts.at(1) += 0.0454644;
				vNvtxWgts.at(2) += 0.115158;
				vNvtxWgts.at(3) += 0.300805;
				vNvtxWgts.at(4) += 0.932384;
				vNvtxWgts.at(5) += 4.37951;
				*/
				assert (false && "syst for ww is not implemented");

			}

		} else 
		{
			std::cout << "NO NVtx weights defined for this sample!" <<std::endl;
			assert( false);
		}

		
		std::cout << green << " ======= Nvtx weights ====== " << std::endl;
		for (int i=0; i < vNvtxWgts.size(); ++i)
		{
			std::cout <<  " vtx, w = " << i+1 << ", " << vNvtxWgts.at(i) << std::endl;
		}
		std::cout << clearatt << std::endl;
	}

	njet15counterr = 0;
	for (int i=0; i<11; ++i) iCount[i]=0;
} // Init

//-------------------------------------------------------------------
void WJetSel::CleanUp()
{
	vMcCount.clear();

	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void WJetSel::Main(TChain *ch, int iRunEvents)
{
	TBenchmark timer;
	timer.Start("phoana_time");


	if (ch) {
		myChain = ch;
		//ch->Print();	
	} else {
		std::cout << "NULL chain!" << std::endl;
		assert(false);
	}
	if (ch->GetEntries()<1)
	{
		std::cout << "No entries found!" << std::endl;
		return;
	}

	gROOT->Reset();
	Init(myChain);
	std::cout << yellow << "******* Job Settings ***********" << std::endl;
	PrintJobSettings();
	std::cout << "********************************" << clearatt << std::endl;

	BookHistograms();

	unsigned int iENTRIES = (unsigned int) myChain->GetEntries();
	unsigned int iEvt2Process = 0; 	//_____________________________________________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;			//______________________________________________ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//______________________________________________ default

	std::cout << " ==> Entries found :: " << iENTRIES << std::endl;

	//do a check for the MC. this is a common mistake I make

	readIn->GetEntry(1);
	if (stuple.evt_McFlag != GetMcFlag()) 
	{
		std::cout << red << "MC Flag," << GetMcFlag()
			<< ", setting does not match what is in Stuple, " 
			<< stuple.evt_McFlag << ". pleas check." 
			<< " returning!" << clearatt << std::endl;
		//std::cout << "Changing the MCFlag to Stuple setting McFlag = " << stuple.evt_McFlag << clearatt << std::endl;
		//SetMcFlag(stuple.evt_McFlag); //for some reason this seems to overwrite counting arrays? 12-07-2009
		//return; // simple return causes trouble when rerun. so drastic measure needed, even exit(1) would not work
		assert (false);
	}

	unsigned int iEvtProc = 0;			//______________________________________________ number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;
	std::cout << "Init done. Ready, Set , GO... " << std::endl; 

	timer.Start("looptimer");

	iCurrTree = -1;
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
		if (iCurrTree != myChain->GetTreeNumber())
		{
			std::cout << green << "Opening file " << myChain->GetFile()->GetName() << clearatt << std::endl;
			iCurrTree = myChain->GetTreeNumber();
		}
		readIn->GetEntry(iEvtProc);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed [" << (int)( iEvtProc/(double)iEvt2Process * 100)<< "%] ";
			timer.Show("looptimer");
			timer.Start("looptimer");
		}

		//___________________________________________________________________________ drop the first 400pb-1
		if (stuple.evt_McFlag == 0) {  //______________________________________________ if data
			if (stuple.evt_RunNumber < 190851) {
				iFirst400++;
				continue;
			}
		}
		++iCount[0];

		if (stuple.ele_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
		else if (stuple.ele_num > 0 && stuple.ele_num < 1) iEleEvts++;
		else if (stuple.ele_num < 1 && stuple.ele_num > 0) iPhoEvts++;
		else iOther++;

		//apply some strip cuts here to speed up things a bit
		//the met cut and eletron is required
		//met is not calculated again in this so this should be fine
		//2-9-2010
		if (stuple.met_Met < fMinMet || stuple.met_Met > fMaxMet) 
		{
			//std::cout << "\t failed met cut " << std::endl;
			continue;
		} else ++iCount[iMet_cut]; //met cut
		if (stuple.ele_num<1) 
		{
			continue;
		} else ++iCount[ielenum_cut]; //ele cut

		CommonVars commVars;
		ElectronList centralEles;

		GetCommonVariables(stuple, commVars);
		CreateElectronLists(stuple, centralEles,0);

		DoMyStuff(commVars, centralEles);

	}

	//print out the job summary
	std::cout << magenta << "======== WJetSel =======================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Entries request  = " << iRunEvents << " (" << iRunEvents/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Events Processed = " << iEvtProc << " (" << iEvtProc/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	std::cout << red << "__________________________________________" << std::endl;

	PrintJobSettings();

	const int itextWidth = 25;
	std::cout << "Njet15counterr Evts  = " << njet15counterr << std::endl;
	std::cout << "W Evts               = " << vMcCount[0] << std::endl;
	std::cout << "W+>=1Jet Evts        = " << vMcCount[1] << std::endl; 
	std::cout << cyan << "Data Written to        = " << GetHistFileName() << clearatt << std::endl;
	std::cout << "======= END WJetSel ====================" << std::endl;

	for (int i=0; i<11; ++i) 
	{
		std::cout << "iCount["<< i << "] = " << iCount[i];
		std::string cut("");
		if (iMet_cut == i) cut += " Met cut";
		else if (ielenum_cut == i) cut += " elenum cut";
		else if (ivtx_cut == i) cut += " vtx cut";
		else if (ivtxz_cut == i) cut += " vtxz cut";
		else if (itele_cut == i) cut += " tele cut";
		else if (iWpt_cut == i) cut += " Wpt cut";
		else if (iWeta_cut == i) cut += " Weta cut";
		else if (iWmass_cut == i) cut += " Wmass cut";
		else if (iWjet_cut == i) cut += " Wjet cut";
		else if (iMetCleanup_cut == i) cut += " MET cleanup cut";
		std::cout << cut << std::endl;
	}

	EVT.clear();

	CleanUp();
	timer.Show("phoana_time");
} // Main


//-------------------------------------------------------------------
void WJetSel::BookHistograms()
{

	/*Hsignal.z1j_zpt	= new TH1F("z1j_zpt","#gamma+>=1Jet :- p_{T}^{Z} ",200,-500,500);
	  Hsignal.z1j_zet	= new TH1F("z1j_zet","#gamma+>=1Jet :- E_{T}^{Z} ",250,0,500);
	  Hsignal.z1j_zmass	= new TH1F("z1j_zmass","#gamma+>=1Jet :- Z mass ",200,0,1000);
	  Hsignal.z1j_zeta	= new TH1F("z1j_zeta","#gamma+>=1Jet :- Event #Eta^{Z} ",50,-5,5);
	  Hsignal.z1j_zj_mass	= new TH1F("z1j_zj1mass","#gamma+>=1Jet :- Invariant Mass (Z, Lead Jet)",200,0,1000);

	  Hsignal.z2j_zpt	= new TH1F("z2j_zpt","#gamma+>=2Jets :- p_{T}^{Z} ",200,-500,500);
	  Hsignal.z2j_zet	= new TH1F("z2j_zet","#gamma+>=2Jets :- E_{T}^{Z} ",250,0,500);
	  Hsignal.z2j_zmass	= new TH1F("z2j_zmass","#gamma+>=2Jets :- Z mass ",200,0,1000);
	  Hsignal.z2j_zeta	= new TH1F("z2j_zeta","#gamma+>=2Jets :- Event #Eta^{Z} ",50,-5,5);
	  Hsignal.z2j_z1j_mass	= new TH1F("z2j_zj1mass","#gamma+>=2Jets :- Invariant Mass (Z, Lead Jet)",200,0,1000);
	  Hsignal.z2j_z2j_mass	= new TH1F("z2j_zj2mass","#gamma+>=2Jets :- Invariant Mass (Z, Second Lead Jet)",200,0,1000);
	  Hsignal.z2j_zj1j2_mass	= new TH1F("z2j_zj1j2mass","#gamma+>=2Jets :- Invariant Mass (Z, Two Lead Jets)",200,0,1000);
	  */

	Histograms HistoMan;
	std::string text2title1, text2title2;

	TDirectory *l2dir, *l3dir, *l4dir;

	//CENTRAL
	text2title1 = "W+Jet :";
	text2title2 = "W+Jet :";

	topDir->cd();
		std::stringstream title1;
		title1 <<  text2title1 << " - Closest Jet Et>15GeV before cut: #delta#phi(jet,#slash{E}_{T})";
		hClosestJetMetDelPhi_a4	= new TH1F("ClosesetJetMetDelPhi_a4", title1.str().c_str(), 700, -3.5, 3.5);
		hClosestJetMetDelPhi_b4	= new TH1F("ClosesetJetMetDelPhi_b4", title1.str().c_str(), 700, -3.5, 3.5);
		hClosestJetPt_a4	= new TH1F("ClosesetJetPt_a4", title1.str().c_str(), 500, 0, 500);
	

	l2dir = gDirectory->mkdir("Event","Event Histograms");
	l2dir->cd();
	HistoMan.GetEventHistograms(Hmc.Evt,text2title1.c_str());

	topDir->cd();
	l2dir = gDirectory->mkdir("Electron1","Electron-1 Histograms");
	l2dir->cd();
	HistoMan.GetPhotonHistograms(Hmc.Ele1,text2title1.c_str());

	topDir->cd();
	l2dir = gDirectory->mkdir("Jet","Jet Histograms");
	l2dir->cd();
	HistoMan.GetJetHistograms(Hmc.Jet,text2title1.c_str());
	
	topDir->cd();
	l2dir = gDirectory->mkdir("Jet2","Second Jet Histograms");
	l2dir->cd();
	HistoMan.GetJetHistograms(Hmc.Jet2,text2title1.c_str());

	topDir->cd();
	l2dir = gDirectory->mkdir("W","W Histograms");
	l2dir->cd();
	HistoMan.GetTwoJetsHistograms(Hmc.W,text2title2.c_str());
	
		std::stringstream title2;
		title2 <<  text2title1 << "#delta#phi(W, #slash{E}_{T})";
		hWMetDelPhi_b4	= new TH1F("WMetDelPhi_b4",title2.str().c_str(), 700, -3.5, 3.5);
		hWMetDelPhi_a4	= new TH1F("WMetDelPhi_a4",title2.str().c_str(), 700, -3.5, 3.5);


	topDir->cd();
	l2dir = gDirectory->mkdir("WJet","W+Jet Histograms");
	l2dir->cd();
	HistoMan.GetTwoJetsHistograms(Hmc.WJet,text2title2.c_str());
} //BookHistograms

void WJetSel::CreateElectronLists(const Stuple stuple, ElectronList& list, const int collection)
{
	/*{{{*/
	switch (collection)
	{
		case 0:
			{
				list.ele_num    = stuple.ele_num;
				list.ele_Ntight = stuple.ele_Ntight;
				list.ele_Nloose = stuple.ele_Nloose;

				for (unsigned int i=0; i < stuple.ele_num; i++)
				{
					list.ele_Index.push_back(stuple.ele_Index[i]);
					list.ele_PhoBlockIndex.push_back(stuple.ele_PhoBlockIndex[i]);
					list.ele_EleBlockIndex.push_back(stuple.ele_EleBlockIndex[i]);
					list.ele_Etc.push_back(stuple.ele_Etc[i]);
					list.ele_E.push_back(stuple.ele_E[i]);
					list.ele_Px.push_back(stuple.ele_Px[i]);
					list.ele_Py.push_back(stuple.ele_Py[i]);
					list.ele_Pz.push_back(stuple.ele_Pz[i]);
					list.ele_Detector.push_back(stuple.ele_Detector[i]);
					list.ele_DetEta.push_back(stuple.ele_DetEta[i]);
					list.ele_DetPhi.push_back(stuple.ele_DetPhi[i]);
					list.ele_XCes.push_back(stuple.ele_XCes[i]);
					list.ele_ZCes.push_back(stuple.ele_ZCes[i]);
					list.ele_HadEm.push_back(stuple.ele_HadEm[i]);
					list.ele_Chi2Mean.push_back(stuple.ele_Chi2Mean[i]);
					list.ele_N3d.push_back(stuple.ele_N3d[i]);
					list.ele_Iso4.push_back(stuple.ele_Iso4[i]);
					list.ele_TrkIso.push_back(stuple.ele_TrkIso[i]);		
					list.ele_CesWireE2.push_back(stuple.ele_CesWireE2[i]);
					list.ele_CesStripE2.push_back(stuple.ele_CesStripE2[i]);
					list.ele_PhiWedge.push_back(stuple.ele_PhiWedge[i]);
					list.ele_NMuonStubs.push_back(stuple.ele_NMuonStubs[i]);
					list.ele_EmTime.push_back(stuple.ele_EmTime[i]);
					list.ele_PhoenixId.push_back(stuple.ele_PhoenixId[i]);
					list.ele_Halo_seedWedge.push_back(stuple.ele_Halo_seedWedge[i]);
					list.ele_Halo_eastNhad.push_back(stuple.ele_Halo_eastNhad[i]);
					list.ele_Halo_westNhad.push_back(stuple.ele_Halo_westNhad[i]);
					list.ele_matchJetIndex.push_back(stuple.ele_matchJetIndex[i]);
					list.ele_Ntracks.push_back(stuple.ele_Ntracks[i]);	
					list.ele_Emfr.push_back(stuple.ele_Emfr[i]);
					list.ele_EoverP.push_back(stuple.ele_EoverP[i]);
					list.ele_TrackPt.push_back(stuple.ele_TrackPt[i]);
					list.ele_TrackBcPt.push_back(stuple.ele_TrackBcPt[i]);
					list.ele_TrackPhi.push_back(stuple.ele_TrackPhi[i]);
					list.ele_Nssl.push_back(stuple.ele_Nssl[i]);
					list.ele_Nasl.push_back(stuple.ele_Nasl[i]);
					list.ele_TightId.push_back(stuple.ele_TightId[i]);
					list.ele_LooseId.push_back(stuple.ele_LooseId[i]);
					list.ele_ConversionId.push_back(stuple.ele_ConversionId[i]);
				}

				assert( (int)list.ele_num == (int)list.ele_Index.size());
				break;
			}

		case 1:
			{
				list.ele_num    = stuple.ele_up_num;
				list.ele_Ntight = stuple.ele_up_Ntight;
				list.ele_Nloose = stuple.ele_up_Nloose;

				for (unsigned int i=0; i < stuple.ele_up_num; i++)
				{
					list.ele_Index.push_back(stuple.ele_up_Index[i]);
					list.ele_PhoBlockIndex.push_back(stuple.ele_up_PhoBlockIndex[i]);
					list.ele_EleBlockIndex.push_back(stuple.ele_up_EleBlockIndex[i]);
					list.ele_Etc.push_back(stuple.ele_up_Etc[i]);
					list.ele_E.push_back(stuple.ele_up_E[i]);
					list.ele_Px.push_back(stuple.ele_up_Px[i]);
					list.ele_Py.push_back(stuple.ele_up_Py[i]);
					list.ele_Pz.push_back(stuple.ele_up_Pz[i]);
					list.ele_Detector.push_back(stuple.ele_up_Detector[i]);
					list.ele_DetEta.push_back(stuple.ele_up_DetEta[i]);
					list.ele_DetPhi.push_back(stuple.ele_up_DetPhi[i]);
					list.ele_XCes.push_back(stuple.ele_up_XCes[i]);
					list.ele_ZCes.push_back(stuple.ele_up_ZCes[i]);
					list.ele_HadEm.push_back(stuple.ele_up_HadEm[i]);
					list.ele_Chi2Mean.push_back(stuple.ele_up_Chi2Mean[i]);
					list.ele_N3d.push_back(stuple.ele_up_N3d[i]);
					list.ele_Iso4.push_back(stuple.ele_up_Iso4[i]);
					list.ele_TrkIso.push_back(stuple.ele_up_TrkIso[i]);		
					list.ele_CesWireE2.push_back(stuple.ele_up_CesWireE2[i]);
					list.ele_CesStripE2.push_back(stuple.ele_up_CesStripE2[i]);
					list.ele_PhiWedge.push_back(stuple.ele_up_PhiWedge[i]);
					list.ele_NMuonStubs.push_back(stuple.ele_up_NMuonStubs[i]);
					list.ele_EmTime.push_back(stuple.ele_up_EmTime[i]);
					list.ele_PhoenixId.push_back(stuple.ele_up_PhoenixId[i]);
					list.ele_Halo_seedWedge.push_back(stuple.ele_up_Halo_seedWedge[i]);
					list.ele_Halo_eastNhad.push_back(stuple.ele_up_Halo_eastNhad[i]);
					list.ele_Halo_westNhad.push_back(stuple.ele_up_Halo_westNhad[i]);
					list.ele_matchJetIndex.push_back(stuple.ele_up_matchJetIndex[i]);
					list.ele_Ntracks.push_back(stuple.ele_up_Ntracks[i]);	
					list.ele_Emfr.push_back(stuple.ele_up_Emfr[i]);
					list.ele_EoverP.push_back(stuple.ele_up_EoverP[i]);
					list.ele_TrackPt.push_back(stuple.ele_up_TrackPt[i]);
					list.ele_TrackBcPt.push_back(stuple.ele_up_TrackBcPt[i]);
					list.ele_TrackPhi.push_back(stuple.ele_up_TrackPhi[i]);
					list.ele_Nssl.push_back(stuple.ele_up_Nssl[i]);
					list.ele_Nasl.push_back(stuple.ele_up_Nasl[i]);
					list.ele_TightId.push_back(stuple.ele_up_TightId[i]);
					list.ele_LooseId.push_back(stuple.ele_up_LooseId[i]);
					list.ele_ConversionId.push_back(stuple.ele_up_ConversionId[i]);
				}

				assert(list.ele_num == list.ele_Index.size());
				break;

			}

		case -1:
			{
				list.ele_num    = stuple.ele_down_num;
				list.ele_Ntight = stuple.ele_down_Ntight;
				list.ele_Nloose = stuple.ele_down_Nloose;

				for (unsigned int i=0; i < stuple.ele_down_num; i++)
				{
					list.ele_Index.push_back(stuple.ele_down_Index[i]);
					list.ele_PhoBlockIndex.push_back(stuple.ele_down_PhoBlockIndex[i]);
					list.ele_EleBlockIndex.push_back(stuple.ele_down_EleBlockIndex[i]);
					list.ele_Etc.push_back(stuple.ele_down_Etc[i]);
					list.ele_E.push_back(stuple.ele_down_E[i]);
					list.ele_Px.push_back(stuple.ele_down_Px[i]);
					list.ele_Py.push_back(stuple.ele_down_Py[i]);
					list.ele_Pz.push_back(stuple.ele_down_Pz[i]);
					list.ele_Detector.push_back(stuple.ele_down_Detector[i]);
					list.ele_DetEta.push_back(stuple.ele_down_DetEta[i]);
					list.ele_DetPhi.push_back(stuple.ele_down_DetPhi[i]);
					list.ele_XCes.push_back(stuple.ele_down_XCes[i]);
					list.ele_ZCes.push_back(stuple.ele_down_ZCes[i]);
					list.ele_HadEm.push_back(stuple.ele_down_HadEm[i]);
					list.ele_Chi2Mean.push_back(stuple.ele_down_Chi2Mean[i]);
					list.ele_N3d.push_back(stuple.ele_down_N3d[i]);
					list.ele_Iso4.push_back(stuple.ele_down_Iso4[i]);
					list.ele_TrkIso.push_back(stuple.ele_down_TrkIso[i]);		
					list.ele_CesWireE2.push_back(stuple.ele_down_CesWireE2[i]);
					list.ele_CesStripE2.push_back(stuple.ele_down_CesStripE2[i]);
					list.ele_PhiWedge.push_back(stuple.ele_down_PhiWedge[i]);
					list.ele_NMuonStubs.push_back(stuple.ele_down_NMuonStubs[i]);
					list.ele_EmTime.push_back(stuple.ele_down_EmTime[i]);
					list.ele_PhoenixId.push_back(stuple.ele_down_PhoenixId[i]);
					list.ele_Halo_seedWedge.push_back(stuple.ele_down_Halo_seedWedge[i]);
					list.ele_Halo_eastNhad.push_back(stuple.ele_down_Halo_eastNhad[i]);
					list.ele_Halo_westNhad.push_back(stuple.ele_down_Halo_westNhad[i]);
					list.ele_matchJetIndex.push_back(stuple.ele_down_matchJetIndex[i]);
					list.ele_Ntracks.push_back(stuple.ele_down_Ntracks[i]);	
					list.ele_Emfr.push_back(stuple.ele_down_Emfr[i]);
					list.ele_EoverP.push_back(stuple.ele_down_EoverP[i]);
					list.ele_TrackPt.push_back(stuple.ele_down_TrackPt[i]);
					list.ele_TrackBcPt.push_back(stuple.ele_down_TrackBcPt[i]);
					list.ele_TrackPhi.push_back(stuple.ele_down_TrackPhi[i]);
					list.ele_Nssl.push_back(stuple.ele_down_Nssl[i]);
					list.ele_Nasl.push_back(stuple.ele_down_Nasl[i]);
					list.ele_TightId.push_back(stuple.ele_down_TightId[i]);
					list.ele_LooseId.push_back(stuple.ele_down_LooseId[i]);
					list.ele_ConversionId.push_back(stuple.ele_down_ConversionId[i]);
				}

				assert(list.ele_num == list.ele_Index.size());
				break;
			}

		default:
			{
				StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central/em up/em down  lists only!");
				assert(false);
			}
	}
	/*}}}*/
}	//CreateElectronLists

void WJetSel::CreatePhotonLists(const Stuple stuple, PhotonList& pholist, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
			pholist.pho_num 	 =  stuple.pho_num;
			pholist.pho_Ntight = stuple.pho_Ntight;
			pholist.pho_Nloose = stuple.pho_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_E[i]);
				pholist.pho_Px.push_back(stuple.pho_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_matchJetIndex[i]);
				
				pholist.pho_CprWgt.push_back(stuple.pho_CprWgt[i]);
				pholist.pho_CprSys1.push_back(stuple.pho_CprSys1[i]);
				pholist.pho_CprSys2.push_back(stuple.pho_CprSys2[i]);
				pholist.pho_CprSys3.push_back(stuple.pho_CprSys3[i]);
				pholist.pho_CprSys4.push_back(stuple.pho_CprSys4[i]);
				pholist.pho_CprSys5.push_back(stuple.pho_CprSys5[i]);
				pholist.pho_CprSys6.push_back(stuple.pho_CprSys6[i]);
				pholist.pho_CprSys7.push_back(stuple.pho_CprSys7[i]);
				pholist.pho_CprSys8.push_back(stuple.pho_CprSys8[i]);
				
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case 1:
		{
			pholist.pho_num 	 =  stuple.pho_up_num;
			pholist.pho_Ntight = stuple.pho_up_Ntight;
			pholist.pho_Nloose = stuple.pho_up_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_up_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_up_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_up_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_up_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_up_E[i]);
				pholist.pho_Px.push_back(stuple.pho_up_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_up_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_up_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_up_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_up_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_up_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_up_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_up_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_up_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_up_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_up_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_up_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_up_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_up_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_up_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_up_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_up_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_up_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_up_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_up_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_up_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_up_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_up_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_up_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_up_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_up_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case -1:
		{
			pholist.pho_num 	 =  stuple.pho_down_num;
			pholist.pho_Ntight = stuple.pho_down_Ntight;
			pholist.pho_Nloose = stuple.pho_down_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_down_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_down_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_down_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_down_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_down_E[i]);
				pholist.pho_Px.push_back(stuple.pho_down_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_down_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_down_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_down_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_down_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_down_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_down_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_down_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_down_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_down_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_down_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_down_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_down_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_down_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_down_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_down_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_down_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_down_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_down_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_down_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_down_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_down_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_down_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_down_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_down_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_down_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid photon collection requested! can make central/em up/em down photon lists only!");
			assert(false);

		}
		
	}	//switch

/*}}}*/
}	//CreatePhotonLists
void WJetSel::CreateJetLists(const Stuple stuple, JetList& list, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
		  list.jet_num  = stuple.jet_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_NJet15;
			for (unsigned int i=0; i < stuple.jet_num; i++)
			{
				if ( fabs(stuple.jet_Px[i]) > 1500 || 
					  fabs(stuple.jet_Py[i]) > 1500 || 
					  fabs(stuple.jet_Pz[i]) > 1500 || 
					  fabs(stuple.jet_E[i]) > 1500)
				{
					stuple.DumpJetBlock();
					assert(false);
				}

				list.jet_Index.push_back(stuple.jet_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_E[i]);
				list.jet_Px.push_back(stuple.jet_Px[i]);
				list.jet_Py.push_back(stuple.jet_Py[i]);
				list.jet_Pz.push_back(stuple.jet_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_SeedIEta[i]);
			}

			assert(list.jet_num == list.jet_Index.size());

			list.jet_raw_num  = stuple.jet_raw_num;							// number of raw jets
			for (unsigned int i=0; i < stuple.jet_raw_num; i++)
			{
				list.jet_raw_Index.push_back(stuple.jet_raw_Index[i]);				//original index in jet block
				list.jet_raw_Pt.push_back(stuple.jet_raw_Pt[i]);
				list.jet_raw_E.push_back(stuple.jet_raw_E[i]);
				list.jet_raw_Px.push_back(stuple.jet_raw_Px[i]);
				list.jet_raw_Py.push_back(stuple.jet_raw_Py[i]);
				list.jet_raw_Pz.push_back(stuple.jet_raw_Pz[i]);
			}
			
			break;
		}

		case 1:
		{
		  list.jet_num  = stuple.jet_up_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_up_NJet15;
			for (unsigned int i=0; i < stuple.jet_up_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_up_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_up_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_up_E[i]);
				list.jet_Px.push_back(stuple.jet_up_Px[i]);
				list.jet_Py.push_back(stuple.jet_up_Py[i]);
				list.jet_Pz.push_back(stuple.jet_up_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_up_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_up_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_up_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_up_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_up_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_up_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_up_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_up_SeedIEta[i]);
				
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

			
		}

		case -1:
		{
		  list.jet_num  = stuple.jet_down_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_down_NJet15;
			for (unsigned int i=0; i < stuple.jet_down_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_down_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_down_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_down_E[i]);
				list.jet_Px.push_back(stuple.jet_down_Px[i]);
				list.jet_Py.push_back(stuple.jet_down_Py[i]);
				list.jet_Pz.push_back(stuple.jet_down_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_down_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_down_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_down_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_down_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_down_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_down_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_down_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_down_SeedIEta[i]);
				
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central and JES up/down lists only!");
			assert(false);

		}

	}
/*}}}*/
} // CreateJetLists
void WJetSel::GetCommonVariables(const Stuple stuple, CommonVars& cVars)
{
	/*{{{*/
	cVars.evt_McFlag 		= stuple.evt_McFlag;
	cVars.evt_RunNumber 	= stuple.evt_RunNumber;
	cVars.evt_EventNumber= stuple.evt_EventNumber;
	cVars.tri_pho25iso 	= stuple.tri_pho25iso;
	cVars.tri_pho50 		= stuple.tri_pho50;
	cVars.tri_pho70 		= stuple.tri_pho70;

	cVars.vtx_N 			= stuple.vtx_N;
	cVars.vtx_NClass12 	= stuple.vtx_NClass12;
	cVars.vtx_z 			= stuple.vtx_z;
	cVars.vtx_Ntracks 	= stuple.vtx_Ntracks;
	cVars.vtx_SumPt 		= stuple.vtx_SumPt;

	cVars.met_Met 			= stuple.met_Met;	
	cVars.met_MetX 		= stuple.met_MetX;	
	cVars.met_MetY 		= stuple.met_MetY;	
	cVars.met_SumEt 		= stuple.met_SumEt;
	cVars.met_Ht 			= stuple.met_Ht;	
	cVars.met_MetPhi 		= stuple.met_MetPhi;
	//	cVars.met_Gen_d 		= stuple.met_Gen_d;	
	//	cVars.met_Gen_m 		= stuple.met_Gen_m;
	//	cVars.met_Gen_p 		= stuple.met_Gen_p;
	//	cVars.met_Gen_mUn 	= stuple.met_Gen_mUn;
	//	cVars.met_Gen_pUn 	= stuple.met_Gen_pUn;
	/*}}}*/
}


void WJetSel::DoMyStuff(CommonVars cVars, ElectronList cEles)
{	
	//_________________________________________________________________ require a one good vertex 		  
	if (cVars.vtx_NClass12 < 1) return; 
	else ++iCount[ivtx_cut];
	if (fabs(cVars.vtx_z) > 60.) return;
	else ++iCount[ivtxz_cut];


	std::vector<int> UsedEle, UsedPho;

	//______________________________________________________________ select the electron
	for (unsigned int i = 0; i < cEles.ele_num; ++i )
	{
		//if (cEles.ele_Etc[i]<20) continue;
		//if (fabs(cEles.ele_DetEta.at(i))>1.1) continue;
		if (cEles.ele_ConversionId[i] != 0) continue;
		//if (UseTightPholikeEle())  //signal ele+jets  // can't apply cuts to select sideband pho-like ele
																		// as I do not have 2nd track info?
		//{
			if (cEles.ele_TightId[i] == 0 && UsedEle.size()<1)
			{
				UsedEle.push_back(i);
				break;
			}
		
		//} else   //sideband ele+jets selection
		//{
		//	if (cEles.ele_TightId[i] != 0 && IsSidebandPhoLikeElectron(cEles, i) && UsedEle.size()<1)
		//	{
		//		UsedEle.push_back(i);
		//		break;
		//	}
		//}
	}	

	int used = 0;				//____ to make sure evts are not fallen into multiple bins

	EvtTag_t evt;
	evt.Run 		= cVars.evt_RunNumber;
	evt.Evt 		= cVars.evt_EventNumber;

	if (UsedEle.size()>=1)
	{
		++iCount[itele_cut];
		const int iEle1Index = UsedEle.at(0);
		const TLorentzVector tlEle1(cEles.ele_Px[iEle1Index], cEles.ele_Py[iEle1Index], cEles.ele_Pz[iEle1Index], cEles.ele_E[iEle1Index]);
		//met cut is already applied as a strip cut
		const TLorentzVector tlMet(cVars.met_MetX,cVars.met_MetY,0,cVars.met_Met);
		const TLorentzVector tlWvec = tlEle1+tlMet;

		if (tlWvec.Pt()< GetMinWPt()) return;  //to match the photon Pt
		else ++iCount[iWpt_cut]; 
		if (fabs(tlWvec.Eta())<1.1) return;  //keep it central like the photon
		else ++iCount[iWeta_cut]; 
		if (tlWvec.M()<60. ||tlWvec.M()>100.) return; //mass cut to reduce background
		else ++iCount[iWmass_cut]; 

		++vMcCount[0];
		JetList cJets;
		PhotonList cPhos;
		CreatePhotonLists(stuple, cPhos,0);
		CreateJetLists(stuple, cJets,0);

		NewJetList nj;		// this need to be worked out to fit the new changes and new jet lists
		nj.AddUnused(cPhos, cEles, cJets, UsedPho, UsedEle, fMinJetEt, fMaxJetEta);
		if (cJets.jet_NJet15 <2) return;
		else ++iCount[iWjet_cut]; 

		//test debug
		int njet15 =0;
		for (int i=0; i < cJets.jet_Pt.size(); ++i)
		{		
			if (cJets.jet_Pt.at(i)>15.) ++njet15;
		}
		std::stringstream msg;
		msg <<  "njet15 does not match " << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber;
		if (njet15 != cJets.jet_NJet15)
		{
			std::cout << " njet15 does not match"; cVars.PrintHeader();
			stuple.Dump(0);
			cJets.Dump();
			//assert(false);
			++njet15counterr;
			return;
		}
		//assert (njet15 == cJets.jet_NJet15 && msg.str().c_str());

		//MET clean up cut. require MET to be away from any jet above 15GeV.
		const float fDelPhiCut = 0.4;
		bool bad = false;
		float dphi_jm_closest = 10.0;
		float closest_jetPt = -99999.;
		for (int i=0; i <cJets.jet_num ; ++i)
		{
			TLorentzVector tlJetVec(cJets.jet_Px[i], cJets.jet_Py[i], cJets.jet_Pz[i], cJets.jet_E[i]);
			if (tlJetVec.Pt()< GetMinJetEt() || fabs(tlJetVec.Eta())> GetMaxJetEta()) continue;
			const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tlMet.Phi())-TVector2::Phi_0_2pi(tlJetVec.Phi())));
			//std::cout << " ijet, dphi = " << i << ", " << dphi_jm << std::endl;
			if (dphi_jm <0.4 || dphi_jm>(TMath::Pi()-0.4))
			{
				//std::cout << red << "Failed MET CLEAN UP : jet Pt = " << tlJetVec.Pt() << clearatt << std::endl;
				bad = true;
			}
			if (dphi_jm < dphi_jm_closest) 
			{
				dphi_jm_closest = dphi_jm;
				closest_jetPt = tlJetVec.Pt();
			}
		}

		const float dphi_wm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tlMet.Phi())-TVector2::Phi_0_2pi(tlWvec.Phi())));
		hWMetDelPhi_b4->Fill(dphi_wm);
		hClosestJetMetDelPhi_b4->Fill(dphi_jm_closest);

		if (bad) return;
		else ++iCount[iMetCleanup_cut];

		hWMetDelPhi_a4->Fill(dphi_wm);
		hClosestJetMetDelPhi_a4->Fill(dphi_jm_closest);
		hClosestJetPt_a4->Fill(closest_jetPt);
		used = JetSelAndHistFill(cVars, cEles, cJets, UsedEle, &Hmc, vMcCount);
	}

}

int WJetSel::JetSelAndHistFill(const CommonVars Vars, 
		const ElectronList eles, const JetList jets, 
		std::vector<int> vEleIndex, 
		Hist_t* hist,
		std::vector<unsigned int>& count, 
		const float fWgt,
		const bool debug
		)
{
	
	int iEle1Index = vEleIndex[0];	//assuming the electrons are sorted in Et

	++count[1];
	float fWeight = 1;
	if (UseNvtxWeights())
	{
		if (Vars.vtx_NClass12 <= vNvtxWgts.size()) fWeight = vNvtxWgts.at(Vars.vtx_NClass12 - 1);
	}

	
	FillElectronHists(Vars, eles, &(hist->Ele1), iEle1Index, fWeight);
	FillJetHists(Vars, jets, &(hist->Jet),0, fWeight);		//lead jet

	if (jets.jet_NJet15>=2)	FillJetHists(Vars, jets, &(hist->Jet2),1, fWeight);		//second lead jet
	FillEventHists(Vars, eles, jets, &(hist->Evt), fWeight);
	FillWHists(Vars, eles, &(hist->W), vEleIndex.at(0), fWeight);
	FillWJetHists(Vars, eles, jets, &(hist->WJet), vEleIndex.at(0), fWeight);

	return 1;		//1=this event has been used!
}


void WJetSel::FillEventHists(const CommonVars& vars,
								const ElectronList& eles, const JetList& jets, 
								Histograms::EventHists_t* hist,
								const float fWeight)
{ 
	//if (!hist)
	//{
	//	StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
	//	assert(false);
	//}

	hist->Met->Fill(vars.met_Met, fWeight);
	hist->MetX->Fill(vars.met_MetX, fWeight);
	hist->MetY->Fill(vars.met_MetY, fWeight);
	hist->Met_Gen_d->Fill(vars.met_Gen_d, fWeight);
	hist->Sumet->Fill(vars.met_SumEt, fWeight);
	hist->NVertices->Fill(vars.vtx_N, fWeight);
	hist->N12Vertices->Fill(vars.vtx_NClass12, fWeight);
	hist->VertexZ->Fill(vars.vtx_z, fWeight);
	hist->MetRun->Fill(vars.evt_RunNumber, vars.met_Met, fWeight);
	hist->SumetRun->Fill(vars.evt_RunNumber, vars.met_SumEt, fWeight);
	hist->NVerticesRun->Fill(vars.evt_RunNumber, vars.vtx_N, fWeight);
	hist->N12VerticesRun->Fill(vars.evt_RunNumber, vars.vtx_NClass12, fWeight);
	hist->N12VerticesVsMet->Fill(vars.met_Met, vars.vtx_NClass12, fWeight);
	hist->N12VerticesVsSumet->Fill(vars.met_SumEt, vars.vtx_NClass12, fWeight);
	hist->N12VerticesMetProf->Fill(vars.met_Met, vars.vtx_NClass12, fWeight);
	hist->N12VerticesSumetProf->Fill(vars.met_SumEt, vars.vtx_NClass12, fWeight);
	if (vars.tri_pho25iso == 1) hist->Triggers->Fill(0);
	if (vars.tri_pho50 == 1) hist->Triggers->Fill(1);
	if (vars.tri_pho70 == 1) hist->Triggers->Fill(2);

	hist->Ht->Fill(vars.met_Ht, fWeight);
	hist->NJet15->Fill(jets.jet_NJet15, fWeight);

}

void WJetSel::FillElectronHists(const CommonVars& cVars, const ElectronList& eles, Histograms::PhotonHists_t* hist,
								const int iEleIndex, const float fWeight, const bool debug)
{
	hist->Detector->Fill(eles.ele_Detector[iEleIndex], fWeight);
	hist->DetEta->Fill(eles.ele_DetEta[iEleIndex], fWeight);
	hist->DetPhi->Fill(eles.ele_DetPhi[iEleIndex], fWeight);
	hist->EtCorr->Fill(eles.ele_Etc[iEleIndex], fWeight);
	hist->EventWeights->Fill(eles.ele_Etc[iEleIndex], fWeight);
	hist->XCes->Fill(eles.ele_XCes[iEleIndex], fWeight);
	hist->ZCes->Fill(eles.ele_ZCes[iEleIndex], fWeight);
	hist->HadEm->Fill(eles.ele_HadEm[iEleIndex], fWeight);
	hist->Chi2Mean->Fill(eles.ele_Chi2Mean[iEleIndex], fWeight);
	hist->N3d->Fill(eles.ele_N3d[iEleIndex], fWeight);
	hist->Iso4->Fill(eles.ele_Iso4[iEleIndex], fWeight);
	hist->EoverP->Fill(eles.ele_EoverP[iEleIndex], fWeight);
	hist->TrkPt->Fill(eles.ele_TrackPt[iEleIndex], fWeight);
	hist->TrkIso->Fill(eles.ele_TrkIso[iEleIndex], fWeight);
	hist->Ces2Wire->Fill(eles.ele_CesWireE2[iEleIndex], fWeight);
	hist->Ces2Strip->Fill(eles.ele_CesStripE2[iEleIndex], fWeight);
	hist->EmTime->Fill(eles.ele_EmTime[iEleIndex], fWeight);
	hist->PhiWedge->Fill(eles.ele_PhiWedge[iEleIndex], fWeight);
	hist->EmTimeVsRun->Fill(cVars.evt_RunNumber, eles.ele_EmTime[iEleIndex], fWeight);
	hist->EtCorrVsRun->Fill(cVars.evt_RunNumber, eles.ele_Etc[iEleIndex], fWeight);
	TLorentzVector tlEleVec(0,0,0,0);
	TVector2 tv2MetVec(cVars.met_MetX, cVars.met_MetY);
	tlEleVec.SetPxPyPzE(eles.ele_Px[iEleIndex], eles.ele_Py[iEleIndex], eles.ele_Pz[iEleIndex], eles.ele_E[iEleIndex]);

	//from V7 stuples and above has complete MET vector to do this
	const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlEleVec.Phi())));
	hist->PhoMetDelPhi->Fill(dphi_pm);
	hist->ClosestPhoMetDelPhi_a4->Fill(dphi_pm);
	
	hist->PhoEtMetProf->Fill(cVars.met_Met, eles.ele_Etc[iEleIndex], fWeight);
	hist->PhoEMetProf->Fill(cVars.met_Met, tlEleVec.E(), fWeight);
	hist->PhoEtaMetProf->Fill(cVars.met_Met, tlEleVec.Eta(), fWeight);
}

void WJetSel::PrintHeader(const CommonVars& cVars) const
{
	std::cout << " ============== " << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
}

void WJetSel::PrintJobSettings()
{	

	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "MC Flag setting        = ";
	if (bMcFlag) std::cout << "MC" << std::endl;
	else std::cout << "DATA" << std::endl;
	std::cout << "min,max Vtx for Signal = " << std::setw(4) << GetMinClass12Vtx() << ", " << std::setw(4) << GetMaxClass12Vtx() << std::endl;
	std::cout << "Ele Min/max Et, Eta    = " << std::setw(4) << GetMinEleEt() << ", " << std::setw(4) << GetMaxEleEt() << ", " << GetMaxEleEta()<< std::endl;
	std::cout << "Min/max MEt            = " << std::setw(4) << GetMinMet() << ", " << std::setw(4) << GetMaxMet()<< std::endl;
	std::cout << "Min W pt               = " << std::setw(4) << GetMinWPt() << std::endl;
	std::cout << "Jet Min Et/ Max Eta    = " << std::setw(4) << GetMinJetEt() << ", " << std::setw(4) << GetMaxJetEta()<< std::endl;
	std::cout << "Nvtx Weghting?         = " << std::setw(4) << UseNvtxWeights();
	if (UseNvtxWeights() == 0)	std::cout << " (0=NO)" << std::endl;
	else if (UseNvtxWeights() == 1)	std::cout << " (YES, Central Weights)" << std::endl;
	else if (UseNvtxWeights() == 2)	std::cout << " (YES, Central+1sigma Weights)" << std::endl;
	else std::cout << "UNKNOWN SETTING??" << std::endl;
}
void WJetSel::FillWHists(const CommonVars& cVars, const ElectronList& eles, 
						Histograms::TwoJetsHists_t* hist,
						const int iEle1Index,
						const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	const TLorentzVector tlEle(eles.ele_Px[iEle1Index], eles.ele_Py[iEle1Index], eles.ele_Pz[iEle1Index], eles.ele_E[iEle1Index]);
	const TLorentzVector tlMet(cVars.met_MetX,cVars.met_MetY,0,cVars.met_Met);
	//float DelEta = fabs(tlEle.DeltaEta(tlMet));
	
	const TLorentzVector tlW = tlEle + tlMet; 
	float fEtRatio = (tlMet.Pt()/tlEle.Pt());
	float DelPhi = fabs(tlEle.DeltaPhi(tlMet));
	float DelR = tlEle.DeltaR(tlMet);
	
	
	if (tlW.M()<60. ||tlW.M()>100.) {
		cVars.PrintHeader();
		//cVars.Dump();
		std::cout << "ele ind = " << iEle1Index << std::endl;
		eles.Dump();
		assert(false);
	}

	hist->InvMass->Fill(tlW.M(), fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	//hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(fEtRatio, fWeight);
	hist->Pt->Fill(tlW.Pt(), fWeight);
}

void WJetSel::FillWJetHists(const CommonVars& cVars, const ElectronList& eles, const JetList& jets,
						Histograms::TwoJetsHists_t* hist,
						const int iEle1Index,
						const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	const TLorentzVector tlEle(eles.ele_Px[iEle1Index], eles.ele_Py[iEle1Index], eles.ele_Pz[iEle1Index], eles.ele_E[iEle1Index]);
	const TLorentzVector tlMet(cVars.met_MetX,cVars.met_MetY,0,cVars.met_Met);
	const TLorentzVector tlWvec = tlEle + tlMet; 
	//float DelEta = fabs(tlEle.DeltaEta(tlMet));
	const int iJetIndex = 0;//leading jet
	const TLorentzVector tlJet(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);
	const TLorentzVector tlWjet = tlWvec + tlJet;
	
	float fEtRatio = (tlJet.Pt()/tlWvec.Pt());
	float DelPhi = fabs(tlWvec.DeltaPhi(tlJet));
	float DelR = tlWvec.DeltaR(tlJet);
	
	
	hist->InvMass->Fill(tlWjet.M(), fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	//hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(fEtRatio, fWeight);
	hist->Pt->Fill(tlWjet.Pt(), fWeight);
}

void WJetSel::FillJetHists(const CommonVars& cVars, const JetList& jets, Histograms::JetHists_t* hist, 
								const int iJetIndex,	const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	hist->Emfr->Fill(jets.jet_Emfr[iJetIndex], fWeight);
	hist->DetEta->Fill(jets.jet_DetEta[iJetIndex], fWeight);
	hist->NTracks->Fill(jets.jet_Ntracks[iJetIndex], fWeight);
	if (jets.jet_Ntowers[iJetIndex]==0)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ << "A jet with ntower==0 found!!! ";
		//cVars.PrintHeader();
		//jets.DumpAll();
	}
	hist->NTowers->Fill(jets.jet_Ntowers[iJetIndex], fWeight);
	hist->EtCorr->Fill(jets.jet_Pt[iJetIndex], fWeight);

	TLorentzVector tlVec(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);
	hist->EvtPhi->Fill(tlVec.Phi(), fWeight);
	hist->EvtEta->Fill(tlVec.Eta(), fWeight);
	//hist->SeedIEta->Fill(jets.jet_SeedIEta[iJetIndex], fWeight);
	//hist->SeedIPhi->Fill(jets.jet_SeedIPhi[iJetIndex], fWeight);
	
	TVector2 tv2MetVec(cVars.met_MetX, cVars.met_MetY);

	//from V7 stuples and above has complete MET vector to do this
	const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlVec.Phi())));
	hist->JetMetDelPhi->Fill(dphi_jm);
	hist->ClosestJetMetDelPhi_a4->Fill(dphi_jm);
	
	hist->JetEtMetProf->Fill(cVars.met_Met, tlVec.Pt(), fWeight);
	hist->JetEMetProf->Fill(cVars.met_Met, tlVec.E(), fWeight);
	hist->JetEtaMetProf->Fill(cVars.met_Met, tlVec.Eta(), fWeight);
}

//can't do this do not have the 2nd track info -- 08-16-2010
/*int WJetSel::IsSidebandPhoLikeElectron(const ElectronList& cEles, const int index)
{
	
	if (fabs(myStuple.pho_XCes[iPhoIndex]) > 21) return false;
	
	if (fabs(myStuple.pho_ZCes[iPhoIndex]) < 9 ||  fabs(myStuple.pho_ZCes[iPhoIndex]) > 230) return false;
	
	if (myStuple.pho_HadEm[iPhoIndex] > 0.125) return false;

	if (myStuple.pho_Etc[iPhoIndex] < 20) {
		if (myStuple.pho_Iso4[iPhoIndex] > 0.15*myStuple.pho_Etc[iPhoIndex]) return false;
	} else {
		if (myStuple.pho_Iso4[iPhoIndex] > 3.0 + 0.02 * (myStuple.pho_Etc[iPhoIndex] - 20.0)) return false;
	}


	if (myStuple.pho_N3d[iPhoIndex] > 2)


	//E/P of 1st track
	if (cEle.ele_N3d[index]>=1)
	{
		
	}
		

		if (myStuple.pho_TrkPt[iPhoIndex] > (0.25*myStuple.pho_Etc[iPhoIndex])	) {
			IDWord |= kTrkPtBit;
		}


		if (myStuple.pho_TrkIso[iPhoIndex] >5.0 ) {
			IDWord |= kTrkIsoBit;
		}
	}
	
	return IDWord;

}
*/
