#include<iostream>
#include<string>
#include<TH1.h>
#include<TF1.h>
#include<TStyle.h>
#include<TCanvas.h>
#include<TROOT.h>
#include<TFile.h>
#include<TFrame.h>
#include<sstream>

//=====================================================================
//-------------- fit function: p0*exp(-p1*x+p2*x*x)
double Func1(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0)
	{
		//fitval=par[0]*exp(-1.0*par[1]*x[0]+par[2]*x[0]*x[0])+par[3]*exp(-1.0*par[4]*x[0]);
		fitval=par[0]*exp(-1.0*par[1]*x[0]+par[2]*x[0]*x[0])+par[3]*exp(-1.0*par[4]*x[0]);

		//fitval=(par[0]+par[1]*x[0])*exp(par[2]*x[0])+(par[3]+par[4]*x[0])*exp(par[5]*x[0]);
		// p0*exp(-p1*x+p2*x*x)+p3*exp(-1*p4*x);
	}
	else fitval=-1.0E6;
	return fitval;
	
}

double CorrFun(double *x, double *par)
{
	float corr=1.0;
	float cut1=0.0;
	float cut2=20.0;
	TF1 *fitFun=new TF1("fit",Func1,cut1,cut2,5);
	fitFun->SetParameter(0,par[0]);
	fitFun->SetParameter(1,par[1]);
	fitFun->SetParameter(2,par[2]);
	fitFun->SetParameter(3,par[3]);
	fitFun->SetParameter(4,par[4]);

	//float delta=(cut2-cut1)/1000.0;
	float bin_int=0.0;
	float bin_int0=0.0;
	float last_bin_slope=1.0;
	float tot_int=fitFun->Integral(0.0,20.0);

	bin_int=fitFun->Integral(0.0,x[0]);
	bin_int=bin_int/tot_int;
	if(bin_int>0.0 && bin_int<1.0 && x[0]<par[5]) corr=-log10(1-bin_int);
	else 
	{
		bin_int0=fitFun->Integral(0.0,par[5]);
		bin_int0=bin_int0/tot_int;
		last_bin_slope=-log10(1-bin_int0)/par[5];      
		corr=x[0]*last_bin_slope;
	}
	delete fitFun;
	return corr;
}


//=====================================================================
//-------------- Main script ------------------------------------------
// titleCode= 0--data signal; 1--sideband; 2-- PhotonMc 3-- Z->ee data; 
// njet: -1==all, 0==for Njet15=0, 1== for Njet15=1; 2== for Njet15=2 
// 3==for Njet15=3, 4== for Njet15>=4 
// UnclPara = 0 (gg sideband), 1 (Zee)
void MetSigFitter(const std::string fName, const int titleCode=0,
						const int UnclPara=0, const int njet=-1,
						const int rebincode=1, float cut2 = 14,
						TH1* H1 = 0, TH1* H2=0)
{

	assert (fName.length() > 0 && "file name not specified");
	assert (titleCode >=0 && "titleCode cannot be negative!");
	assert (( (UnclPara == 0) || (UnclPara == 1)) && "Uncl. Paramete is either 0,1 !");
	assert (njet>= -1 && "njet cannot be <-1!");
	assert (rebincode >0 && "rebin size cannot be < 1!");

	TFile* myFile= new TFile(fName.c_str()); //connect file, everything stays in memory till "delete"
	if (myFile->IsZombie())
	{
		std::cout << "File " << fName << " did not open! exiting" <<std::endl;
		exit(1);
	}


	gROOT->Reset();
	gStyle->SetOptFit(1111);
	TH1::AddDirectory(kFALSE);

	unsigned int ww = 850, hh=670;
	TCanvas *c1 = new TCanvas("c1","Met Model Study",200,10,700,500);
	//TCanvas *c1 = new TCanvas();
	c1->GetFrame()->SetBorderSize(12);
	c1->SetCanvasSize(ww-20,hh-30);
	c1->SetWindowSize(ww, hh);
	c1->Divide(1,2);

	std::string sUnclPara;
	if (UnclPara == 0) sUnclPara = "0 (#gamma#gamma sideband)";
	else if (UnclPara == 1) sUnclPara = "1 (Zee MC)";


	std::string fObject;


	TH1F *h_MetSig_calib;
	if(njet==-1) fObject = "fMetSig_estimate"; // raw met sig 
	if(njet==0) fObject = "fMetSigCalib_estimate_njet0"; // get object name  
	if(njet==1) fObject = "fMetSigCalib_estimate_njet1"; // get object name  
	if(njet==2) fObject = "fMetSigCalib_estimate_njet2"; // get object name  
	if(njet==3) fObject = "fMetSigCalib_estimate_njet3"; // get object name  
	if(njet==4) fObject = "fMetSigCalib_estimate_njet4"; // get object name  

	std::string snjet;
	if (njet == -1) snjet = "Njet>=0 (Raw #slash{E}_{T}-sig)";
	if (njet == 0) snjet = "Njet=0";
	if (njet == 1) snjet = "Njet=1";
	if (njet == 2) snjet = "Njet=2";
	if (njet == 3) snjet = "Njet=3";
	if (njet == 4) snjet = "Njet>=4";


	h_MetSig_calib= (TH1F*) gDirectory->FindObjectAny(fObject.c_str());   // get histogram
	assert(h_MetSig_calib != NULL && "h_MetSig_calib is null");

	//* I need to do this in order to clean memory from ROOT prompt
	TFile* mytmpfile= new TFile("tmp.tmp","recreate"); // create auxiliary file
	mytmpfile->cd();

	//-------------------------------------- Met histograms
	int nbins = h_MetSig_calib->GetNbinsX();
	int Nevnt = (int) h_MetSig_calib->GetEntries();
	float binwidth = h_MetSig_calib->GetBinWidth(1);
	float lowedge  = (h_MetSig_calib->GetBinCenter(1)) - binwidth / 2.0;
	float highedge = lowedge + nbins * binwidth;
	float mean     = h_MetSig_calib->GetMean();

	if(mean<=0.0) mean=0.5;			//sam: not letting it go below 0. why?

	// setup titles
	std::stringstream title;
	std::string name("Met");
	if(titleCode==0) title << "Estimated upper limit on #slash{E}_{T}-significance in #gamma+jets: data signal" << ", " << snjet << ", Uncl Para=" << sUnclPara;
	if(titleCode==1) title << "Estimated upper limit on #slash{E}_{T}-significance in #gamma+jets: data sideband" << ", " << snjet << ", Uncl Para=" << sUnclPara;
	if(titleCode==2) title << "Estimated upper limit on #slash{E}_{T}-significance in #gamma+jets: #gamma MC signal" << ", " << snjet << ", Uncl Para=" << sUnclPara;
	if(titleCode==3) title << "Estimated upper limit on #slash{E}_{T}-significance in Z#rightarrowe^{+}e^{-}" << ", " << snjet << ", Uncl Para=" << sUnclPara;
	if(titleCode==4) title << "Estimated upper limit on #slash{E}_{T}-significance in W#rightarrowe+#nu MC signal" << ", " << snjet << ", Uncl Para=" << sUnclPara;
	if(titleCode==5) title << "Estimated upper limit on #slash{E}_{T}-significance in W+jets#rightarrow e+jets (#gamma DATA)" << ", " << snjet << ", Uncl Para=" << sUnclPara;
	if(titleCode==6) title << "Estimated upper limit on #slash{E}_{T}-significance in W+jets#rightarrow e+jets (We#nu Inc. MC)" << ", " << snjet << ", Uncl Para=" << sUnclPara;

	TH1F* myFnlMet_gn= new TH1F(name.c_str(),title.str().c_str(),nbins,lowedge,highedge);
	myFnlMet_gn->Sumw2();

	//normalizing the hist to unity
	// why not just use scale(1/Nentries);
	// cut2 is the bin center of the last bin with contents 
	float bin_value;
	float bin_error;
	float cut1 = 0.0;  	//this is the low fit values
	//float cut2 = 0.0;		// this is the upper fit value
	float event_count = 0.0;

	for (int i=0; i<nbins; i++)
	{
		bin_value = h_MetSig_calib->GetBinContent(i+1);
		bin_error = h_MetSig_calib->GetBinError(i+1);
		bin_error = sqrt(bin_value) / (1.0 * Nevnt);
		event_count = event_count + bin_value;
		// I replaced this. Now the cut2 is an input value to the function.
		//if (bin_value > 0.0 && event_count <= 1.0 * Nevnt) cut2 = h_MetSig_calib->GetBinCenter(i+1); 
		myFnlMet_gn->SetBinContent(i+1, bin_value/(1.0 * Nevnt));
		myFnlMet_gn->SetBinError(i+1, bin_error);
	}

	if (rebincode>1) myFnlMet_gn->Rebin(rebincode);

	//this is where I can set the fit function parameters!
	TF1 *fitFun1=new TF1("fit1",Func1,0,20.0,5);
	////////////// parameters used to make first set of calib fit plots ////////
	
	/*	fitFun1->SetParameter(0,0.03 * rebincode);  //default value
	//fitFun1->SetParameter(1,1.478); // default value, controls the curvature around 4
	fitFun1->SetParameter(1,3.0);
	fitFun1->SetParameter(2,-0.3052);  //default value
	//fitFun1->SetParameter(2,-3.0);
	fitFun1->SetParameter(3,0.0014); 	//default value
	//fitFun1->SetParameter(3,0.02);
	fitFun1->SetParameter(4,0.9129);		// default value
	 */
	/////////////////////////////////////////////////////////////////////////////////

	fitFun1->SetParameter(0,0.03 * rebincode);
	//fitFun1->SetParameter(1,2.991);// rebin=15 till=12 is the best for jet=4
	fitFun1->SetParameter(1,10.051);

	fitFun1->SetParameter(1,1.991);
	fitFun1->SetParameter(2,-0.009);
	fitFun1->SetParameter(3,0.0005);
	fitFun1->SetParameter(4,0.0001);
	//fitFun1->SetParameter(5,0.0000);
	
	//fitFun1->SetParLimits(1,2.5,3.0);
	//fitFun1->SetParameter(5,1);
	//fitFun1->SetParameter(0,-1.23390e+00);
	//fitFun1->SetParameter(1,4.82919e-01);
   //fitFun1->SetParameter(2,-1.23390e+00/2);
	//fitFun1->SetParameter(3,0);
	//fitFun1->SetParLimits(1,0,1);
	//fitFun1->SetParLimits(3,0,1);


	c1->cd(1);
	gPad->SetLogy();
	myFnlMet_gn->SetMarkerStyle(20);
	myFnlMet_gn->SetMarkerColor(4);
	myFnlMet_gn->SetTitle(title.str().c_str());
	myFnlMet_gn->GetXaxis()->SetTitle("-log(1-P_{#slash{E}_{T}^{fluc}<#slash{E}_{T}^{meas}})");
	std::stringstream ytitle;
	ytitle << "Events / " << myFnlMet_gn->GetBinWidth(1);
	myFnlMet_gn->GetYaxis()->SetTitle(ytitle.str().c_str());
	myFnlMet_gn->GetXaxis()->SetLabelSize(0.04);
	myFnlMet_gn->GetYaxis()->SetLabelSize(0.04);
	myFnlMet_gn->SetFillColor(4);
	myFnlMet_gn->Fit(fitFun1,"E","",cut1,cut2);  //"E"=perform better error estimate using Minos // 
	myFnlMet_gn->SetDirectory(0);
	myFnlMet_gn->Draw("ep");

	TF1* fitted = myFnlMet_gn->GetFunction("fit1");
	//for (unsigned iter = 0; iter != 20; ++ iter)
	//	std::cout << "eval(" << iter << ")=" << fitted->Eval (iter) << std::endl;
	assert (fitted != NULL && "fitted function returned NULL!");
	double chi2 = fitted->GetChisquare();
	double ndf = fitted->GetNDF();
	double  freepar = fitted->GetNumberFreeParameters();
	std::cout << "********************************************************" << std::endl;
	std::cout << "*** CHI2/NDF = " << chi2 << "/" << ndf << " = " << chi2/ndf << "\tNFreePara = " << freepar << std::endl;
	std::cout << "********************************************************" << std::endl;
	

	H1 = (TH1*) myFnlMet_gn->Clone("met_sig_copy");
	H1->SetDirectory(0);
	assert (H1 != NULL);

	//std::cout<<"____ cut2="<<cut2<<std::endl;

	c1->cd(2);
	TF1 *corFun1 = new TF1("corr",CorrFun,0.0,20.0,6);
	corFun1->SetParameter(0,fitFun1->GetParameter(0));
	corFun1->SetParameter(1,fitFun1->GetParameter(1));
	corFun1->SetParameter(2,fitFun1->GetParameter(2));
	corFun1->SetParameter(3,fitFun1->GetParameter(3));
	corFun1->SetParameter(4,fitFun1->GetParameter(4));
	corFun1->SetParameter(5,cut2);
	corFun1->GetYaxis()->SetTitle("CorrFun");
	corFun1->Draw();

	//H2 = (TH1*)corFun1->Clone("corFun1_copy");
	//H2->SetDirectory(0);

	c1->cd();
	std::stringstream gifName;
	if (titleCode == 0) gifName << "PhoDataSignal" << fObject <<"_Unlc"<< UnclPara << ".gif"; 
	if (titleCode == 1) gifName << "PhoDataSideband" << fObject <<"_Unlc"<< UnclPara << ".gif"; 
	if (titleCode == 2) gifName << "PhoMCSignal" << fObject <<"_Unlc"<< UnclPara << ".gif"; 
	if (titleCode == 3) gifName << "Zee" << fObject <<"_Unlc"<< UnclPara << ".gif"; 
	if (titleCode == 4) gifName << "WenMCSignal" << fObject <<"_Unlc"<< UnclPara << ".gif"; 
	if (titleCode == 5) gifName << "PhoData_ejets" << fObject <<"_Unlc"<< UnclPara << ".gif"; 
	if (titleCode == 6) gifName << "WenMC_ejets" << fObject <<"_Unlc"<< UnclPara << ".gif"; 
	c1->Print(gifName.str().c_str());


	//	delete myFile;
	// in ROOT prompt do > delete gDirectory->GetFile()  

	/*new TCanvas();
	  H1->Draw();
	  new TCanvas();
	  H2->Draw();
	 */
}

void MetSigFitter(const int titlecode, const int njet, const int rebin, const float cut2)
{
	TH1 *H1, *H2;
	//MetSigFitter("Merged.root",2,0,1,2,8, H1, H2);
	//MetSigFitter("Merged.root",titlecode,0,njet,rebin,cut2, H1, H2);
	MetSigFitter("Merged_MetSigCalibHists.root",titlecode,0,njet,rebin,cut2, H1, H2);
	//assert (H1 != NULL && "H1 is empty");
	//assert (H2 != NULL && "H2 is empty");
	//new TCanvas();
	//H1->Draw();
	//new TCanvas();
	//H2->Draw();
}
