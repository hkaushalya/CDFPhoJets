#include "TMath.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TBox.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TRandom.h"


/*
	.L TPeaksHiggs.C+
	Plot p;
	p.SetName("13Aug07");
// PRL plots
//p.SetPRL(1);

// mode=0 log, 1 lin, 2=int
// fitMode=0 for floating DiPhox contribution, 1 for fixing it
//this can take a while due to the trials thrown for 
//finding systematic uncertainty
p.FitSignal(1,0);


input file
gamgam_x_2fb.root

And they are titled:

pair_mass_2GeV1  (for CC)
pair_mass_2GeV3  (for CP)


shape experiments

TF1* xx = new TF1("xxf","[0]*pow((x-30.0),[1])*exp(-[2]*(x-30.0))",30.0,1000.0);
xx->SetParameter(0,0.005)
xx->SetParameter(1,0.5)
//xx->SetParLimits(1,1,-1.)
xx->SetParameter(2,0.03)
pair_mass_2GeV1->Fit(xx,"","",30.0,1000.0)

TF1* xx = new TF1("xxf","[0]*pow((x-30.0),[1])*(exp(-[2]*(x-30.0))+[3]*exp(-[4]*(x-30.0)))",0.0,500.0);
xx->SetParameter(0,0.005)
xx->SetParameter(1,2)
//xx->SetParLimits(1,1,-1.)
xx->SetParameter(2,0.03)
pair_mass_2GeV1->Fit(xx,"L","",30.0,500.0)

TF1* xx = new TF1("xxf","[0]*(pow((x-30.0),[1])+[3]*pow((x-30.0),0.5))*exp(-[2]*(x-30.0))",0.0,500.0);
xx->SetParameter(0,0.005)
xx->SetParameter(1,2)
//xx->SetParLimits(1,1,-1.)
xx->SetParameter(2,0.03)
xx->SetParameter(3,0.01)
pair_mass_2GeV1->Fit(xx,"L","",30.0,1000.0)

//possible winner
TF1* xx = new TF1("xxf","([0]*pow((x-30.0),[1])+[3]*pow((x-30.0),0.2))*exp(-[2]*(x-30.0))",30.0,500.0);
xx->SetParameter(0,0.005)
xx->SetParameter(1,2)
//xx->SetParLimits(1,1,-1.)
xx->SetParameter(2,0.03)
xx->SetParameter(3,0.005)
pair_mass_2GeV1->Fit(xx,"L","",31.0,1000.0)

//best
TF1* xx = new TF1("xxf","([0]*pow((x-30.0),[1])+[3]*pow((x-30.0),0.2))*(exp(-[2]*(x-30.0))+[4]*exp(-[5]*(x-30.0)))",30.0,500.0);
xx->SetParameter(0,0.005)
xx->SetParameter(1,2)
//xx->SetParLimits(1,1,-1.)
xx->SetParameter(2,0.02)
xx->SetParameter(3,0.005)
xx->SetParameter(4,1.005)
xx->SetParameter(5,0.05)
pair_mass_2GeV1->Fit(xx,"L","",31.0,1000.0)
pair_mass_2GeV3->Fit(xx,"L","",31.0,1000.0)

*/


class Plot {
	public:
		Plot();
		~Plot();
		TCanvas* NewCanvas();
		// use this to tag the name of the output eps files
		void SetName(char* n) { name=new char[100]; sprintf(name,"%s",n);}

		//mode 0=log, 1=lin
		void FitSignal(int mode, int fitMode=0);

		// for making PRL plots
		void SetPRL(int p = 0) { prl = p;}

		// internal use
		void savePlot(TCanvas* ccc, char* plot);

		TFile *fd,*fp;

		TCanvas* ccc[100];
		int nccc;
		char* name;

		int prl;
		bool debug;

};


Plot::Plot() {

	fd = new TFile("gamgam_x_2fb.root");    fd->ReadAll();

	name = NULL;
	nccc=0;
	prl = 0;
	debug = false;
}

Plot::~Plot() {
	fd->Close();
}


void Plot::FitSignal(int mode, int fitMode) {

	const int nPar = 6;
	TRandom ran;

	if(mode==0) {
		gStyle->SetOptLogy(1);
	} else {
		gStyle->SetOptLogy(0);
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	const float limitBinSize = 2.0; // **** bin size here
	TCanvas* c = NewCanvas();
	c->Divide(1,2);

	gROOT->cd();
	TH1F* cc = new TH1F("CCSignal","CC Signal",500,0.0,1000.0);
	TH1F* ccBg = new TH1F("CCBgFit","CC Bg Fit",500,0.0,1000.0);
	TH1F* ccBgErr = new TH1F("CCBgErr","CC Bg Err",500,0.0,1000.0);
	TH1F* ccBgP = new TH1F("CCBgPlus","CC Bg Plus",500,0.0,1000.0);
	TH1F* ccBgM = new TH1F("CCBgMinus","CC Bg Minus",500,0.0,1000.0);
	TH1F* cp = new TH1F("CPSignal","CP Signal",500,0.0,1000.0);
	TH1F* cpBg = new TH1F("CPBgFit","CP Bg Fit",500,0.0,1000.0);
	TH1F* cpBgErr = new TH1F("CPBgErr","CP Bg Err",500,0.0,1000.0);
	TH1F* cpBgP = new TH1F("CPBgPlus","CP Bg Plus",500,0.0,1000.0);
	TH1F* cpBgM = new TH1F("CPBgMinus","CP Bg Minus",500,0.0,1000.0);
	TMatrixD matrix(nPar,nPar);

	fd->cd();

	TH1F* hInt,*hBgInt;
	char fitname[100];
	for(int ind=0; ind<2; ind++) {
		if(debug) printf("starting ind %i\n",ind);
		c->cd(ind+1);
		gStyle->SetOptLogy(1);
		printf("Starting %i ######################################\n",ind);

		TH1F* h;
		//char cind[20];
		//char handle[100];
		//sprintf(handle,"side_1exp_%02i_%02i_%02i",ind,mode,fitMode);
		TF1* fits[4];
		//TF1* dpx[4];
		if(debug) printf("looking for h %i\n",int(fd));
		if(ind==0) {
			h = (TH1F*)fd->FindObjectAny("pair_mass_2GeV1");
		} else if(ind==1) {
			h = (TH1F*)fd->FindObjectAny("pair_mass_2GeV3");
		}
		if(debug) printf("new h %i\n",int(h));

		if(debug) printf("new fit\n");
		sprintf(fitname,"hfit_%1i",ind);
		fits[ind] = new TF1(fitname,"([0]*pow((x-30.0),[1])+[3]*pow((x-30.0),0.2))*([2]*exp(-[2]*(x-30.0))+[4]*[5]*exp(-[5]*(x-30.0)))",30.0,500.0);
		//fits[ind] = new TF1(fitname,"([0]*((1-[3])*pow((x-30.0),[1])+[3]*pow((x-30.0),0.2)))*(exp(-[2]*(x-30.0))+[4]*exp(-[5]*(x-30.0)))",30.0,500.0);
		fits[ind]->SetParameter(0,0.0004);
		fits[ind]->SetParameter(1,2);
		fits[ind]->SetParameter(2,0.02);
		fits[ind]->SetParameter(3,0.005);
		//fits[ind]->SetParameter(3,0.5);
		fits[ind]->SetParameter(4,1.005);
		fits[ind]->SetParameter(5,0.05);

		float llim = 30.0;
		h->Fit(fits[ind],"LN","",llim,1000.0);

		double par[20],parMin[20],fval,fvalMin;
		for(int i=0; i<nPar; i++) parMin[i] = fits[ind]->GetParameter(i);

		gMinuit->Eval(nPar,0,fvalMin,parMin,0);
		//printf("got back %10.5f\n",fvalMin);

		// save the fit results in a histogram, for limit program
		for(int ibin=16; ibin<250; ibin++) {
			float xx = h->GetBinCenter(ibin);
			float yy = fits[ind]->Eval(xx);
			if(ind==0) {
				cc->SetBinContent(ibin,h->GetBinContent(ibin));
				ccBg->SetBinContent(ibin,yy);
				ccBgErr->SetBinContent(ibin,0.0);
				ccBgP->SetBinContent(ibin,0.0);
				ccBgM->SetBinContent(ibin,99999.0);
			} else {
				cp->SetBinContent(ibin,h->GetBinContent(ibin));
				cpBg->SetBinContent(ibin,yy);
				cpBgErr->SetBinContent(ibin,0.0);
				cpBgP->SetBinContent(ibin,0.0);
				cpBgM->SetBinContent(ibin,99999.0);
			}
		}

		//vary the parameters to find an error envelope
		double par2[20],fval2=1e10;
		int pslim = (ind==0?25000:150000);
		for(int ips=0; ips<pslim; ips++) {
			if(ips%10000==0) printf("Processing %d\n",ips);
			for(int i=0; i<nPar; i++) {
				par[i] = parMin[i];
			}
			for(int i=0; i<nPar; i++) {
				//int i = (ips%2==0?0:3);
				par[i] = parMin[i]+(2.0*(ran.Uniform()-0.5))*fits[ind]->GetParError(i);
			}
			fval = 0.0;
			gMinuit->Eval(nPar,0,fval,par,0);
			if((fval-fvalMin)<1.0) {
				printf("Found nearby min %10.5f\n",fval-fvalMin);
				float eOld,eNew;
				for(int ibin=16; ibin<250; ibin++) {
					float xx = h->GetBinCenter(ibin);
					for(int i=0; i<nPar; i++) fits[ind]->SetParameter(i,par[i]);
					float yy = fits[ind]->Eval(xx);
					for(int i=0; i<nPar; i++) fits[ind]->SetParameter(i,parMin[i]);
					float yyMin = fits[ind]->Eval(xx);
					TH1F *hBgErr,*hBgP,*hBgM;
					if(ind==0) {
						hBgErr = ccBgErr; hBgP = ccBgP; hBgM = ccBgM;
					} else {
						hBgErr = cpBgErr; hBgP = cpBgP; hBgM = cpBgM;
					}

					eOld = hBgErr->GetBinContent(ibin);
					eNew = yy - yyMin;
					if(eOld>fabs(eNew)) hBgErr->SetBinContent(ibin,fabs(eNew));
					eOld = hBgP->GetBinContent(ibin);
					if(yy>eOld)  hBgP->SetBinContent(ibin,yy);
					eOld = hBgM->GetBinContent(ibin);
					if(yy<eOld)  hBgM->SetBinContent(ibin,yy);
				}

			} // end if near maximum

			/*
				if(fval<fval2) {
				for(int i=0; i<nPar; i++) par2[i] = par[i];
				fval2 = fval;
				}
				*/
		}

		/*
			printf("forcing new fit..\n");
			for(int i=0; i<nPar; i++) {
			printf("old,new = %10.5f %10.5f\n",parMin[i],par2[i]);
			fits[ind]->SetParameter(i,par2[i]);
			}
			*/

		// restore original fit
		fval = 0.0;
		gMinuit->Eval(nPar,0,fval,parMin,0);
		for(int i=0; i<nPar; i++) fits[ind]->SetParameter(i,parMin[i]);



		//extract fit error matrix
		gMinuit->mnemat(matrix.GetMatrixArray(),nPar);
		matrix.Print();

		for(int i=0; i<nPar; i++) {
			for(int j=0; j<nPar; j++) {
				printf("%10.5f",matrix(i,j)/sqrt(matrix(i,i)*matrix(j,j)));
			}
			printf("\n");
		}
		//matrix.Draw("text");

		float hm = h->GetMaximum();
		if(mode==0) {
			//TAxis* ax = h->GetXaxis();
			//ax->SetRangeUser(24.1,199.9);
			h->SetMaximum(1.2*hm);
			//h->SetMinimum(0.0);
		} else if(mode==1) {
			TAxis* ax = h->GetXaxis();
			ax->SetRangeUser(20.0,500.0);
			h->SetMaximum(1.15*hm);
			h->SetMinimum(0.0);
		}


		h->Draw();
		fits[ind]->SetLineColor(1);
		fits[ind]->SetLineWidth(2.0);
		fits[ind]->Draw("SAME");
		// find chi2's and KS's
		//AnaChiKs(h,fits[ind]);



		TAxis* ax,*ay;
		ax = h->GetXaxis(); 
		ay = h->GetYaxis();
		ax->SetTitle("m(#gamma#gamma) (GeV/c^{2})"); 
		ay->SetTitle("Entries/2 GeV/c^{2}");
		ax->CenterTitle(); ay->CenterTitle();
		ax->SetTitleOffset(0.9);
		ay->SetTitleOffset(1.0);
		ax->SetTitleSize(0.08);
		ay->SetTitleSize(0.07);
		ax->SetLabelSize(0.07);
		ay->SetLabelSize(0.07);

		gPad->SetLeftMargin(0.16);
		gPad->SetBottomMargin(0.16);

		TText* text;
		text = new TLatex(0.5,0.8,"Diphoton Data");
		text->SetNDC(true);
		text->SetTextSize(0.06);
		text->Draw();
		if(ind==0)      text = new TLatex(0.5,0.72,"Central-Central");
		else if(ind==1) text = new TLatex(0.5,0.72,"Central-Plug");
		text->SetNDC(true);
		text->SetTextSize(0.06);
		text->Draw();
		if(ind==0) {
			text = new TLatex(0.15,0.92,"W/Z H#rightarrow X(#gamma#gamma)");
			text->SetNDC(true);
			text->SetTextSize(0.08);
			text->Draw();

			text = new TLatex(0.5,0.92,"CDF Run II Preliminary, 2.0 fb^{-1}");
			text->SetNDC(true);
			text->SetTextSize(0.06);
			text->Draw();
		}    

		/*
			if(debug) printf("start loop\n");
			int ibin;
			for(ibin=16; ibin<=250; ibin++) {
			if(debug) printf("start bin            %i\n",ibin);
			float xx = (ibin-0.5)*2.0; // *** bin width here
			if(debug) printf("-1 test ibin %i\n",ibin);
			float yy = fits[ind]->Eval(xx);
		//printf("%f  yy= %f \n",xx,yy);
		// the derivative of this yield wrt parameters
		if(debug) printf("0 test ibin %i\n",ibin);
		double y0 = yy;
		if(debug) printf("1 test ibin %i\n",ibin);
		TMatrixD vv(nPar,1);
		float dirSize = 0.5;
		for(int i=0; i<nPar; i++){
		int ipar = i;
		double par = fits[ind]->GetParameter(ipar);
		double spar = fits[ind]->GetParError(ipar);
		double parp = par + dirSize*spar;
		fits[ind]->SetParameter(ipar,parp);
		double yp = fits[ind]->Eval(xx);
		vv(i,0) = limitBinSize*(yp-y0)/(dirSize*spar);
		fits[ind]->SetParameter(ipar,par);
		//printf("%f %f %f\n",yp,y0,spar);
		}
		//vv.Print();
		if(debug) printf("start matrix %i\n",ibin);
		TMatrixD tempM(matrix, TMatrixDBase::kMult, vv);
		//matrix.Print();
		TMatrixD tempN(vv, TMatrixDBase::kTransposeMult, tempM);
		//tempN.Print();
		float bgSig = 0.0;
		if(tempN(0,0)>0.0) bgSig = sqrt(tempN(0,0));
		// ****** hack temp  **********
		bgSig = 0.3*y0;


		// file hists to be saved
		if(debug) printf("start fill %i\n",ibin);
		if(ind==0) {
		//printf("filling cc %i %f\n",ibin,h->GetBinContent(ibin));
		cc->SetBinContent(ibin,h->GetBinContent(ibin));
		//printf("getting cc %i %f\n",ibin,cc->GetBinContent(ibin));
		ccBg->SetBinContent(ibin,yy);
		ccBgErr->SetBinContent(ibin,bgSig);
		ccBgP->SetBinContent(ibin,yy+bgSig);
		ccBgM->SetBinContent(ibin,TMath::Max(yy-bgSig,float(0.0)));
		//if(ibin==27) {
		//printf("bg %f %f \n",yy,bgSig);
		//}
		} else {
		cp->SetBinContent(ibin,h->GetBinContent(ibin));
		cpBg->SetBinContent(ibin,yy);
		cpBgErr->SetBinContent(ibin,bgSig);
		cpBgP->SetBinContent(ibin,yy+bgSig);
		cpBgM->SetBinContent(ibin,TMath::Max(yy-bgSig,float(0.0)));
		}
		if(debug) printf("end fill %i\n",ibin);
		}
		*/

	}

	printf("cc plus  BG=%f\n",ccBgP->GetSum());
	printf("cc minus BG=%f\n",ccBgM->GetSum());
	printf("cp plus  BG=%f\n",cpBgP->GetSum());
	printf("cp minus BG=%f\n",cpBgM->GetSum());

	char fn[100];
	if(mode==0) {
		sprintf(fn,"FitSignal_%d",fitMode);
		savePlot(c,fn);
	} else if(mode==1) {
		sprintf(fn,"FitSignalLin_%d",fitMode);
		savePlot(c,fn);
	}

	//if(mode!=0) return;

	// plot of fit results
	gStyle->SetOptLogy(0);
	c = NewCanvas();
	c->Divide(1,2);

	c->cd(1);
	cc->Draw();
	ccBg->Draw("SAME");
	c->cd(2);
	ccBgErr->SetMinimum(0.0); ccBgErr->SetMaximum(4.0); 
	ccBgErr->Draw();
	ccBgP->SetLineStyle(2); ccBgP->Draw("SAME");
	ccBgM->SetLineStyle(2); ccBgM->Draw("SAME");

	savePlot(c,"FitSignalResultsCC");

	c = NewCanvas();
	c->Divide(1,2);

	c->cd(1);
	cp->Draw();
	cpBg->Draw("SAME");
	c->cd(2);
	cpBgErr->SetMinimum(0.0); cpBgErr->SetMaximum(4.0); 
	cpBgErr->Draw();
	cpBgP->SetLineStyle(2); cpBgP->Draw("SAME");
	cpBgM->SetLineStyle(2); cpBgM->Draw("SAME");

	savePlot(c,"FitSignalResultsCP");

	char title[100];
	if(name) {
		sprintf(title,"TPeaksHiggs_FitSignalHist_%s.root",name);
		TFile* ff = new TFile(title,"RECREATE");
		gROOT->GetList()->Write();
		ff->Close();
	}

}


/*****************************************/

void Plot::savePlot(TCanvas* ccc, char* plot) {
	char title[100];
	if(name) {
		sprintf(title,"TPeaksHiggs_%s_%s.eps",plot,name);
		ccc->Print(title);
		//sprintf(title,"/cdf/scratch/rlc/TDiPho_%s_%s.gif",plot,name);
		//ccc->Print(title);
	}

}

/*****************************************/

TCanvas* Plot::NewCanvas() {
	char title[100];
	sprintf(title,"ccc%d",nccc);
	TCanvas* c = new TCanvas(title,title,100+20*nccc,100+20*nccc,600,600);
	ccc[nccc] = c;
	nccc++;
	return c;
}
