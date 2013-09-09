////////////////////////////////////////////////////////////
// This is my toy MC to study JET smearing  according to  //
// JER. I started with a monocromatic jet.                //
////////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/*
 * $Log: JetSmearStudy.cc,v $
 * Revision 1.2  2009/07/21 16:54:28  samantha
 * MINOR EDITS.
 *
 * Revision 1.1  2009/06/29 20:34:47  samantha
 * This is a toy MC to study the jet energy smearing accroding to JER fucntions.
 * I am starting out with a monochromatic jet and standard gausing/landau JER
 * functions.
 *
 *
 */


#include "samantha/Pho/JetSmearStudy.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "TF1.h"
#include "TLorentzVector.h"
#include "Stntuple/obj/TStnJet.hh"
#include "TRandom3.h"
#include "samantha/utils/FreeFunctions.hh"
#include "TCanvas.h"

ClassImp(JetSmearStudy)

//-------- Jet Energy Resolution: centered at zero
// //-------------- fit function: [Const*Gaus(-x/(1+x))+Landau(-x/(1+x))]/(1+Const)
double MetModelJER(double* x, double* par) {
	double value=0.0;
	double arg=-x[0]/(x[0]+1.0);
	double arg_L=-x[0]/(x[0]+1.0);
	double mean=par[0];
	double sigmaG=par[1];
	double mpv=par[2];
	double sigmaL=par[3];
	double normL=par[4];
	if(normL<0.0) normL=0.0;
	double f1=normL/(1.0+normL)*TMath::Gaus(arg,mean,sigmaG);
	double f2=TMath::Landau(arg_L,mpv,sigmaL)/(1.0+normL);
	value=f1+f2;
	if(value<0.0) return 0.0;
	return value;
}


//_____________________________________________________________________________
void JetSmearStudy::BookHistograms()
{
	DeleteHistograms();

	TFolder *new_folder;// *sub_folder;

	//new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	//vCounterHistLabels.push_back("Events Processed");				// 0
	//vCounterHistLabels.push_back("Events Passed");					// 1
	//vCounterHistLabels.push_back("All Tight #gamma Events");		// 2
	//HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	
	new_folder = GetHistFolder(this, "SmearStuff","SmearStuff");
	hJetPt_Orig = new TH1F("JetPtOrig","Original Jet Pt",1000,0,100);
	hJetPt_Orig->SetMarkerColor(kBlue);
	hJetPt_Orig->SetLineColor(kBlue);
	hJetPt_Orig->SetMarkerSize(4);
	hJetPt_Orig->SetLineWidth(2);
	hJetPt_Smear = new TH1F("JetPtSmear","Smeared Jet Pt",1000,0,100);
	hJetPt_Smear->SetMarkerColor(kRed);
	hJetPt_Smear->SetMarkerSize(4);
	hJetPt_Smear->SetLineColor(kRed);
	hJetPt_Smear->SetLineWidth(2);

	new_folder->Add(hJetPt_Orig);
	new_folder->Add(hJetPt_Smear);
}
	


//-------------------------------------------------------------------
int JetSmearStudy::BeginJob()
{
	//_____________________________________________________ register the data block
	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlock);
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();
	
	return 0;
}


//-------------------------------------------------------------------
int JetSmearStudy::BeginRun()
{
	return 0;
}


//-------------------------------------------------------------------
int JetSmearStudy::Event   (int ientry)
{

//	fJetBlock->GetEntry(ientry);
//	std::cout << "i\tPx\tPy\tPz\tE\tEta" << std::endl;

	
//	for (int i=0; i<fJetBlock->NJets(); ++i)
//	{
//		TStnJet* jet = fJetBlock->Jet(i);
//		TLorentzVector jetvec = *(jet->Momentum()); 
		//std::cout << i << "\t" << jetvec.Px()<< "\t" << jetvec.Py()<< "\t" << jetvec.Pz()<< "\t" 
		//				<< jetvec.E() <<"\t" << jetvec.Eta() << std::endl;
//	}

	TLorentzVector jetvec(28.4494,-9.44013,-79.3183,84.8118);
	hJetPt_Orig->Fill(jetvec.Pt());
	SmearJet(jetvec);
	
	return 0;
}
//-------------------------------------------------------------------
int JetSmearStudy::EndJob()
{

	std::cout << "----- EndJob :: " << GetName() <<std::endl;
	new TCanvas();
	hJetPt_Smear->Draw("E");
	new TCanvas();
	hJetPt_Orig->Draw("E");

	return 0;
}


void JetSmearStudy::SmearJet(const TLorentzVector jetvec)
{

	TRandom3 *rn = new TRandom3();		//should not use gRandom->GetSeed() twice as it returns the same seed.
	int rnd_seed_1 = rn->GetSeed();
	gRandom->SetSeed(rnd_seed_1);


	double parJER[5];
	int systcode = 0;
	parJER[0]=Jer_meanG(jetvec.E(),jetvec.Eta(),0,systcode); // for now systcode is used as jer_stat_code
	parJER[1]=Jer_sigmaG(jetvec.E(),jetvec.Eta(),0,systcode);
	parJER[2]=Jer_mpvL(jetvec.E(),jetvec.Eta(),0,systcode);
	parJER[3]=Jer_sigmaL(jetvec.E(),jetvec.Eta(),0,systcode);
	parJER[4]=Jer_normG(jetvec.E(),jetvec.Eta(),0,systcode);
	double jer_limit=IntegralUpLimit(jetvec.E(),jetvec.Eta(),0,systcode);
	TF1* JerFun=new TF1("jer",MetModelJER,-1.0,jer_limit,5);
	//------ setting Gaussian
	JerFun->SetParameter(0,parJER[0]);
	JerFun->SetParameter(1,parJER[1]);
	//------ setting Landau
	JerFun->SetParameter(2,parJER[2]);
	JerFun->SetParameter(3,parJER[3]);
	//------ setting Gaussian normalization
	JerFun->SetParameter(4,parJER[4]);
	//double off_set=JerFun->GetMaximumX(-1.0,1.0); // returns mpv of JER function;

	for (unsigned int i = 0; i < GetNptsToGenerate(); ++i)
	{
		//double jer_rand = JerFun->GetRandom(); // get jet energy scale factor; default=assume jet scale is 1.0


		//float scale_factor=1.0-off_set+jer_rand; // get jet energy scale factor; default=assume jet scale is 1.0
		//float scale_factor = rn->Gaus();		//mean=0, sigma=1
		float scale_factor = rn->Landau(0,1);  	//mean=0, sigma=1

		if(scale_factor<0.0) {
			scale_factor=0.0;
		}
		TLorentzVector newjetvec;
		newjetvec = scale_factor * jetvec;
		//std::cout << "Jet Pt Rat new/old= " << newjetvec.Pt()/jetvec.Pt() << std::endl; 
		//std::cout << "Jet Pt new/old= " << newjetvec.Pt() << "\t" << jetvec.Pt() << std::endl; 
		//hJet:q
		hJetPt_Smear->Fill(newjetvec.Pt());

	}

	delete JerFun;


}


//________________ returns JER param or its uncertainty
double JetSmearStudy::JERparam(int code, int i, int j)
{
	double jer_param[15][13]; // [eta][param]
	double jer_param_er[15][13]; // [eta][param]
	double value=0.0;

	jer_param[0][0]=-0.00505533;  jer_param_er[0][0]=0.00366313;
	jer_param[0][1]=-7.37332E-5;  jer_param_er[0][1]=1.02366E-5;
	jer_param[0][2]=5.06548;  jer_param_er[0][2]=0.198112;
	jer_param[0][3]=1.03604;  jer_param_er[0][3]=0.0162871;
	jer_param[0][4]=0.000363984;  jer_param_er[0][4]=9.99763E-5;
	jer_param[0][5]=0.0404354;  jer_param_er[0][5]=0.00394954;
	jer_param[0][6]=7.70406E-5;  jer_param_er[0][6]=1.80066E-5;
	jer_param[0][7]=-1.05943;  jer_param_er[0][7]=0.148324;
	jer_param[0][8]=0.157408;  jer_param_er[0][8]=0.00948233;
	jer_param[0][9]=0.00745406;  jer_param_er[0][9]=0.000132492;
	jer_param[0][10]=-4.48426;  jer_param_er[0][10]=0.997815;
	jer_param[0][11]=1.14401;  jer_param_er[0][11]=0.269919;
	jer_param[0][12]=0.0459092;  jer_param_er[0][12]=0.0168848;
	jer_param[1][0]=0.0928237;  jer_param_er[1][0]=0.00248223;
	jer_param[1][1]=-0.000201847;  jer_param_er[1][1]=8.82907E-6;
	jer_param[1][2]=0.177168;  jer_param_er[1][2]=0.137301;
	jer_param[1][3]=0.670063;  jer_param_er[1][3]=0.0119732;
	jer_param[1][4]=0.00308811;  jer_param_er[1][4]=9.34083E-5;
	jer_param[1][5]=0.000336138;  jer_param_er[1][5]=0.00152807;
	jer_param[1][6]=5.34059E-6;  jer_param_er[1][6]=5.37503E-6;
	jer_param[1][7]=-1.15979;  jer_param_er[1][7]=0.0761452;
	jer_param[1][8]=0.137285;  jer_param_er[1][8]=0.00491652;
	jer_param[1][9]=0.0020658;  jer_param_er[1][9]=3.46802E-5;
	jer_param[1][10]=-14.7515;  jer_param_er[1][10]=0.612517;
	jer_param[1][11]=4.80581;  jer_param_er[1][11]=0.195418;
	jer_param[1][12]=-0.127518;  jer_param_er[1][12]=0.0109911;
	jer_param[2][0]=0.100634;  jer_param_er[2][0]=0.00241696;
	jer_param[2][1]=-0.000253429;  jer_param_er[2][1]=8.43763E-6;
	jer_param[2][2]=-0.688927;  jer_param_er[2][2]=0.125514;
	jer_param[2][3]=0.655545;  jer_param_er[2][3]=0.0103474;
	jer_param[2][4]=0.00300742;  jer_param_er[2][4]=8.73698E-5;
	jer_param[2][5]=0.0105033;  jer_param_er[2][5]=0.00178783;
	jer_param[2][6]=-3.76786E-5;  jer_param_er[2][6]=5.83178E-6;
	jer_param[2][7]=-1.85956;  jer_param_er[2][7]=0.096736;
	jer_param[2][8]=0.139571;  jer_param_er[2][8]=0.00510727;
	jer_param[2][9]=0.00227237;  jer_param_er[2][9]=3.75247E-5;
	jer_param[2][10]=-19.2885;  jer_param_er[2][10]=1.05958;
	jer_param[2][11]=5.66883;  jer_param_er[2][11]=0.262952;
	jer_param[2][12]=-0.161499;  jer_param_er[2][12]=0.013508;
	jer_param[3][0]=0.0952545;  jer_param_er[3][0]=0.00263113;
	jer_param[3][1]=-0.000189574;  jer_param_er[3][1]=9.63887E-6;
	jer_param[3][2]=0.236287;  jer_param_er[3][2]=0.12846;
	jer_param[3][3]=0.731283;  jer_param_er[3][3]=0.0128143;
	jer_param[3][4]=0.00341218;  jer_param_er[3][4]=0.00011154;
	jer_param[3][5]=0.0258263;  jer_param_er[3][5]=0.0016674;
	jer_param[3][6]=-4.83423E-5;  jer_param_er[3][6]=5.80565E-6;
	jer_param[3][7]=-2.22975;  jer_param_er[3][7]=0.0737835;
	jer_param[3][8]=0.179229;  jer_param_er[3][8]=0.00515981;
	jer_param[3][9]=0.0024204;  jer_param_er[3][9]=4.113E-5;
	jer_param[3][10]=-10.4553;  jer_param_er[3][10]=0.561997;
	jer_param[3][11]=3.36745;  jer_param_er[3][11]=0.171045;
	jer_param[3][12]=-0.0618917;  jer_param_er[3][12]=0.00980421;
	jer_param[4][0]=0.115229;  jer_param_er[4][0]=0.0027871;
	jer_param[4][1]=-0.000157746;  jer_param_er[4][1]=9.95943E-6;
	jer_param[4][2]=-0.0905143;  jer_param_er[4][2]=0.134539;
	jer_param[4][3]=0.795373;  jer_param_er[4][3]=0.016243;
	jer_param[4][4]=0.00635758;  jer_param_er[4][4]=0.000146956;
	jer_param[4][5]=0.0345064;  jer_param_er[4][5]=0.00216488;
	jer_param[4][6]=-5.72569E-5;  jer_param_er[4][6]=7.07075E-6;
	jer_param[4][7]=-2.21206;  jer_param_er[4][7]=0.0958281;
	jer_param[4][8]=0.212192;  jer_param_er[4][8]=0.00725209;
	jer_param[4][9]=0.00349351;  jer_param_er[4][9]=6.06418E-5;
	jer_param[4][10]=-10.6021;  jer_param_er[4][10]=0.660483;
	jer_param[4][11]=3.23948;  jer_param_er[4][11]=0.191871;
	jer_param[4][12]=-0.0236045;  jer_param_er[4][12]=0.0106671;
	jer_param[5][0]=0.0900639;  jer_param_er[5][0]=0.00267652;
	jer_param[5][1]=-8.17704E-5;  jer_param_er[5][1]=7.7302E-6;
	jer_param[5][2]=-0.567021;  jer_param_er[5][2]=0.160166;
	jer_param[5][3]=0.706037;  jer_param_er[5][3]=0.0137402;
	jer_param[5][4]=0.0066424;  jer_param_er[5][4]=0.000110101;
	jer_param[5][5]=0.0161473;  jer_param_er[5][5]=0.00313166;
	jer_param[5][6]=-3.72664E-5;  jer_param_er[5][6]=8.95092E-6;
	jer_param[5][7]=-3.40086;  jer_param_er[5][7]=0.154954;
	jer_param[5][8]=0.163174;  jer_param_er[5][8]=0.00755835;
	jer_param[5][9]=0.00291679;  jer_param_er[5][9]=6.87186E-5;
	jer_param[5][10]=-10.8209;  jer_param_er[5][10]=1.39492;
	jer_param[5][11]=2.05925;  jer_param_er[5][11]=0.352187;
	jer_param[5][12]=0.184539;  jer_param_er[5][12]=0.0192903;
	jer_param[6][0]=0.110669;  jer_param_er[6][0]=0.00336557;
	jer_param[6][1]=-0.000302301;  jer_param_er[6][1]=1.15245E-5;
	jer_param[6][2]=-1.70381;  jer_param_er[6][2]=0.17128;
	jer_param[6][3]=0.84678;  jer_param_er[6][3]=0.0199186;
	jer_param[6][4]=0.00418531;  jer_param_er[6][4]=0.00022283;
	jer_param[6][5]=-0.00347175;  jer_param_er[6][5]=0.0027646;
	jer_param[6][6]=8.50977E-5;  jer_param_er[6][6]=1.10633E-5;
	jer_param[6][7]=-3.2283;  jer_param_er[6][7]=0.126136;
	jer_param[6][8]=0.0758951;  jer_param_er[6][8]=0.00708642;
	jer_param[6][9]=0.00415232;  jer_param_er[6][9]=8.11982E-5;
	jer_param[6][10]=-16.6687;  jer_param_er[6][10]=1.09606;
	jer_param[6][11]=5.18449;  jer_param_er[6][11]=0.288507;
	jer_param[6][12]=-0.190928;  jer_param_er[6][12]=0.016249;
	jer_param[7][0]=0.0658253;  jer_param_er[7][0]=0.00272349;
	jer_param[7][1]=-0.000196113;  jer_param_er[7][1]=1.14286E-5;
	jer_param[7][2]=-1.07207;  jer_param_er[7][2]=0.137955;
	jer_param[7][3]=0.889453;  jer_param_er[7][3]=0.015712;
	jer_param[7][4]=0.00152891;  jer_param_er[7][4]=0.000140114;
	jer_param[7][5]=-0.0287057;  jer_param_er[7][5]=0.00306348;
	jer_param[7][6]=0.000128678;  jer_param_er[7][6]=1.50299E-5;
	jer_param[7][7]=-2.67266;  jer_param_er[7][7]=0.129521;
	jer_param[7][8]=0.170338;  jer_param_er[7][8]=0.00800203;
	jer_param[7][9]=0.00119121;  jer_param_er[7][9]=7.89355E-5;
	jer_param[7][10]=-8.77829;  jer_param_er[7][10]=1.94651;
	jer_param[7][11]=2.5475;  jer_param_er[7][11]=0.581462;
	jer_param[7][12]=0.0673843;  jer_param_er[7][12]=0.0379263;
	jer_param[8][0]=0.058632;  jer_param_er[8][0]=0.00246133;
	jer_param[8][1]=-8.6882E-5;  jer_param_er[8][1]=7.27552E-6;
	jer_param[8][2]=-2.81222;  jer_param_er[8][2]=0.177198;
	jer_param[8][3]=0.87612;  jer_param_er[8][3]=0.0135608;
	jer_param[8][4]=0.000862571;  jer_param_er[8][4]=8.72723E-5;
	jer_param[8][5]=-0.0272999;  jer_param_er[8][5]=0.00336864;
	jer_param[8][6]=5.7274E-5;  jer_param_er[8][6]=1.19158E-5;
	jer_param[8][7]=-3.16373;  jer_param_er[8][7]=0.187887;
	jer_param[8][8]=0.218132;  jer_param_er[8][8]=0.00804389;
	jer_param[8][9]=-0.000115495;  jer_param_er[8][9]=5.26047E-5;
	jer_param[8][10]=6.59552;  jer_param_er[8][10]=4.14997;
	jer_param[8][11]=-3.91415;  jer_param_er[8][11]=1.16124;
	jer_param[8][12]=0.689584;  jer_param_er[8][12]=0.0728198;
	jer_param[9][0]=0.047587;  jer_param_er[9][0]=0.00405247;
	jer_param[9][1]=-3.67703E-5;  jer_param_er[9][1]=1.32232E-5;
	jer_param[9][2]=-2.82013;  jer_param_er[9][2]=0.277916;
	jer_param[9][3]=0.90659;  jer_param_er[9][3]=0.0207463;
	jer_param[9][4]=0.00138775;  jer_param_er[9][4]=0.000128634;
	jer_param[9][5]=-0.0150813;  jer_param_er[9][5]=0.00500665;
	jer_param[9][6]=4.34831E-5;  jer_param_er[9][6]=1.98595E-5;
	jer_param[9][7]=-4.24538;  jer_param_er[9][7]=0.258383;
	jer_param[9][8]=0.235095;  jer_param_er[9][8]=0.0112121;
	jer_param[9][9]=-0.000176939;  jer_param_er[9][9]=7.11021E-5;
	jer_param[9][10]=17.7377;  jer_param_er[9][10]=6.10388;
	jer_param[9][11]=-6.60975;  jer_param_er[9][11]=1.55793;
	jer_param[9][12]=0.757554;  jer_param_er[9][12]=0.0890909;
	jer_param[10][0]=0.0702999;  jer_param_er[10][0]=0.00524409;
	jer_param[10][1]=-0.0001322;  jer_param_er[10][1]=1.44554E-5;
	jer_param[10][2]=-4.64174;  jer_param_er[10][2]=0.433494;
	jer_param[10][3]=1.06671;  jer_param_er[10][3]=0.0316977;
	jer_param[10][4]=0.00105544;  jer_param_er[10][4]=0.000168547;
	jer_param[10][5]=-0.02437;  jer_param_er[10][5]=0.00587413;
	jer_param[10][6]=2.26747E-5;  jer_param_er[10][6]=1.97589E-5;
	jer_param[10][7]=-3.8046;  jer_param_er[10][7]=0.382791;
	jer_param[10][8]=0.268162;  jer_param_er[10][8]=0.0148949;
	jer_param[10][9]=-6.26511E-5;  jer_param_er[10][9]=8.37528E-5;
	jer_param[10][10]=10.2852;  jer_param_er[10][10]=7.77639;
	jer_param[10][11]=-4.34279;  jer_param_er[10][11]=1.80774;
	jer_param[10][12]=0.564928;  jer_param_er[10][12]=0.0947709;
	jer_param[11][0]=0.0276176;  jer_param_er[11][0]=0.00700865;
	jer_param[11][1]=-1.08182E-5;  jer_param_er[11][1]=1.77066E-5;
	jer_param[11][2]=-4.30341;  jer_param_er[11][2]=0.642783;
	jer_param[11][3]=1.15381;  jer_param_er[11][3]=0.0466366;
	jer_param[11][4]=0.000660203;  jer_param_er[11][4]=0.000219101;
	jer_param[11][5]=-0.0365288;  jer_param_er[11][5]=0.007012;
	jer_param[11][6]=5.34815E-5;  jer_param_er[11][6]=2.21572E-5;
	jer_param[11][7]=-4.83064;  jer_param_er[11][7]=0.460712;
	jer_param[11][8]=0.28416;  jer_param_er[11][8]=0.0195007;
	jer_param[11][9]=-4.88533E-7;  jer_param_er[11][9]=0.000101491;
	jer_param[11][10]=14.821;  jer_param_er[11][10]=8.92697;
	jer_param[11][11]=-5.32396;  jer_param_er[11][11]=2.01108;
	jer_param[11][12]=0.568491;  jer_param_er[11][12]=0.10153;
	jer_param[12][0]=-0.00456382;  jer_param_er[12][0]=0.00797029;
	jer_param[12][1]=5.93327E-5;  jer_param_er[12][1]=2.2512E-5;
	jer_param[12][2]=-6.05263;  jer_param_er[12][2]=0.641518;
	jer_param[12][3]=1.01168;  jer_param_er[12][3]=0.0580327;
	jer_param[12][4]=0.00124107;  jer_param_er[12][4]=0.000266261;
	jer_param[12][5]=-0.0672398;  jer_param_er[12][5]=0.0103374;
	jer_param[12][6]=0.000144969;  jer_param_er[12][6]=3.13729E-5;
	jer_param[12][7]=-7.44458;  jer_param_er[12][7]=0.643695;
	jer_param[12][8]=0.236166;  jer_param_er[12][8]=0.0233694;
	jer_param[12][9]=0.000272798;  jer_param_er[12][9]=0.000115376;
	jer_param[12][10]=15.9149;  jer_param_er[12][10]=16.3677;
	jer_param[12][11]=-5.06773;  jer_param_er[12][11]=3.26237;
	jer_param[12][12]=0.515062;  jer_param_er[12][12]=0.14614;
	jer_param[13][0]=-0.0364159;  jer_param_er[13][0]=0.0230243;
	jer_param[13][1]=0.000103308;  jer_param_er[13][1]=6.3942E-5;
	jer_param[13][2]=-2.04414;  jer_param_er[13][2]=1.61387;
	jer_param[13][3]=1.16241;  jer_param_er[13][3]=0.154151;
	jer_param[13][4]=0.00449765;  jer_param_er[13][4]=0.000680983;
	jer_param[13][5]=-0.0713694;  jer_param_er[13][5]=0.0128;
	jer_param[13][6]=6.17634E-5;  jer_param_er[13][6]=3.50159E-5;
	jer_param[13][7]=-6.69936;  jer_param_er[13][7]=0.883717;
	jer_param[13][8]=0.401773;  jer_param_er[13][8]=0.0401819;
	jer_param[13][9]=0.00048766;  jer_param_er[13][9]=0.000171866;
	jer_param[13][10]=-6.84164;  jer_param_er[13][10]=15.4788;
	jer_param[13][11]=1.14003;  jer_param_er[13][11]=2.85757;
	jer_param[13][12]=0.063226;  jer_param_er[13][12]=0.116727;
	jer_param[14][0]=0;  jer_param_er[14][0]=0;
	jer_param[14][1]=0;  jer_param_er[14][1]=0;
	jer_param[14][2]=0;  jer_param_er[14][2]=0;
	jer_param[14][3]=0;  jer_param_er[14][3]=0;
	jer_param[14][4]=0;  jer_param_er[14][4]=0;
	jer_param[14][5]=-0.127208;  jer_param_er[14][5]=0.0191958;
	jer_param[14][6]=5.65717E-5;  jer_param_er[14][6]=3.57586E-5;
	jer_param[14][7]=-4.42422;  jer_param_er[14][7]=1.87829;
	jer_param[14][8]=0.9097;  jer_param_er[14][8]=0.0515384;
	jer_param[14][9]=-0.000149165;  jer_param_er[14][9]=0.000158924;
	jer_param[14][10]=0;  jer_param_er[14][10]=0;
	jer_param[14][11]=0;  jer_param_er[14][11]=0;
	jer_param[14][12]=0;  jer_param_er[14][12]=0;

	if(i<15 && j<13)
	{
		if(code==0) value=jer_param[i][j];
		if(code==1) value=jer_param_er[i][j];
	}
	return value;
}




//___ This function returns upper limit of allowed jet energy fluctuation.
double JetSmearStudy::IntegralUpLimit(double jet_E, double eta_det, int stat_code, int syst_code) {
	double value=5.0; // maximum allowed: E(det)/E(true)-1<=5.0
	double parJER[5];
	parJER[0]=Jer_meanG(jet_E,eta_det,0,syst_code); // for now systcode is used as jer_stat_code
	parJER[1]=Jer_sigmaG(jet_E,eta_det,0,syst_code);
	parJER[2]=Jer_mpvL(jet_E,eta_det,0,syst_code);
	parJER[3]=Jer_sigmaL(jet_E,eta_det,0,syst_code);
	parJER[4]=Jer_normG(jet_E,eta_det,0,syst_code);

	double lim_L=5.0;
	double lim_G=5.0;
	//------ calculating limit for Landau part
	//--- This limit is roughly equivalent to 1.0E-6 significance level of JER Landau part if it is found in range [-1.0;10.0]
	if(fabs(1+parJER[2]-3.7*parJER[3])>0.0) lim_L=(3.7*parJER[3]-parJER[2])/(1+parJER[2]-3.7*parJER[3]);
	if(lim_L>5.0 || lim_L<parJER[3]) lim_L=5.0;
	//------ calculating limit for Gaussian part taking into account relative normalization of Gaussian
	if(parJER[4]>0.0 && parJER[4]<=1.0)
	{
		double sig_tmp=-log10(parJER[4]/(1+parJER[4]));
		double sig=4.2-sig_tmp; // 4.2 roughly corresponds to 4*Sigma significance level (~1.0E-5)
		if(sig<0.0) sig=0.0;
		double n_sigma=0.6807+1.035*sig-0.05574*sig*sig; // crude parametrization of N_sigma vs. significance
		if(fabs(1+parJER[0]-n_sigma*parJER[1])>0.0) lim_G=(n_sigma*parJER[1]-parJER[0])/(1+parJER[0]-n_sigma*parJER[1]);
		else // just to make a "smooth" transition between Landau and Gauss
		{
			if(parJER[0]<1.0 && parJER[0]>0.0) lim_G=1.0/parJER[0]-1.0;
			else lim_G=-1.0;
		}
		if(lim_G>5.0 || lim_G<parJER[1]) lim_G=5.0;
	}
	//------ deciding what limit should be used: Landau or Gauss
	if(parJER[4]>0.0) // gauss is non-zero
	{
		if(lim_L<lim_G) value=lim_G; // Landau lim smaller
		else value=lim_L;// Landau lim larger
	}
	else value=lim_L; // gauss is zero
	if(value>5.0) value=5.0;
	return value;
}




//_________________ returns mean of Gauss in JER: A0+A1*Pt+A2/E
double JetSmearStudy::Jer_meanG(double jet_E, double eta_det, int stat_code, int syst_code)
{
	int bin=WhatJEReta(eta_det);
	double p1=JERparam(0,bin,0);
	double p2=JERparam(0,bin,1);
	double p3=JERparam(0,bin,2);
	double val_err=0.0;
	double val=0.0;
	if(jet_E>0.0)
	{
		if(abs(syst_code)==1)
		{
			double t1=JERparam(1,bin,0);
			double t2=JERparam(1,bin,1)*jet_E;
			double t3=JERparam(1,bin,2)/jet_E;
			double t12=2.0*t1*t2*JERCorrCoeff(bin,0);
			double t13=2.0*t1*t3*JERCorrCoeff(bin,1);
			double t23=2.0*t2*t3*JERCorrCoeff(bin,2);
			double t=t1*t1+t2*t2+t3*t3+t12+t13+t23;
			if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
		}
		val=p1+p2*jet_E+p3/jet_E+val_err;
		if(val<-1.0 || val>5.0) val=0.0;
	}
	return val;
}
//_________________ returns sigma of Gauss in JER: sqrt(A3/E+A4)
double JetSmearStudy::Jer_sigmaG(double jet_E, double eta_det, int stat_code, int syst_code)
{
	int bin=WhatJEReta(eta_det);
	double p1=JERparam(0,bin,3);
	double p2=JERparam(0,bin,4);
	double val=0.0;
	double val_err=0.0;
	if(jet_E>0.0)
	{
		val=p1/jet_E+p2;
		if(val>0.0)
		{
			if(abs(syst_code)==2)
			{
				double t1=JERparam(1,bin,3)/jet_E;
				double t2=JERparam(1,bin,4);
				double t12=2.0*t1*t2*JERCorrCoeff(bin,3);
				double t=0.25*(t1*t1+t2*t2+t12)/val;
				if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
			}
			val=sqrt(val)+val_err;
		}
		if(val<0.0) val=0.0;
	}
	return val;
}
//_________________ returns mpv of Landau in JER: A5+A6*Pt+A7/E
double JetSmearStudy::Jer_mpvL(double jet_E, double eta_det, int stat_code, int syst_code)
{
	int bin=WhatJEReta(eta_det);
	double p1=JERparam(0,bin,5);
	double p2=JERparam(0,bin,6);
	double p3=JERparam(0,bin,7);
	double val_err=0.0;
	double val=0.0;
	if(jet_E>0.0)
	{
		if(abs(syst_code)==3)
		{
			double t1=JERparam(1,bin,5);
			double t2=JERparam(1,bin,6)*jet_E;
			double t3=JERparam(1,bin,7)/jet_E;
			double t12=2.0*t1*t2*JERCorrCoeff(bin,4);
			double t13=2.0*t1*t3*JERCorrCoeff(bin,5);
			double t23=2.0*t2*t3*JERCorrCoeff(bin,6);
			double t=t1*t1+t2*t2+t3*t3+t12+t13+t23;
			if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
		}
		val=p1+p2*jet_E+p3/jet_E+val_err;
		if(val<-1.0 || val>5.0) val=0.0;
	}
	return val;
}
//_________________ returns sigma of Landau in JER: sqrt(A8/E+A9)
double JetSmearStudy::Jer_sigmaL(double jet_E, double eta_det, int stat_code, int syst_code)
{
	int bin=WhatJEReta(eta_det);
	double p1=JERparam(0,bin,8);
	double p2=JERparam(0,bin,9);
	double val=0.0;
	double val_err=0.0;
	if(jet_E>0.0)
	{
		val=p1/jet_E+p2;
		if(val>0.0)
		{
			if(abs(syst_code)==4)
			{
				double t1=JERparam(1,bin,8)/jet_E;
				double t2=JERparam(1,bin,9);
				double t12=2.0*t1*t2*JERCorrCoeff(bin,7);
				double t=0.25*(t1*t1+t2*t2+t12)/val;
				if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
			}
			val=sqrt(val)+val_err;
		}
		if(val<0.0) val=0.0;
	}
	return val;
}

//_________________ returns normalization of Gauss: [A10+A11*sqrt(E)]/E+A12
double JetSmearStudy::Jer_normG(double jet_E, double eta_det, int stat_code, int syst_code)
{
	int bin=WhatJEReta(eta_det);
	double p1=JERparam(0,bin,10);
	double p2=JERparam(0,bin,11);
	double p3=JERparam(0,bin,12);
	double val=0.0;
	double val_err=0.0;
	if(jet_E>0.0)
	{
		if(abs(syst_code)==5)
		{
			double t1=JERparam(1,bin,10)/jet_E;
			double t2=JERparam(1,bin,11)/sqrt(jet_E);
			double t3=JERparam(1,bin,12);
			double t12=2.0*t1*t2*JERCorrCoeff(bin,8);
			double t13=2.0*t1*t3*JERCorrCoeff(bin,9);
			double t23=2.0*t2*t3*JERCorrCoeff(bin,10);
			double t=t1*t1+t2*t2+t3*t3+t12+t13+t23;
			if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
		}
		val=(p1+p2*sqrt(jet_E))/jet_E+p3+val_err;
	}
	if(val<0.0) val=0.0;
	return val;
}




JetSmearStudy::JetSmearStudy(const char* name, const char* title):
	TStnModule(name,title)
{
}

JetSmearStudy::~JetSmearStudy()
{
}

//________________ returns index for eta_bin in JER
int JetSmearStudy::WhatJEReta(double eta_det)
{
	int eta_bin;
	for(int i=0; i<14; i++)
	{
		if(fabs(eta_det)>=i*0.2 && fabs(eta_det)<(i+1)*0.2)
		{
			eta_bin=i;
			break;
		}
	}
	if(fabs(eta_det)>=2.8) eta_bin=14;
	return eta_bin;
}

//________________ returns correlation coefficients for JER params
double JetSmearStudy::JERCorrCoeff(int i, int j)
{
	double corr_param[15][11]; // [eta][param]
	corr_param[0][0]=-0.953359; // Gauss mean: c01
	corr_param[0][1]=-0.883424; // Gauss mean: c02
	corr_param[0][2]=0.783976; // Gauss mean: c12
	corr_param[0][3]=-0.769484; // Gauss sigma: c01
	corr_param[0][4]=-0.916005; // Landau mpv: c01
	corr_param[0][5]=-0.920646; // Landau mpv: c02
	corr_param[0][6]=0.774913; // Landau mpv: c12
	corr_param[0][7]=-0.798127; // Landau sigma: c01
	corr_param[0][8]=-0.979048; // norm: c01
	corr_param[0][9]=0.914083; // norm: c02
	corr_param[0][10]=-0.972488; // norm: c12
	corr_param[1][0]=-0.930448; // Gauss mean: c01
	corr_param[1][1]=-0.893613; // Gauss mean: c02
	corr_param[1][2]=0.743464; // Gauss mean: c12
	corr_param[1][3]=-0.855253; // Gauss sigma: c01
	corr_param[1][4]=-0.93506; // Landau mpv: c01
	corr_param[1][5]=-0.790079; // Landau mpv: c02
	corr_param[1][6]=0.633928; // Landau mpv: c12
	corr_param[1][7]=-0.791488; // Landau sigma: c01
	corr_param[1][8]=-0.979059; // norm: c01
	corr_param[1][9]=0.935759; // norm: c02
	corr_param[1][10]=-0.978165; // norm: c12
	corr_param[2][0]=-0.933467; // Gauss mean: c01
	corr_param[2][1]=-0.903664; // Gauss mean: c02
	corr_param[2][2]=0.769636; // Gauss mean: c12
	corr_param[2][3]=-0.834006; // Gauss sigma: c01
	corr_param[2][4]=-0.941456; // Landau mpv: c01
	corr_param[2][5]=-0.863316; // Landau mpv: c02
	corr_param[2][6]=0.730447; // Landau mpv: c12
	corr_param[2][7]=-0.795485; // Landau sigma: c01
	corr_param[2][8]=-0.973561; // norm: c01
	corr_param[2][9]=0.922559; // norm: c02
	corr_param[2][10]=-0.979269; // norm: c12
	corr_param[3][0]=-0.92146; // Gauss mean: c01
	corr_param[3][1]=-0.875047; // Gauss mean: c02
	corr_param[3][2]=0.720247; // Gauss mean: c12
	corr_param[3][3]=-0.820492; // Gauss sigma: c01
	corr_param[3][4]=-0.932934; // Landau mpv: c01
	corr_param[3][5]=-0.806167; // Landau mpv: c02
	corr_param[3][6]=0.666462; // Landau mpv: c12
	corr_param[3][7]=-0.752656; // Landau sigma: c01
	corr_param[3][8]=-0.97604; // norm: c01
	corr_param[3][9]=0.916931; // norm: c02
	corr_param[3][10]=-0.970166; // norm: c12
	corr_param[4][0]=-0.921373; // Gauss mean: c01
	corr_param[4][1]=-0.865997; // Gauss mean: c02
	corr_param[4][2]=0.713135; // Gauss mean: c12
	corr_param[4][3]=-0.812212; // Gauss sigma: c01
	corr_param[4][4]=-0.933552; // Landau mpv: c01
	corr_param[4][5]=-0.826206; // Landau mpv: c02
	corr_param[4][6]=0.696959; // Landau mpv: c12
	corr_param[4][7]=-0.731748; // Landau sigma: c01
	corr_param[4][8]=-0.972365; // norm: c01
	corr_param[4][9]=0.900603; // norm: c02
	corr_param[4][10]=-0.96338; // norm: c12
	corr_param[5][0]=-0.945056; // Gauss mean: c01
	corr_param[5][1]=-0.918674; // Gauss mean: c02
	corr_param[5][2]=0.808282; // Gauss mean: c12
	corr_param[5][3]=-0.811252; // Gauss sigma: c01
	corr_param[5][4]=-0.941345; // Landau mpv: c01
	corr_param[5][5]=-0.908828; // Landau mpv: c02
	corr_param[5][6]=0.805962; // Landau mpv: c12
	corr_param[5][7]=-0.737379; // Landau sigma: c01
	corr_param[5][8]=-0.974342; // norm: c01
	corr_param[5][9]=0.878683; // norm: c02
	corr_param[5][10]=-0.952568; // norm: c12
	corr_param[6][0]=-0.863692; // Gauss mean: c01
	corr_param[6][1]=-0.898828; // Gauss mean: c02
	corr_param[6][2]=0.683443; // Gauss mean: c12
	corr_param[6][3]=-0.864674; // Gauss sigma: c01
	corr_param[6][4]=-0.910948; // Landau mpv: c01
	corr_param[6][5]=-0.878636; // Landau mpv: c02
	corr_param[6][6]=0.709797; // Landau mpv: c12
	corr_param[6][7]=-0.790837; // Landau sigma: c01
	corr_param[6][8]=-0.976028; // norm: c01
	corr_param[6][9]=0.920123; // norm: c02
	corr_param[6][10]=-0.976345; // norm: c12
	corr_param[7][0]=-0.901329; // Gauss mean: c01
	corr_param[7][1]=-0.868814; // Gauss mean: c02
	corr_param[7][2]=0.664341; // Gauss mean: c12
	corr_param[7][3]=-0.879818; // Gauss sigma: c01
	corr_param[7][4]=-0.910909; // Landau mpv: c01
	corr_param[7][5]=-0.848758; // Landau mpv: c02
	corr_param[7][6]=0.654596; // Landau mpv: c12
	corr_param[7][7]=-0.854805; // Landau sigma: c01
	corr_param[7][8]=-0.98475; // norm: c01
	corr_param[7][9]=0.944593; // norm: c02
	corr_param[7][10]=-0.983111; // norm: c12
	corr_param[8][0]=-0.93546; // Gauss mean: c01
	corr_param[8][1]=-0.919387; // Gauss mean: c02
	corr_param[8][2]=0.775089; // Gauss mean: c12
	corr_param[8][3]=-0.885525; // Gauss sigma: c01
	corr_param[8][4]=-0.902932; // Landau mpv: c01
	corr_param[8][5]=-0.858334; // Landau mpv: c02
	corr_param[8][6]=0.664962; // Landau mpv: c12
	corr_param[8][7]=-0.83488; // Landau sigma: c01
	corr_param[8][8]=-0.989491; // norm: c01
	corr_param[8][9]=0.946333; // norm: c02
	corr_param[8][10]=-0.980041; // norm: c12
	corr_param[9][0]=-0.958622; // Gauss mean: c01
	corr_param[9][1]=-0.92783; // Gauss mean: c02
	corr_param[9][2]=0.813899; // Gauss mean: c12
	corr_param[9][3]=-0.918357; // Gauss sigma: c01
	corr_param[9][4]=-0.946212; // Landau mpv: c01
	corr_param[9][5]=-0.874383; // Landau mpv: c02
	corr_param[9][6]=0.732471; // Landau mpv: c12
	corr_param[9][7]=-0.890866; // Landau sigma: c01
	corr_param[9][8]=-0.993186; // norm: c01
	corr_param[9][9]=0.965285; // norm: c02
	corr_param[9][10]=-0.987298; // norm: c12
	corr_param[10][0]=-0.959116; // Gauss mean: c01
	corr_param[10][1]=-0.934235; // Gauss mean: c02
	corr_param[10][2]=0.822224; // Gauss mean: c12
	corr_param[10][3]=-0.929578; // Gauss sigma: c01
	corr_param[10][4]=-0.945125; // Landau mpv: c01
	corr_param[10][5]=-0.889286; // Landau mpv: c02
	corr_param[10][6]=0.741957; // Landau mpv: c12
	corr_param[10][7]=-0.900695; // Landau sigma: c01
	corr_param[10][8]=-0.992615; // norm: c01
	corr_param[10][9]=0.967849; // norm: c02
	corr_param[10][10]=-0.989424; // norm: c12
	corr_param[11][0]=-0.961204; // Gauss mean: c01
	corr_param[11][1]=-0.933701; // Gauss mean: c02
	corr_param[11][2]=0.822536; // Gauss mean: c12
	corr_param[11][3]=-0.938163; // Gauss sigma: c01
	corr_param[11][4]=-0.95; // Landau mpv: c01
	corr_param[11][5]=-0.875959; // Landau mpv: c02
	corr_param[11][6]=0.737455; // Landau mpv: c12
	corr_param[11][7]=-0.895806; // Landau sigma: c01
	corr_param[11][8]=-0.993612; // norm: c01
	corr_param[11][9]=0.968624; // norm: c02
	corr_param[11][10]=-0.988782; // norm: c12
	corr_param[12][0]=-0.961494; // Gauss mean: c01
	corr_param[12][1]=-0.888003; // Gauss mean: c02
	corr_param[12][2]=0.756614; // Gauss mean: c12
	corr_param[12][3]=-0.940026; // Gauss sigma: c01
	corr_param[12][4]=-0.971454; // Landau mpv: c01
	corr_param[12][5]=-0.907369; // Landau mpv: c02
	corr_param[12][6]=0.819359; // Landau mpv: c12
	corr_param[12][7]=-0.890781; // Landau sigma: c01
	corr_param[12][8]=-0.995606; // norm: c01
	corr_param[12][9]=0.976755; // norm: c02
	corr_param[12][10]=-0.991439; // norm: c12
	corr_param[13][0]=-0.979353; // Gauss mean: c01
	corr_param[13][1]=-0.901815; // Gauss mean: c02
	corr_param[13][2]=0.819153; // Gauss mean: c12
	corr_param[13][3]=-0.906452; // Gauss sigma: c01
	corr_param[13][4]=-0.977947; // Landau mpv: c01
	corr_param[13][5]=-0.902621; // Landau mpv: c02
	corr_param[13][6]=0.821221; // Landau mpv: c12
	corr_param[13][7]=-0.907155; // Landau sigma: c01
	corr_param[13][8]=-0.996171; // norm: c01
	corr_param[13][9]=0.98616; // norm: c02
	corr_param[13][10]=-0.996056; // norm: c12
	corr_param[14][0]=0; // Gauss mean: c01
	corr_param[14][1]=0; // Gauss mean: c02
	corr_param[14][2]=0; // Gauss mean: c12
	corr_param[14][3]=0; // Gauss sigma: c01
	corr_param[14][4]=-0.993129; // Landau mpv: c01
	corr_param[14][5]=-0.979168; // Landau mpv: c02
	corr_param[14][6]=0.957156; // Landau mpv: c12
	corr_param[14][7]=-0.895607; // Landau sigma: c01
	corr_param[14][8]=0; // norm: c01
	corr_param[14][9]=0; // norm: c02
	corr_param[14][10]=0; // norm: c12

	double value=0.0;
	if(i<15 && j<11) value=corr_param[i][j];

	return value;
}


