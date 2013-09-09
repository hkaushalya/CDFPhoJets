#include "samantha/utils/JERFunction.hh"
#include<iostream>


//this is a helper class for JetFilterModuleV2
//It handles all the JER function methods(requests)
//create a CLASSof JER and access your info via 
//set of common methods
//


//-----------------------------------------------------------------------------
int JERFunction::EtaBin(const double fEtaDet) const 
{
	int iEtaBin;
	
	if (fabs(fEtaDet)>=2.8) iEtaBin=14;
	else 
	{
		for (int i=0; i<14; i++)
		{
			if (fabs(fEtaDet)>=i*0.2 && fabs(fEtaDet)<(i+1)*0.2)
			{
				iEtaBin=i;
				break;
			}
		}
	}
	return iEtaBin;

}

//-----------------------------------------------------------------------------
JERFunction::JERFunction()
{
}


//-----------------------------------------------------------------------------
JERFunction::JERFunction(JERCLASS_t jerclass, const double fJetE, const double fDetEta, const int syst)
{
	//first pick the custom JER function based on the class choosen
	JERCLASS = jerclass;
	int iEtaBin = EtaBin(fDetEta);
	switch (jerclass)
	{
		case DEFAULT:
			{
				//FuncTf1 = new TF1("jer_tf1",c1_JER,0,500,5);

				break;
			}
		case SYST:
			{
				//FuncTf1 = new TF1("jer_tf1",c2_JER,0,500,5);
				break;
			}
		default:
			{
				std::cout << __FILE__ << "::" << __LINE__ << "::JER CLASS NOT DEFINED!" << std::endl;
				exit(1);
			}
	}

	//assign initial parameters to the function
	//SetParameters(fJetE, iEtaBin);

	
	
}


//-----------------------------------------------------------------------------
double JERFunction::RandomPoint() const
{
	return FuncTf1->GetRandom();
}

//-----------------------------------------------------------------------------
double JERFunction::Integral(const double xlow, const double xmax) const
{
	return FuncTf1->Integral(xlow,xmax);
}


//-----------------------------------------------------------------------------
double JERFunction::JERParameter(const int iEtaBin, const int iPar, const int iSystCode)
{
	switch (JERCLASS)
	{
		case DEFAULT:
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

				if (iEtaBin<15 && iPar<13)
				{
					if (iSystCode == 0) value=jer_param[iEtaBin][iPar];
					if (iSystCode == 1) value=jer_param_er[iEtaBin][iPar];
				}
				return value;
			} //class 1=DEFAULT
	}
}


