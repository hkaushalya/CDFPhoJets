{
//Before running you need to modify the lines with ////// comments in them. - Steve
/*
This is the script Steve used to make the final calibration 
and occupancy plots.
*/

gROOT->Reset();
gStyle.SetCanvasColor(10); 
gStyle.SetOptStat(1111111);
//gStyle.SetOptStat(0011000);
//gStyle.SetOptStat(0000000);
gStyle.SetHistLineWidth(2); gStyle.SetFuncWidth(3); 
gStyle.SetHistLineColor(1);
gStyle.SetMarkerStyle(22); gStyle.SetMarkerSize(1.0);
gStyle.SetMarkerColor(2); gStyle.SetTitleOffset(1.3,"Y");
gStyle.SetPadLeftMargin(0.15);
gStyle.SetPadTickX(1); gStyle.SetPadTickY(1); 
gStyle.SetOptFit(111); 
gStyle->SetPalette(1,0); 
gStyle.SetTitleTextColor(4); gStyle.SetTitleSize(0.05,"X");
gStyle.SetTitleSize(0.05,"Y"); gStyle.SetPadBottomMargin(0.15);
gROOT->ForceStyle();

                                              // database ler1 was teststand
  float lerdec = 219.87;  // December             database ler2
  //float ler2 = 211.0635;  // 195620 March 22nd
  float ler2 = 211.0635;  // 195624 March 23rd  database ler3
  float ler3 = 209.4476;  // 196367  April 10th
  float ler4 = 207.6802;  // 196664  April 17th
  float ler5 = 206.9211;  // 196989  April 26th   database ler4
  float ler6 = 207.4479;  // 197657  May 7th
  float ler7 = 207.1190;  // 198613  May 26th
  float ler8 = 203.1341;  // 199189  June 10th 2005   database ler5
  float ler9 = 204.5386;  // 200027  June 26th
  float ler10 = 209.1750;  // 200710  July 8th
  float ler11 = 202.4663;  // 201132  July 16th
  float ler12 = 199.9296;   // 201658  July 24th  ler6 made but never put in database
  float ler13 = 199.9296;   // 202367  August 7th
  float ler14 = 204.9342;   // 203348  August 28th
  float ler15 = 204.3491;  // 204550  Sep 20th
  float ler16 = 203.2912;   // 206723   Nov 7th JKim version
  //float ler16 = 204.47;    //206723  Ai Nagano version
  float ler17 = 202.5118;   //218658 June 19 2006
  float ler18 = 200.9812;   //222957 Oct 1 2006   ler7 in database
  float ler19 = 198.36;    //232522 Jan 15 2007 
  float ler20 = 206.555;   //233764 Feb 7  2007
  float ler21 = 208.1363;   //246125 Aug 1 2007
  float ler22 = 192.61;   //255202 Dec 17 2007
  float ler23 = 189.827;   //256904 Jan 29, 2008   ler8 in database
  float ler24 = 188.000;   //258365 
//  ***  Add ler25 and beyond global mean here.
  float ler25 = 190.415;	//new test run 258365 - sam
  float ler26 = 197.742;	//this is one of the old runs to refer 236806
  float ler27 = 212.025;	//run 266163 online integration test

	float ler28 = 191.364;	// p19:264101-264500
	float ler29 = 191.116;  // p19:265214_265463
	float ler30 = 190.275;  // p19:265464_265813
  	float ler31 = 190.347;	// p19:265814_266513

	float ler32 = 190.528;  //p20


  Float_t cp2lerin[48][54]; 
  float bgg,error;
  int err,id0,ipad,iwed,iside;
  std::string def_file("cp2traler-256904_LastDBEntry.txt"); // *** Change "default" run here
  FILE *fpp=fopen(def_file.c_str(),"r");
  assert(fpp != NULL && "default file did not open!");
  std::cout << "reading default file " << def_file  << std::endl;
  for (int i = 0; i < 2592; i++) {
    err = fscanf(fpp,"%d %f %f\n", &id0, &bgg, &error);
    if (err < 0) break;
    ipad  = (id0 & 0x0000003f);
    iwed  = ((id0 >> 8) & 0x0000001f);
    iside = ((id0 >> 6) & 0x00000001);
    cp2lerin[iwed+iside*24][ipad] = bgg*ler23/lerdec;        // *** change ler18 to default run
  }
  fclose(fpp);
  

  Float_t cp2lerin1[48][54];
  float bgg1,error1;
  int err1,id01,ipad1,iwed1,iside1;
  //int err1;
  std::string new_file("P19Gain_AB264101_264500.txt");         //*** change new run here
  FILE *fpp1=fopen(new_file.c_str(),"r");         //*** change new run here
  
  std::cout << "reading new file " << new_file << std::endl;
  assert(fpp1 != NULL && "new file did not open!");
  for (int i=0;i<2592;i++) {
    err1 = fscanf(fpp1,"%d %f %f\n",&id0,&bgg1,&error1);
    if (err1<0) break;
    ipad  = (id0 & 0x0000003f);
    iwed  = ((id0 >> 8) & 0x0000001f);
    iside = ((id0 >> 6) & 0x00000001);
    cp2lerin1[iwed+iside*24][ipad] = bgg1*ler32/lerdec;    //*** change ler23 to new run
	 																			// we are looking a shift wrt to the last database entry - sam
  }
  fclose(fpp1);
  

  TH1F *ratio;
  ratio = new TH1F("","",100,0,3);
  ratio->SetXTitle("New / Old");
  ratio->SetTitle("Ratio of New to Old LERs");

  TH1F *oldler;
  oldler = new TH1F("","",100,0.0,3.0);
  oldler->SetXTitle("Old LERs");
  oldler->SetLineColor(2);

  TH1F *newler;
  newler = new TH1F("","",100,0.0,3.0);
  newler->SetXTitle("New LERs");

//add a 2D plot of wedge vs pad # weighted by this ratio
  TH2F* cp2Weight = new TH2F("cp2Weight","",48,-0.5,47.5,54,-0.5,53.5);
  cp2Weight->SetTitle("CPR2 Pad Hits Weighted by the ratio");
  cp2Weight->SetXTitle("Wedge Number");
  cp2Weight->SetYTitle("CPR2 Pad Number");
//  cp2Weight->SetMaximum(1.3); 
//  cp2Weight->SetMinimum(0.7); 
  cp2Weight->SetMaximum(2.); 
  cp2Weight->SetMinimum(0.); 
  cp2Weight->SetStats(kFALSE);
  
//add a 1D plot of wedge weighted by this ratio
  TH1F* cp2wedwei = new TH1F("cp2wedwei","",48,-0.5,47.5);
  cp2wedwei->SetTitle("CPR2 Wedge Hits Weighted by the LER Ratio");
  cp2wedwei->SetXTitle("Wedge Number");
  cp2wedwei->SetMaximum(1.5); 
  cp2wedwei->SetMinimum(0.5); 
  cp2wedwei->SetStats(kFALSE);

//add a 1D plot of pad weighted by this ratio
  TH1F* cp2padwei = new TH1F("cp2padwei","",54,-0.5,53.5);
  cp2padwei->SetTitle("CPR2 Pad Hits Weighted by the LER Ratio");
  cp2padwei->SetXTitle("Pad Number");
  cp2padwei->SetMaximum(1.5); 
  cp2padwei->SetMinimum(0.5); 
  cp2padwei->SetStats(kFALSE);

//add a 1D plot of pad weighted by this ratio for a special wedge cut
  TH1F* cp2padwei2 = new TH1F("cp2padwei2","",54,-0.5,53.5);
  cp2padwei2->SetTitle("CPR2 Pad Hits Weighted by the LER Ratio");
  cp2padwei2->SetXTitle("Pad Number");
  cp2padwei2->SetMaximum(2.0); 
  cp2padwei2->SetMinimum(0.); 
  cp2padwei2->SetStats(kFALSE);


  Float_t rat;
	for (int iside=0; iside<2; iside++) {
		for (int iwed=0; iwed<24; iwed++) {
			for (int ipad=0; ipad<54; ipad++) {
         	rat = cp2lerin1[iwed + iside * 24][ipad] / cp2lerin[iwed + iside * 24][ipad];
          	if (iwed + iside * 24 ==  0 && ipad == 2) continue;
	      	if (iwed + iside * 24 == 31 && ipad == 5) continue;
	      	if (iwed + iside * 24 == 39 && ipad == 3) continue;
          	ratio->Fill(rat,1);
          	oldler->Fill(cp2lerin[iwed+iside*24][ipad],1);
          	newler->Fill(cp2lerin1[iwed+iside*24][ipad],1);
	      	//cp2pad[iwed+iside*24][ipad]->Fill(rat,1);
  	      	cp2Weight->Fill(iwed+iside*24,ipad,rat);
  	      	//cp2Weight->Fill(iwed+iside*24,ipad,cp2lerin1[iwed+iside*24][ipad]);
	      	//cp2wedwei->Fill(iwed+iside*24,cp2lerin1[iwed+iside*24][ipad]);
	      	//cp2padwei->Fill(ipad,cp2lerin1[iwed+iside*24][ipad]);
	      	//if ((iwed+iside*24)==29) cp2padwei2->Fill(ipad,cp2lerin1[iwed+iside*24][ipad]);
	      	cp2wedwei->Fill(iwed+iside*24,rat/54.);
	      	cp2padwei->Fill(ipad,rat/48.);
	      	if ((iwed+iside*24)==29) cp2padwei2->Fill(ipad,rat);
			}
		}
	}


TCanvas* c2 = new TCanvas("c2","",0,0,900,700);
c2.Divide(2,2); c2.cd(1);
oldler->Fit("gaus");
oldler->Draw();
c2.cd(2);
cp2Weight->Draw();   //is this is correct hist for this place? compare with cvs version
newler->Fit("gaus");
newler->Draw();
c2.cd(3);
ratio->Fit("gaus"); 
ratio->Draw();
oldler->Draw("same");  
c2.cd(4);
cp2wedwei->Draw();
TLine * ln0 = new TLine (-0.5, 1., 47.5, 1.0); ln0->Draw();
//c2.Print("ratio-222957-256904-fits.eps");   /////// change output plot name here
std::string::size_type dotpos = new_file.find(".",0);
std::string c1name = "ratio-fits_"+ new_file.substr(0,dotpos) + ".gif";
assert(c1name.length()>0 && "canvas name not specified!");
//c2.Print("ratio-266163-fits.gif");   /////// change output plot name here
c2.Print(c1name.c_str());   /////// change output plot name here

TCanvas* c3 = new TCanvas("c3","",0,0,900,700);
c3.Divide(2,2); c3.cd(1);
cp2padwei->Draw();
TLine * ln0 = new TLine (-0.5, 1., 53.5, 1.0); ln0->Draw();
c3.cd(2);
cp2Weight->Draw("colz");
c3.cd(3);
cp2padwei2->SetTitle("Wedge 5East");
cp2padwei2->Draw();
//c3.Print("ratio-222957-256904-pads.eps");    ///////// change output plot name here
std::string c3name = "ratio-pads_"+ new_file.substr(0,dotpos) + ".gif";
//c3.Print("ratio-266163-pads.gif");    ///////// change output plot name here
c3.Print(c3name.c_str());    ///////// change output plot name here

}
