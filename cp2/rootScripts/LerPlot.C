#include <vector>
#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <string>
void LerPlot()
{

TH1F* ler1 = new TH1F("LER1","CP2 Global AVG LERs Vs Run Peiod;Run Period;Global Avg. LER", 30,0,30);
TH1F* histler2= new TH1F("LER2","CP2 Global AVG LERs/Test Stand  LER  Vs Run Peiod;Run Period; #frac{Avg. LER}{LER Teststand(==lerdec)}", 30,0,30);

  float lerdec = 219.87;   // December             database ler2
  float ler2 = 211.0635;   //P2 195624 March 23rd    database ler3
  float ler3 = 209.4476;   //P2 196367  April 10th
  float ler4 = 207.6802;   //P2 196664  April 17th
  float ler5 = 206.9211;   //P2 196989  April 26th   database ler4
  float ler6 = 207.4479;   //P2 197657  May 7th
  float ler7 = 207.1190;   //P3 198613  May 26th
  float ler8 = 203.1341;   //P3 199189  June 10th    database ler5
  float ler9 = 204.5386;   //P3 200027  June 26th
  float ler10 = 209.1750;  //P3 200710  July 8th
  float ler11 = 202.4663;  //P3 201132  July 16th
  float ler12 = 199.9296;  //P4 201658  July 24th   database ler6 never put in
  float ler14 = 204.9342;  //P4- 203348  August 28th
  float ler15 = 204.3491;  //P5- 204550  Sep 20th
  float ler16 = 203.2912;  //P5- 206723   Nov 7th JKim version
  //float ler16 = 204.47;  //P5- 206723  Ai Nagano version
  float ler17 = 202.5118;  //P8- 218658 June 19 2006
  float ler18 = 200.9812;  //P10- 222957 Oct 1 2006  database ler7
  float ler19 = 198.36;    //P10- 232522 Jan 15 2007 
  float ler20 = 206.555;   //P11- 233764 Feb 7  2007
  float ler21 = 208.1363;  //P13- 246125 Aug 1 2007
  float ler22 = 192.61;    //P15- 255202 Dec 17 2007
  float ler23 = 189.827;   //P16 - 256904 Jan 29, 2008  database ler8
  float ler24 = 192.254;   //P21:267719 - 271000  database ler 9

  float ler25 = 191.447;   //P22:271071_272214 database ler 10
  float ler26 = 192.359; 	//P23:271071 - 272214 database ler 11
  float ler27 = 190.771;	//Pp24:274123-275848 database ler 12
  float ler28 = 189.888;	//P25:275873-277511  database ler 13
  float ler29 = 190.787; 	//P26: 282976-284843 (2009 post sutdown) database ler 14
  
  std::vector < std::string > labels;
  labels.push_back("*");
  labels.push_back("2");
  labels.push_back("2");
  labels.push_back("2");
  labels.push_back("2");
  labels.push_back("2");
  labels.push_back("3");
  labels.push_back("3");
  labels.push_back("3");
  labels.push_back("3");
  labels.push_back("3");
  labels.push_back("4");
  labels.push_back("4");
  labels.push_back("5");
  labels.push_back("5");
  labels.push_back("8");
  labels.push_back("10");
  labels.push_back("10");
  labels.push_back("11");
  labels.push_back("13");
  labels.push_back("15");
  labels.push_back("16");
  labels.push_back("21");
  labels.push_back("22");
  labels.push_back("23");
  labels.push_back("24");
  labels.push_back("25");
  labels.push_back("26");
  //assert (labels.size() == ler1->GetNbinsX() && "labels and N bins do not match");
  for (int i=1; i<=labels.size(); ++i)
  {
	  ler1->GetXaxis()->SetBinLabel(i,labels.at(i-1).c_str());
	  histler2->GetXaxis()->SetBinLabel(i,labels.at(i-1).c_str());
  }
 

  ler1->SetBinContent(1,lerdec);
  ler1->SetBinContent(2,ler2);
  ler1->SetBinContent(3,ler3);
  ler1->SetBinContent(4,ler4);
  ler1->SetBinContent(5,ler5);
  ler1->SetBinContent(6,ler6);
  ler1->SetBinContent(7,ler7);
  ler1->SetBinContent(8,ler8);
  ler1->SetBinContent(9,ler9);
  ler1->SetBinContent(10,ler10);
  ler1->SetBinContent(11,ler11);
  ler1->SetBinContent(12,ler12);
  ler1->SetBinContent(13,ler14);
  ler1->SetBinContent(14,ler15);
  ler1->SetBinContent(15,ler16);
  ler1->SetBinContent(16,ler17);
  ler1->SetBinContent(17,ler18);
  ler1->SetBinContent(18,ler19);
  ler1->SetBinContent(19,ler20);
  ler1->SetBinContent(20,ler21);
  ler1->SetBinContent(21,ler22);
  ler1->SetBinContent(22,ler23);
  ler1->SetBinContent(23,ler24);
  ler1->SetBinContent(24,ler25);
  ler1->SetBinContent(25,ler26);
  ler1->SetBinContent(26,ler27);
  ler1->SetBinContent(27,ler28);
  ler1->SetBinContent(28,ler29);
  

  histler2->SetBinContent(1,lerdec/lerdec);
  histler2->SetBinContent(2,ler2/lerdec);
  histler2->SetBinContent(3,ler3/lerdec);
  histler2->SetBinContent(4,ler4/lerdec);
  histler2->SetBinContent(5,ler5/lerdec);
  histler2->SetBinContent(6,ler6/lerdec);
  histler2->SetBinContent(7,ler7/lerdec);
  histler2->SetBinContent(8,ler8/lerdec);
  histler2->SetBinContent(9,ler9/lerdec);
  histler2->SetBinContent(10,ler10/lerdec);
  histler2->SetBinContent(11,ler11/lerdec);
  histler2->SetBinContent(12,ler12/lerdec);
  histler2->SetBinContent(13,ler14/lerdec);
  histler2->SetBinContent(14,ler15/lerdec);
  histler2->SetBinContent(15,ler16/lerdec);
  histler2->SetBinContent(16,ler17/lerdec);
  histler2->SetBinContent(17,ler18/lerdec);
  histler2->SetBinContent(18,ler19/lerdec);
  histler2->SetBinContent(19,ler20/lerdec);
  histler2->SetBinContent(20,ler21/lerdec);
  histler2->SetBinContent(21,ler22/lerdec);
  histler2->SetBinContent(22,ler23/lerdec);
  histler2->SetBinContent(23,ler24/lerdec);
  histler2->SetBinContent(24,ler25/lerdec);
  histler2->SetBinContent(25,ler26/lerdec);
  histler2->SetBinContent(26,ler27/lerdec);
  histler2->SetBinContent(27,ler28/lerdec);
  histler2->SetBinContent(28,ler29/lerdec);
 
  ler1->SetMarkerStyle(20);
  ler1->SetMarkerColor(kBlue);
  ler1->SetMarkerSize(1);
  histler2->SetMarkerStyle(20);
  histler2->SetMarkerColor(kBlue);
  histler2->SetMarkerSize(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  ler1->GetXaxis()->CenterTitle(1);
  ler1->GetXaxis()->LabelsOption("h");
  ler1->GetYaxis()->CenterTitle(1);
  ler1->Fit("pol1");
  new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTickx();
  gPad->SetTicky();
  ler1->Draw("P");
  new TCanvas();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTickx();
  gPad->SetTicky();
  histler2->GetXaxis()->CenterTitle(1);
  histler2->GetXaxis()->LabelsOption("h");
  histler2->GetYaxis()->CenterTitle(1);
  histler2->Draw("P");

}
